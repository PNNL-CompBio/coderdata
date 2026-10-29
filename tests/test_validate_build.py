"""
Tests for coderbuild/utils/validate_build.py.

Each test writes a tiny build with one planted problem modelled on what
shipped in CoderData 2.4, and checks that the validator reports it (and that
a clean build passes).
"""
import gzip
import json
import os
import sys

import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'coderbuild', 'utils'))

import validate_build as vb  # noqa: E402

EXP_COLS = ['source', 'improve_sample_id', 'improve_drug_id', 'study', 'time', 'time_unit',
            'dose_response_metric', 'dose_response_value']
DRUG_COLS = ['improve_drug_id', 'chem_name', 'pubchem_id', 'canSMILES', 'InChIKey', 'formula', 'weight']
SAMPLE_COLS = ['other_id', 'improve_sample_id', 'other_names', 'common_name', 'cancer_type',
               'model_type', 'other_id_source', 'species']


def _tsv(path, rows, cols, sep='\t'):
    pd.DataFrame(rows, columns=cols).to_csv(path, sep=sep, index=False)


def _samples(d, ds, ids, model='patient derived organoid', cancer='Colon cancer'):
    _tsv(d / f'{ds}_samples.csv', [[f'X{i}', i, f'n{i}', f'c{i}', cancer, model, 'src', 'Homo sapiens']
                                    for i in ids], SAMPLE_COLS, sep=',')


def _drugs(d, ds, rows):
    _tsv(d / f'{ds}_drugs.tsv', rows, DRUG_COLS)


def _exp(d, ds, rows):
    _tsv(d / f'{ds}_experiments.tsv', rows, EXP_COLS)


def _mappings(d, samples=None, drugs=None):
    samples = samples or [{'stable_id': '1', 'quadruples': [
        {'other_id': 'X1', 'other_id_source': 'src', 'model_type': 'patient derived organoid', 'other_names': 'n1'}]}]
    drugs = drugs or [{'stable_id': 'SMI_1', 'canSMILES': 'CCO'}]
    with gzip.open(d / 'improve_sample_mapping.json.gz', 'wt') as f:
        json.dump({'metadata': {'builds': []}, 'samples': samples}, f)
    with gzip.open(d / 'improve_drug_mapping.json.gz', 'wt') as f:
        json.dump({'metadata': {'builds': []}, 'drugs': drugs}, f)


def _clean_build(d):
    _samples(d, 'liver', [1, 2, 3])
    _drugs(d, 'liver', [['SMI_1', 'ethanol', '702', 'CCO', 'K', 'F', '1'],
                        ['SMI_2', 'taxol', '36314', 'CCC', 'K2', 'F', '1']])
    rows = []
    for s, a, b in ((1, '0.5', '0.9'), (2, '0.6', '0.2'), (3, '0.7', '0.4')):
        rows.append(['syn', s, 'SMI_1', 'L', '72', 'hrs', 'fit_auc', a])
        rows.append(['syn', s, 'SMI_2', 'L', '72', 'hrs', 'fit_auc', b])
    _exp(d, 'liver', rows)
    _tsv(d / 'liver_transcriptomics.csv', [[1, 7105, 1.0], [2, 7105, 2.0]],
         ['improve_sample_id', 'entrez_id', 'transcriptomics'], sep=',')
    _mappings(d)


def _checks(report, kind='errors'):
    return {c for c, _, _ in getattr(report, kind)}


def test_clean_build_passes(tmp_path):
    _clean_build(tmp_path)
    report, summary = vb.validate(str(tmp_path))
    assert report.errors == []
    assert summary['datasets']['liver']['pairs'] == 6
    assert vb.main(['--dir', str(tmp_path)]) == 0


def test_duplicate_results_fail(tmp_path):
    # issue #512: two curves for one sample/drug
    _clean_build(tmp_path)
    _exp(tmp_path, 'liver', [['syn', 1, 'SMI_1', 'L', '72', 'hrs', 'fit_auc', '0.82'],
                             ['syn', 1, 'SMI_1', 'L', '72', 'hrs', 'fit_auc', '0.78']])
    report, _ = vb.validate(str(tmp_path))
    assert 'duplicate-result' in _checks(report)
    assert vb.main(['--dir', str(tmp_path)]) == 1


def test_decimal_sample_ids_are_orphans(tmp_path):
    # novartis 2.4: experiments wrote sample ids as '3524.0'
    _clean_build(tmp_path)
    _exp(tmp_path, 'liver', [['syn', '1.0', 'SMI_1', 'L', '72', 'hrs', 'fit_auc', '0.5']])
    report, _ = vb.validate(str(tmp_path))
    msg = [m for c, _, m in report.errors if c == 'orphan-sample-id'][0]
    assert 'written as decimals' in msg


def test_orphan_drug_ids_fail(tmp_path):
    _clean_build(tmp_path)
    _exp(tmp_path, 'liver', [['syn', 1, 'SMI_9', 'L', '72', 'hrs', 'fit_auc', '0.5']])
    report, _ = vb.validate(str(tmp_path))
    assert 'orphan-drug-id' in _checks(report)


def test_copied_responses_warn(tmp_path):
    _clean_build(tmp_path)
    rows = [['syn', s, d, 'L', '72', 'hrs', 'fit_auc', v]
            for d in ('SMI_1', 'SMI_2') for s, v in ((1, '0.5'), (2, '0.6'), (3, '0.7'))]
    _exp(tmp_path, 'liver', rows)
    report, _ = vb.validate(str(tmp_path))
    assert 'copied-responses' in _checks(report, 'warnings')
    assert report.errors == []


def test_sample_ids_shared_by_two_datasets_fail(tmp_path):
    _clean_build(tmp_path)
    _samples(tmp_path, 'bladder', [3, 4])
    report, _ = vb.validate(str(tmp_path))
    assert 'sample-id-overlap' in _checks(report)


def test_cell_lines_may_share_ids(tmp_path):
    for ds in ('ccle', 'prism'):
        _samples(tmp_path, ds, [1, 2], model='cell line', cancer='Lung cancer')
    _mappings(tmp_path)
    report, _ = vb.validate(str(tmp_path))
    assert 'sample-id-overlap' not in _checks(report)


def test_conflicting_model_types_fail(tmp_path):
    _clean_build(tmp_path)
    rows = [['X1', 1, 'n1', 'c1', 'Colon cancer', 'tumor', 'src', 'Homo sapiens'],
            ['X1b', 1, 'n1', 'c1', 'Colon cancer', 'patient derived organoid', 'src2', 'Homo sapiens'],
            ['X2', 2, 'n2', 'c2', 'Colon cancer', 'tumor', 'src', 'Homo sapiens'],
            ['X3', 3, 'n3', 'c3', 'Colon cancer', 'tumor', 'src', 'Homo sapiens']]
    _tsv(tmp_path / 'liver_samples.csv', rows, SAMPLE_COLS, sep=',')
    report, _ = vb.validate(str(tmp_path))
    assert 'sample-model-type' in _checks(report)


def test_drug_with_two_structures_or_split_cid_fails(tmp_path):
    _clean_build(tmp_path)
    _drugs(tmp_path, 'liver', [['SMI_1', 'ethanol', '702', 'CCO', 'K', 'F', '1'],
                               ['SMI_1', 'ethyl alcohol', '702', 'OCC', 'K', 'F', '1'],
                               ['SMI_2', 'taxol', '36314', 'CCC', 'K2', 'F', '1'],
                               ['SMI_3', 'paclitaxel', '36314.0', 'CCCC', 'K2', 'F', '1']])
    report, _ = vb.validate(str(tmp_path))
    assert 'drug-structure' in _checks(report)
    assert 'drug-split' in _checks(report)  # one chemical under two ids (issue #512)


def test_blank_ids_fail(tmp_path):
    _clean_build(tmp_path)
    _exp(tmp_path, 'liver', [['syn', 1, 'SMI_1', 'L', '72', 'hrs', 'fit_auc', '0.5'],
                             ['syn', '', 'SMI_1', 'L', '72', 'hrs', 'fit_auc', '0.6'],
                             ['syn', 2, '', 'L', '72', 'hrs', 'fit_auc', '0.7']])
    report, _ = vb.validate(str(tmp_path))
    blank = [m for c, _, m in report.errors if c == 'missing-id']
    assert any('blank sample id' in m for m in blank) and any('blank drug id' in m for m in blank)
    assert 'orphan-sample-id' not in _checks(report)  # blanks are not reported twice


def test_decimal_ids_in_the_samples_file_fail(tmp_path):
    # ids written as decimals everywhere still match each other, so only the
    # format check catches them
    _clean_build(tmp_path)
    _tsv(tmp_path / 'liver_samples.csv', [[f'X{i}', f'{i}.0', f'n{i}', f'c{i}', 'Colon cancer',
                                           'patient derived organoid', 'src', 'Homo sapiens'] for i in (1, 2, 3)],
         SAMPLE_COLS, sep=',')
    rows = [['syn', f'{s}.0', 'SMI_1', 'L', '72', 'hrs', 'fit_auc', v] for s, v in ((1, '0.5'), (2, '0.6'), (3, '0.7'))]
    _exp(tmp_path, 'liver', rows)
    report, _ = vb.validate(str(tmp_path))
    assert 'decimal-id' in _checks(report)


def test_mapping_identity_with_two_ids_fails(tmp_path):
    _clean_build(tmp_path)
    _mappings(tmp_path, drugs=[{'stable_id': 'SMI_1', 'canSMILES': 'CCO'},
                               {'stable_id': 'SMI_2', 'canSMILES': 'CCO'}])
    report, _ = vb.validate(str(tmp_path))
    assert 'mapping-file' in _checks(report)


def test_categorical_metrics_are_not_range_checked(tmp_path):
    _clean_build(tmp_path)
    _exp(tmp_path, 'liver', [['syn', 1, 'SMI_1', 'L', '72', 'hrs', 'mrecist', 'PD'],
                             ['syn', 2, 'SMI_1', 'L', '72', 'hrs', 'fit_auc', '9.0']])
    report, _ = vb.validate(str(tmp_path))
    assert report.errors == []
    assert any('outside' in m for c, _, m in report.warnings if c == 'value-ranges')


def test_release_changes_and_large_drops_warn(tmp_path):
    # hcmi 2.2 -> 2.4 lost its tumor and cell-line samples
    _clean_build(tmp_path)
    previous = {'datasets': {'liver': {'samples': 10, 'drugs': 2, 'pairs': 6,
                                       'model_types': {'tumor': 5, 'patient derived organoid': 5},
                                       'omics': {'transcriptomics': 2, 'mutations': 2}},
                             'hcmi': {'samples': 5}}}
    report, _ = vb.validate(str(tmp_path), previous=previous)
    messages = [m for c, _, m in report.warnings if c == 'release-change']
    drops = [m for c, _, m in report.warnings if c == 'release-drop']
    assert any(m.startswith('samples: 10 -> 3') for m in drops)  # -70%: flagged as a drop
    assert report.errors == []
    assert any('removed' in m and 'tumor' in m for m in messages)
    assert any('mutations' in m for m in messages)
    assert ('dataset-missing', 'hcmi') in {(c, ds) for c, ds, _ in report.warnings}


def test_moderate_drop_only_warns(tmp_path):
    # nci60 2.4 rebuild lost 98% of its drugs; a 30% change is only a warning
    _clean_build(tmp_path)
    _, summary = vb.validate(str(tmp_path))
    liver = summary['datasets']['liver']
    previous = {'datasets': {'liver': dict(liver, samples=round(liver['samples'] / 0.7))}}
    report, _ = vb.validate(str(tmp_path), previous=previous)
    assert 'release-drop' not in _checks(report, 'warnings')
    assert any(m.startswith('samples:') for c, _, m in report.warnings if c == 'release-change')


def test_datasets_left_out_on_purpose_are_not_missing(tmp_path):
    _clean_build(tmp_path)
    _, summary = vb.validate(str(tmp_path))
    previous = {'datasets': dict(summary['datasets'], hcmi={'samples': 5})}
    report, _ = vb.validate(str(tmp_path), ['liver'], previous)
    assert 'dataset-missing' not in _checks(report, 'warnings')


def test_summary_round_trip(tmp_path):
    _clean_build(tmp_path)
    out = tmp_path / 'summary.json'
    assert vb.main(['--dir', str(tmp_path), '--write-summary', str(out)]) == 0
    assert vb.main(['--dir', str(tmp_path), '--previous-summary', str(out)]) == 0
    assert json.loads(out.read_text())['datasets']['liver']['pairs'] == 6
    # releases publish the summary gzipped
    gz = tmp_path / 'validation_summary.json.gz'
    with gzip.open(gz, 'wt') as f:
        f.write(out.read_text())
    assert vb.load_summary(str(gz)) == json.loads(out.read_text())
    assert vb.main(['--dir', str(tmp_path), '--previous-summary', str(gz)]) == 0


def test_copied_samples_warn(tmp_path):
    # novartis: model NIBRX-2043 was attached to sample ids 3828 and 3829
    _clean_build(tmp_path)
    rows = [['syn', s, d, 'L', '72', 'hrs', 'fit_auc', v]
            for s in (1, 2) for d, v in (('SMI_1', '0.5'), ('SMI_2', '0.6'), ('SMI_1', '0.7'))]
    rows = [r for i, r in enumerate(rows)]
    rows[2][6] = 'auc'
    rows[5][6] = 'auc'
    rows.append(['syn', 3, 'SMI_1', 'L', '72', 'hrs', 'fit_auc', '0.1'])
    _exp(tmp_path, 'liver', rows)
    report, _ = vb.validate(str(tmp_path))
    msg = [m for c, _, m in report.warnings if c == 'copied-samples']
    assert msg and "['1', '2']" in msg[0]


def test_na_names_and_cell_line_cancer_wording_do_not_warn(tmp_path):
    rows = [['CVCL_1', 1, 'RH-30', 'RH-30', 'Alveolar Rhabdomyosarcoma', 'cell line', 'Cellosaurus', 'Homo sapiens'],
            ['ACH-1', 1, 'RH-30', 'NA', 'Rhabdomyosarcoma', 'cell line', 'DepMap', 'Homo sapiens']]
    _tsv(tmp_path / 'ccle_samples.csv', rows, SAMPLE_COLS, sep=',')
    _mappings(tmp_path)
    report, _ = vb.validate(str(tmp_path))
    assert 'sample-metadata-conflict' not in _checks(report, 'warnings')
    rows[1][3] = 'SJRH-30'
    _tsv(tmp_path / 'ccle_samples.csv', rows, SAMPLE_COLS, sep=',')
    report, _ = vb.validate(str(tmp_path))
    assert 'sample-metadata-conflict' in _checks(report, 'warnings')


def test_same_size_samples_with_different_data_are_not_copies(tmp_path):
    # PRISM 2.4 false positive: hundreds of samples with the same number of
    # rows collided when per-sample hash sums overflowed in pandas.
    import numpy as np
    rng = np.random.default_rng(0)
    ids = list(range(1, 301))
    _samples(tmp_path, 'liver', ids)
    _drugs(tmp_path, 'liver', [[f'SMI_{d}', f'd{d}', str(d), f'C{d}', 'K', 'F', '1'] for d in range(1, 41)])
    rows = [['syn', s, f'SMI_{d}', 'L', '72', 'hrs', 'fit_auc', f'{rng.random():.4f}']
            for s in ids for d in range(1, 41)]
    _exp(tmp_path, 'liver', rows)
    _mappings(tmp_path)
    report, _ = vb.validate(str(tmp_path))
    assert 'copied-samples' not in _checks(report, 'warnings')
