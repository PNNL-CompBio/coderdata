"""
Tests for drug id assignment in the build (issue #512).

One chemical must get one improve_drug_id, one experiment must map to one
drug, and the stable id step must refuse to merge ids into more than one
value for a sample/drug/time/metric. Fixtures are modelled on the NCI60 cases that
shipped in 2.4: pentamidine (CID 359323) screened as NSC 620107 and NSC
757400, and NSC 694279 listed under two related PubChem records.
"""
import csv
import gzip
import importlib.util
import os
import sys

import pandas as pd
import polars as pl
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UTILS = os.path.join(ROOT, 'coderbuild', 'utils')
BROAD = os.path.join(ROOT, 'coderbuild', 'broad_sanger')
sys.path.insert(0, UTILS)

import map_improve_drug_ids  # noqa: E402
import pubchem_retrieval  # noqa: E402


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(BROAD, filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


nci60_drugs = _load('nci60_drugs', '03a-nci60Drugs.py')
nci60_exp = _load('nci60_exp', '04b-nci60-updated.py')

NCI_SMILES = 'NC(=N)c1ccc(OCCCCCOc2ccc(cc2)C(=N)N)cc1'
PUBCHEM_SMILES = 'C1=CC(=CC=C1C(=N)N)OCCCCCOC2=CC=C(C=C2)C(=N)N.C(CS(=O)(=O)O)O'
DRUG_COLUMNS = ['improve_drug_id', 'chem_name', 'pubchem_id', 'canSMILES', 'InChIKey', 'formula', 'weight']


def _write_tsv(path, rows, columns):
    pd.DataFrame(rows, columns=columns).to_csv(path, sep='\t', index=False)


# --- pubchem_retrieval cache -------------------------------------------------

def test_cache_uses_pubchem_row_and_drops_nsc_names(tmp_path):
    fill_in = ['SMI_9', 'nsc-620107', 359323, NCI_SMILES, '', '', '']
    pubchem = ['SMI_4', 'pentamidine isethionate', 359323, PUBCHEM_SMILES, 'XDRYMKDFEDOLFX-UHFFFAOYSA-N',
               'C21H32N4O10S2', '592.7']
    first, second = tmp_path / 'a.tsv', tmp_path / 'b.tsv'
    _write_tsv(first, [fill_in], DRUG_COLUMNS)
    _write_tsv(second, [pubchem], DRUG_COLUMNS)

    pubchem_retrieval._pubchem_data_cache.clear()
    pubchem_retrieval.load_pubchem_cache([str(first), str(second)])
    entry = pubchem_retrieval._pubchem_data_cache['359323']

    assert entry['SMILES'] == PUBCHEM_SMILES
    assert 'nsc-620107' not in entry['synonyms']
    assert 'nsc-620107' not in pubchem_retrieval._pubchem_data_cache
    assert pubchem_retrieval._pubchem_data_cache['pentamidine isethionate'] is entry


def test_cache_skips_drugs_named_only_by_nsc(tmp_path):
    # CID 24202157's only PubChem synonym is 'nsc-630659'. Cached without it,
    # the entry had no names, wrote no rows, and the drug was dropped (113
    # NCI60 drugs in the first 2.4 rebuild). It must be fetched live instead.
    path = tmp_path / 'a.tsv'
    _write_tsv(path, [['SMI_7', 'nsc-630659', 24202157, 'N.[Co+3]', 'KEY', 'F', '1']], DRUG_COLUMNS)
    pubchem_retrieval._pubchem_data_cache.clear()
    pubchem_retrieval.load_pubchem_cache([str(path)])
    assert '24202157' not in pubchem_retrieval._pubchem_data_cache


def test_rerun_keeps_drugs_already_written(tmp_path, monkeypatch):
    # 03a reruns after a failed later step find their own earlier output. All
    # CIDs are requested again; already-written ones must be kept, not fetched.
    out = tmp_path / 'nci60_drugs.tsv'
    rows = [['SMI_1', 'pentamidine', '359323', PUBCHEM_SMILES, 'K1', 'F', '1'],
            ['SMI_2', 'nsc-630659', '24202157', 'N.[Co+3]', 'K2', 'F', '1']]
    _write_tsv(out, rows, DRUG_COLUMNS)
    monkeypatch.setattr(pubchem_retrieval, 'fetch_url',
                        lambda url: pytest.fail(f'already written, should not fetch {url}'))
    pubchem_retrieval._pubchem_data_cache.clear()
    pubchem_retrieval.update_dataframe_and_write_tsv(
        unique_names={'359323', '24202157'}, output_filename=str(out),
        ignore_chems=str(tmp_path / 'ignore.txt'), isname=False,
        restrict_to_raw_names={'359323', '24202157'})
    kept = pd.read_csv(out, sep='\t', dtype=str)
    assert set(kept['improve_drug_id']) == {'SMI_1', 'SMI_2'}


def test_cache_entry_does_not_depend_on_file_order(tmp_path):
    rows = [['SMI_9', 'pentamidine', 359323, NCI_SMILES, '', '', ''],
            ['SMI_4', 'pentamidine isethionate', 359323, PUBCHEM_SMILES, 'KEY', 'F', '1']]
    a, b = tmp_path / 'a.tsv', tmp_path / 'b.tsv'
    _write_tsv(a, rows[:1], DRUG_COLUMNS)
    _write_tsv(b, rows[1:], DRUG_COLUMNS)
    smiles = []
    for order in ([a, b], [b, a]):
        pubchem_retrieval._pubchem_data_cache.clear()
        pubchem_retrieval.load_pubchem_cache([str(p) for p in order])
        smiles.append(pubchem_retrieval._pubchem_data_cache['359323']['SMILES'])
    assert smiles == [PUBCHEM_SMILES, PUBCHEM_SMILES]


# --- 03a-nci60Drugs: NSC names missing from PubChem ---------------------------

def test_missing_nsc_reuses_existing_drug_id():
    res = pl.DataFrame({'improve_drug_id': ['SMI_4'], 'chem_name': ['pentamidine isethionate'],
                        'pubchem_id': [359323], 'canSMILES': [PUBCHEM_SMILES], 'InChIKey': ['KEY'],
                        'formula': ['F'], 'weight': ['1']})
    mdf = pl.DataFrame({'NSC': [620107, 757400], 'NAME': ['Pentamidine isethionate', 'Pentamidine'],
                        'CID': [359323, 359323], 'SMILES': [NCI_SMILES, NCI_SMILES]})
    out = nci60_drugs.add_missing_nscs(res, mdf)

    assert set(out['improve_drug_id']) == {'SMI_4'}
    assert {'nsc-620107', 'nsc-757400'} <= set(out['chem_name'])
    assert set(out['canSMILES']) == {PUBCHEM_SMILES}


def test_missing_nsc_handles_float_cids_and_taken_names():
    # pubchem_id read as float (a blank somewhere in the file) must still match,
    # and an NCI name already used by another compound must not be reused.
    res = pl.DataFrame({'improve_drug_id': ['SMI_4', 'SMI_7'], 'chem_name': ['pentamidine', 'cisplatin'],
                        'pubchem_id': [359323.0, 5702198.0], 'canSMILES': ['A', 'B'],
                        'InChIKey': ['K1', 'K2'], 'formula': ['F', 'F'], 'weight': ['1', '1']})
    mdf = pl.DataFrame({'NSC': [757400], 'NAME': ['Cisplatin'], 'CID': [359323], 'SMILES': ['X']})
    out = nci60_drugs.add_missing_nscs(res, mdf)
    added = out.filter(pl.col('chem_name') == 'nsc-757400')
    assert added['improve_drug_id'].to_list() == ['SMI_4']
    assert out.filter(pl.col('chem_name') == 'cisplatin')['improve_drug_id'].to_list() == ['SMI_7']


def test_every_nsc_gets_a_dash_name_from_the_nci_table():
    # PubChem lists the NSC as 'nsc620107' only; 04b needs 'nsc-620107'.
    # NSC 99 has a CID not in the drug table and is skipped.
    res = pl.DataFrame({'improve_drug_id': ['SMI_4', 'SMI_4'], 'chem_name': ['pentamidine', 'nsc620107'],
                        'pubchem_id': [359323, 359323], 'canSMILES': ['A', 'A'], 'InChIKey': ['K', 'K'],
                        'formula': ['F', 'F'], 'weight': ['1', '1']})
    nsc_cid = pl.DataFrame({'NSC': [620107, 757400, 99], 'CID': [359323.0, 359323.0, 1.0]})
    out = nci60_drugs.add_nsc_names(res, nsc_cid)

    dash = out.filter(pl.col('chem_name').str.starts_with('nsc-'))
    assert sorted(dash['chem_name']) == ['nsc-620107', 'nsc-757400']
    assert set(dash['improve_drug_id']) == {'SMI_4'}
    assert set(dash['canSMILES']) == {'A'}
    assert nci60_drugs.add_nsc_names(out, nsc_cid).height == out.height  # rerun adds nothing


def test_cache_keys_are_integer_cids_even_with_blank_ids(tmp_path):
    rows = [['SMI_4', 'pentamidine', 359323, PUBCHEM_SMILES, 'KEY', 'F', '1'],
            ['SMI_9', 'unknown drug', '', 'CCO', '', '', '']]
    path = tmp_path / 'a.tsv'
    _write_tsv(path, rows, DRUG_COLUMNS)
    pubchem_retrieval._pubchem_data_cache.clear()
    pubchem_retrieval.load_pubchem_cache([str(path)])
    assert pubchem_retrieval._pubchem_data_cache['359323']['CID'] == '359323'
    assert '359323.0' not in pubchem_retrieval._pubchem_data_cache


# --- 04b-nci60-updated: one drug per NSC ---------------------------------------

def _nsc_inputs():
    drugmapping = pl.DataFrame({'NSC': [694279, 694279, 1], 'improve_drug_id': ['SMI_1', 'SMI_2', 'SMI_3']})
    drugs = pl.DataFrame({'improve_drug_id': ['SMI_1', 'SMI_2', 'SMI_3'],
                          'chem_name': ['nsc-694279', 'nsc-694279', 'nsc-1'],
                          'pubchem_id': [6712523, 123664, 5]})
    return drugmapping, drugs


def test_ambiguous_nsc_resolved_with_nci_cid_table():
    drugmapping, drugs = _nsc_inputs()
    nci = lambda: pl.DataFrame({'NSC': ['694279'], 'CID': ['123664']})
    out = nci60_exp.resolve_nsc_ids(drugmapping, drugs, nci)
    assert sorted(out.rows()) == [(1, 'SMI_3'), (694279, 'SMI_2')]


def test_nsc_under_two_ids_with_the_same_cid_keeps_the_lowest():
    # One chemical listed twice in the drug table: both ids carry NCI's CID.
    drugmapping = pl.DataFrame({'NSC': [179940, 179940], 'improve_drug_id': ['SMI_56129', 'SMI_6490']})
    drugs = pl.DataFrame({'improve_drug_id': ['SMI_56129', 'SMI_6490'], 'pubchem_id': [42640, 42640]})
    nci = lambda: pl.DataFrame({'NSC': ['179940'], 'CID': ['42640']})
    out = nci60_exp.resolve_nsc_ids(drugmapping, drugs, nci)
    assert out.rows() == [(179940, 'SMI_6490')]


def test_ambiguous_nsc_resolved_with_float_cids():
    drugmapping, drugs = _nsc_inputs()
    drugs = drugs.with_columns(pl.col('pubchem_id').cast(pl.Float64))
    nci = lambda: pl.DataFrame({'NSC': ['694279'], 'CID': ['123664.0']})
    out = nci60_exp.resolve_nsc_ids(drugmapping, drugs, nci)
    assert sorted(out.rows()) == [(1, 'SMI_3'), (694279, 'SMI_2')]


def test_unresolvable_nsc_stops_the_build():
    drugmapping, drugs = _nsc_inputs()
    nci = lambda: pl.DataFrame({'NSC': ['694279'], 'CID': ['999']})
    with pytest.raises(ValueError, match='more than one drug'):
        nci60_exp.resolve_nsc_ids(drugmapping, drugs, nci, fetch_inchikeys=lambda cids: {})


def test_nsc_falls_back_to_inchikey_skeleton():
    # NSC 141993: NCI lists minocycline HCl without stereo (CID 54706380);
    # the drug table has the stereo salt record and the free base.
    drugmapping = pl.DataFrame({'NSC': [141993, 141993], 'improve_drug_id': ['SMI_16296', 'SMI_57059']})
    drugs = pl.DataFrame({'improve_drug_id': ['SMI_16296', 'SMI_57059'],
                          'pubchem_id': [54685925, 54675783],
                          'InChIKey': ['WTJXVDPDEQKTCV-VQAITOIOSA-N', 'FFTVPQUHLQBXQZ-KVUCHLLUSA-N']})
    nci = lambda: pl.DataFrame({'NSC': ['141993'], 'CID': ['54706380']})
    keys = lambda cids: {'54706380': 'WTJXVDPDEQKTCV-UHFFFAOYSA-N'}
    out = nci60_exp.resolve_nsc_ids(drugmapping, drugs, nci, fetch_inchikeys=keys)
    assert out.rows() == [(141993, 'SMI_16296')]


def test_unambiguous_nscs_do_not_load_nci_table():
    drugmapping = pl.DataFrame({'NSC': [1], 'improve_drug_id': ['SMI_3']})
    def fail():
        raise AssertionError('NCI table should not be needed')
    out = nci60_exp.resolve_nsc_ids(drugmapping, pl.DataFrame({'improve_drug_id': ['SMI_3'], 'pubchem_id': [5]}), fail)
    assert out.equals(drugmapping)


# --- map_improve_drug_ids: refuse to publish two values for one key ------------

EXP_COLUMNS = ['source', 'improve_sample_id', 'improve_drug_id', 'study', 'time', 'time_unit',
               'dose_response_metric', 'dose_response_value']


def _exp_file(tmp_path, rows, gz=False):
    path = tmp_path / 'nci60_experiments.tsv'
    _write_tsv(path, rows, EXP_COLUMNS)
    if gz:
        with open(path, 'rb') as f_in, gzip.open(str(path) + '.gz', 'wb') as f_out:
            f_out.write(f_in.read())
        os.remove(path)
        return str(path) + '.gz'
    return str(path)


def test_rewrite_raises_when_two_ids_merge_onto_one_key(tmp_path):
    path = _exp_file(tmp_path, [['NCI60_24', '46', 'SMI_9', 'NCI60', '48', 'hours', 'fit_auc', '0.8265'],
                                ['NCI60_24', '46', 'SMI_8', 'NCI60', '48', 'hours', 'fit_auc', '0.7799']])
    mapping = {('nci60', 'SMI_9'): 'SMI_13316', ('nci60', 'SMI_8'): 'SMI_13316'}
    with pytest.raises(RuntimeError, match='same key'):
        map_improve_drug_ids.rewrite_drugs_file(path, mapping)
    # the file is left untouched
    with open(path) as f:
        assert [r['improve_drug_id'] for r in csv.DictReader(f, delimiter='\t')] == ['SMI_9', 'SMI_8']


def test_rewrite_without_collision_still_maps_ids(tmp_path):
    path = _exp_file(tmp_path, [['NCI60_24', '46', 'SMI_9', 'NCI60', '48', 'hours', 'fit_auc', '0.8'],
                                ['NCI60_24', '47', 'SMI_8', 'NCI60', '48', 'hours', 'fit_auc', '0.7']], gz=True)
    mapping = {('nci60', 'SMI_9'): 'SMI_13316', ('nci60', 'SMI_8'): 'SMI_13316'}
    map_improve_drug_ids.rewrite_drugs_file(path, mapping)
    with gzip.open(path, 'rt') as f:
        assert {r['improve_drug_id'] for r in csv.DictReader(f, delimiter='\t')} == {'SMI_13316'}


def test_rewrite_raises_on_colliding_combination_rows(tmp_path):
    columns = ['improve_drug_1', 'improve_drug_2', 'improve_sample_id', 'time_unit', 'time',
               'drug_combination_metric', 'drug_combination_value', 'source', 'study']
    path = tmp_path / 'mpnst_combinations.tsv'
    _write_tsv(path, [['SMI_9', 'SMI_1', '3971', 'days', '25', 'lmm', '-0.1', 'NF', 'PDX'],
                      ['SMI_8', 'SMI_1', '3971', 'days', '25', 'lmm', '-0.2', 'NF', 'PDX']], columns)
    mapping = {('mpnst', 'SMI_9'): 'SMI_5', ('mpnst', 'SMI_8'): 'SMI_5'}
    with pytest.raises(RuntimeError, match='same key'):
        map_improve_drug_ids.rewrite_drugs_file(str(path), mapping)
