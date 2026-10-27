#!/usr/bin/env python3
"""
Data checks on a CoderData build that the LinkML schema cannot express.

The schema checks that each row is well formed. These checks look across
rows, files and releases, at the problems that shipped in 2.4 without failing
the schema: two results for one sample/drug (NCI60, PRISM), one experiment
reported under two drugs, drug responses attached to the wrong specimens
(BeatAML) and whole sample classes dropped (HCMI).

Errors are definite integrity violations and make the script exit 1.
Warnings need a person to look; they never fail the build.

Usage:
    python coderbuild/utils/validate_build.py --dir local/all_files_dir \\
        [--previous-summary prev_summary.json] [--write-summary this_summary.json]

Results are printed to the terminal.

--previous-summary compares this build with an earlier one (a summary written
by --write-summary; each release publishes one as validation_summary.json.gz) and warns on large
changes, drops of more than half and datasets that are no longer built.
Files may be plain or .gz.
"""
import argparse
import glob
import gzip
import hashlib
import json
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

CELL_LINE_DATASETS = {'ccle', 'ctrpv2', 'fimm', 'gcsi', 'gdscv1', 'gdscv2', 'nci60', 'prism'}
OMICS = ['transcriptomics', 'proteomics', 'phosphoproteomics', 'mutations', 'copy_number']
CHUNK_ROWS = 2_000_000

EXPERIMENT_KEY = ['source', 'improve_sample_id', 'improve_drug_id', 'study',
                  'time', 'time_unit', 'dose_response_metric']
COMBINATION_KEY = ['source', 'improve_sample_id', 'improve_drug_1', 'improve_drug_2',
                   'study', 'time', 'time_unit', 'drug_combination_metric']
# Metrics on a 0..~1 scale; values far outside are worth a look.
BOUNDED_METRICS = {'auc': (0, 1.5), 'fit_auc': (0, 1.5), 'aac': (-0.5, 1.5)}
# Relative change versus the previous build that triggers a warning.
CHANGE_THRESHOLD = 0.10
# A drop this large versus the previous build is reported on its own: it is far
# more likely a lost download or broken join than a real change in the data.
DROP_THRESHOLD = 0.5
# Values that stand for "no value" in sample metadata.
PLACEHOLDERS = {'', 'na', 'nan', 'none', 'not reported'}


# --------------------------------------------------------------------------- io

def _find(directory, dataset, kind):
    for ext in ('.csv', '.tsv', '.csv.gz', '.tsv.gz'):
        path = os.path.join(directory, f'{dataset}_{kind}{ext}')
        if os.path.exists(path):
            return path
    return None


def _sep(path):
    return '\t' if '.tsv' in path else ','


def _columns(path):
    return list(pd.read_csv(path, sep=_sep(path), nrows=0).columns)


def _chunks(path, columns):
    present = [c for c in _columns(path) if c in columns]
    return pd.read_csv(path, sep=_sep(path), usecols=present, dtype=str,
                       keep_default_na=False, chunksize=CHUNK_ROWS)


def _read(path, columns):
    return pd.concat(_chunks(path, columns), ignore_index=True)


def _unique(path, column, keep_blank=False):
    values = set()
    if column not in _columns(path):
        return values
    for chunk in _chunks(path, [column]):
        values.update(chunk[column].unique())
    if not keep_blank:
        values.discard('')
    return values


def discover_datasets(directory):
    names = set()
    for path in glob.glob(os.path.join(directory, '*_samples.csv*')):
        names.add(os.path.basename(path).split('_samples.csv')[0])
    return sorted(names)


# ----------------------------------------------------------------------- report

class Report:
    def __init__(self):
        self.errors = []
        self.warnings = []

    def error(self, check, dataset, message):
        self.errors.append((check, dataset, message))

    def warn(self, check, dataset, message):
        self.warnings.append((check, dataset, message))

    def markdown(self):
        lines = ['# CoderData build validation', '',
                 f'{len(self.errors)} error(s), {len(self.warnings)} warning(s).', '']
        for title, items in (('Errors', self.errors), ('Warnings', self.warnings)):
            lines.append(f'## {title}')
            lines.append('')
            if not items:
                lines.append('None.')
            for check, dataset, message in items:
                lines.append(f'- **{check}** ({dataset}): {message}')
            lines.append('')
        return '\n'.join(lines)


def _examples(values, n=5):
    values = list(values)
    shown = ', '.join(str(v) for v in values[:n])
    return shown + (f', ... ({len(values)} total)' if len(values) > n else '')


# ----------------------------------------------------------------------- checks

def duplicate_keys(path, key, report, check, dataset):
    """Rows sharing every key column (two values for one sample/drug/metric)."""
    cols = _columns(path)
    missing = [c for c in key if c not in cols]
    if missing:
        report.error(check, dataset, f'{os.path.basename(path)} lacks key column(s) {missing}')
        return
    hashes = [pd.util.hash_pandas_object(c[key], index=False).to_numpy() for c in _chunks(path, key)]
    if not hashes:
        return
    values, counts = np.unique(np.concatenate(hashes), return_counts=True)
    suspect = values[counts > 1]
    if not len(suspect):
        return
    rows = pd.concat([c[key][np.isin(pd.util.hash_pandas_object(c[key], index=False).to_numpy(), suspect)]
                      for c in _chunks(path, key)], ignore_index=True)
    dup = rows[rows.duplicated(key, keep=False)]
    if len(dup):
        pairs = dup.drop_duplicates().apply(lambda r: f"{r['improve_sample_id']}/"
                                            f"{r.get('improve_drug_id', r.get('improve_drug_1'))}", axis=1)
        report.error(check, dataset,
                     f'{len(dup)} rows in {os.path.basename(path)} share a key with another row '
                     f'(more than one value for a sample/drug/time/metric); sample/drug e.g. '
                     f'{_examples(pairs.unique())}')


def copied_responses(path, report, dataset, metric='fit_auc', min_rows=3):
    """Different drug ids with identical responses in every sample."""
    cols = ['improve_sample_id', 'improve_drug_id', 'time', 'dose_response_metric', 'dose_response_value']
    parts = [c[c['dose_response_metric'] == metric] for c in _chunks(path, cols)]
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    if df.empty:
        return
    key = df['improve_sample_id'] + '|' + df['time'] + '|' + df['dose_response_value']
    sig = key.groupby(df['improve_drug_id']).agg(
        lambda s: hashlib.md5('\n'.join(sorted(s)).encode()).hexdigest() if len(s) >= min_rows else '')
    sig = sig[sig != '']
    shared = sig[sig.duplicated(keep=False)]
    if len(shared):
        groups = [sorted(g.index) for _, g in shared.groupby(shared)]
        report.warn('copied-responses', dataset,
                    f'{len(groups)} group(s) of drug ids report identical {metric} for every sample '
                    f'(one experiment attached to several drugs?): {_examples(groups)}')


def value_ranges(path, report, dataset):
    """Bounded metrics far outside their range, and metrics that are mostly blank."""
    cols = ['dose_response_metric', 'dose_response_value']
    total, blank, out = defaultdict(int), defaultdict(int), defaultdict(int)
    for chunk in _chunks(path, cols):
        values = pd.to_numeric(chunk['dose_response_value'], errors='coerce')
        for metric, idx in chunk.groupby('dose_response_metric').groups.items():
            v = values.loc[idx]
            total[metric] += len(v)
            blank[metric] += int(v.isna().sum() + np.isinf(v).sum())
            if metric in BOUNDED_METRICS:
                lo, hi = BOUNDED_METRICS[metric]
                out[metric] += int(((v < lo) | (v > hi)).sum())
    for metric in total:
        if blank[metric] == total[metric]:
            continue  # categorical metric (e.g. mrecist: PD/SD/PR/CR)
        if blank[metric] / total[metric] > 0.5:
            report.warn('value-ranges', dataset,
                        f'{metric}: {blank[metric] / total[metric]:.0%} of values are blank or infinite')
        if out.get(metric):
            lo, hi = BOUNDED_METRICS[metric]
            report.warn('value-ranges', dataset,
                        f'{metric}: {out[metric]} of {total[metric]} values outside [{lo}, {hi}]')


def _id_values(report, dataset, kind, where, values):
    """Report blank ids and ids written as decimals; return the ids without blanks."""
    if '' in values:
        report.error('missing-id', dataset, f'{where} has rows with a blank {kind} id')
    decimal = sorted(v for v in values if re.fullmatch(r'\d+\.0', v))
    if decimal:
        report.error('decimal-id', dataset,
                     f'{len(decimal)} {kind} ids in {where} are written as decimals (e.g. "{decimal[0]}")')
    return values - {''}


def _orphans(report, dataset, kind, where, missing, known):
    if not missing:
        return
    decimal = {m for m in missing if re.fullmatch(r'\d+\.0', m) and m[:-2] in known}
    hint = (f'; {len(decimal)} of them are known ids written as decimals (e.g. "{sorted(decimal)[0]}")'
            if decimal else '')
    report.error(f'orphan-{kind}-id', dataset,
                 f'{len(missing)} {kind} ids in {where} are not in the {kind}s file: '
                 f'{_examples(sorted(missing))}{hint}')


def sample_checks(dataset, path, report):
    """Per-sample metadata problems; returns the dataset's sample ids."""
    wanted = ['improve_sample_id', 'other_id', 'other_id_source', 'model_type',
              'cancer_type', 'common_name', 'other_names']
    s = _read(path, wanted).reindex(columns=wanted, fill_value='')
    ids = _id_values(report, dataset, 'sample', 'samples', set(s['improve_sample_id']))
    mt = s.groupby('improve_sample_id')['model_type'].nunique()
    if (mt > 1).any():
        report.error('sample-model-type', dataset,
                     f'{int((mt > 1).sum())} sample ids carry more than one model_type: '
                     f'{_examples(mt[mt > 1].index)}')
    # One sample, conflicting descriptions. other_names are aliases and may
    # differ; 'NA' is a placeholder, not a name. Cell lines combine DepMap and
    # Sanger records whose cancer_type vocabularies differ, so only their names
    # are compared.
    cols = ('common_name',) if dataset in CELL_LINE_DATASETS else ('common_name', 'cancer_type')
    for col in cols:
        named = s[~s[col].str.strip().str.lower().isin(PLACEHOLDERS)]
        n = named.groupby('improve_sample_id')[col].nunique()
        if (n > 1).any():
            report.warn('sample-metadata-conflict', dataset,
                        f'{int((n > 1).sum())} sample ids have more than one {col}: {_examples(n[n > 1].index)}')
    return ids, s


def copied_samples(path, value_cols, report, dataset, what, min_rows=3):
    """
    Different sample ids whose rows in a file are identical.

    Each sample gets an order-independent fingerprint: two independent 64-bit
    row hashes summed with wrap-around (numpy uint64; pandas sums can spill
    into float64 and collide), plus the row count. The file is streamed, and
    samples whose fingerprints match are then compared row by row, so only
    exact copies are reported.
    """
    cols = ['improve_sample_id'] + value_cols
    if not all(c in _columns(path) for c in cols):
        return
    keys = ('0123456789123456', 'fedcba9876543210')
    acc = {}  # sample -> [sum1, sum2, rows]
    for chunk in _chunks(path, cols):
        sid = chunk['improve_sample_id'].str.replace(r'\.0$', '', regex=True)
        codes, uniq = pd.factorize(sid)
        n = np.bincount(codes, minlength=len(uniq))
        sums = []
        for key in keys:
            h = pd.util.hash_pandas_object(chunk[value_cols], index=False, hash_key=key).to_numpy()
            s = np.zeros(len(uniq), dtype=np.uint64)
            np.add.at(s, codes, h)
            sums.append(s)
        for k, sample in enumerate(uniq):
            a = acc.setdefault(sample, np.zeros(3, dtype=np.uint64))
            a += np.array((sums[0][k], sums[1][k], n[k]), dtype=np.uint64)  # array add wraps silently
    if not acc:
        return
    sig = pd.DataFrame.from_dict(acc, orient='index', columns=['h1', 'h2', 'n'])
    sig = sig[sig['n'] >= min_rows]
    suspect = sig[sig.duplicated(['h1', 'h2', 'n'], keep=False)]
    if suspect.empty:
        return
    # Confirm exactly: re-read only the suspect samples and compare their rows.
    wanted = set(suspect.index)
    rows = pd.concat([c.assign(improve_sample_id=c['improve_sample_id'].str.replace(r'\.0$', '', regex=True))
                      .loc[lambda x: x['improve_sample_id'].isin(wanted)]
                      for c in _chunks(path, cols)], ignore_index=True)
    content = {s: tuple(map(tuple, g[value_cols].sort_values(value_cols).to_numpy()))
               for s, g in rows.groupby('improve_sample_id')}
    groups = []
    for _, g in suspect.groupby(['h1', 'h2', 'n']):
        by_content = defaultdict(list)
        for s in g.index:
            by_content[content[s]].append(s)
        groups += [sorted(v) for v in by_content.values() if len(v) > 1]
    if groups:
        report.warn('copied-samples', dataset,
                    f'{len(groups)} group(s) of sample ids have identical {what}: {_examples(sorted(groups))}')


def drug_checks(dataset, path, report):
    """One structure per drug id and one drug id per CID; returns the drug ids.

    One CID under two drug ids is one chemical listed twice; its responses
    are split between the ids and merging them later gives duplicate results
    (issue #512), so it is an error.
    """
    d = _read(path, ['improve_drug_id', 'canSMILES', 'pubchem_id'])
    ids = _id_values(report, dataset, 'drug', 'drugs', set(d['improve_drug_id']))
    smiles = d[d['canSMILES'] != ''].groupby('improve_drug_id')['canSMILES'].nunique()
    if (smiles > 1).any():
        report.error('drug-structure', dataset,
                     f'{int((smiles > 1).sum())} drug ids have more than one canSMILES: '
                     f'{_examples(smiles[smiles > 1].index)}')
    if 'pubchem_id' in d.columns:
        cid = d[d['pubchem_id'] != ''].assign(pubchem_id=lambda x: x['pubchem_id'].str.replace(r'\.0$', '', regex=True))
        per_cid = cid.groupby('pubchem_id')['improve_drug_id'].agg(lambda x: sorted(set(x)))
        split = per_cid[per_cid.apply(len) > 1]
        if len(split):
            report.error('drug-split', dataset,
                        f'{len(split)} PubChem CID(s) sit under more than one drug id: '
                        f'{_examples(f"{k}: {v}" for k, v in split.items())}')
    return ids


def mapping_checks(directory, report):
    for name, items_key, key_fn in (
            ('improve_sample_mapping.json', 'samples',
             lambda e: [tuple(q.get(k, '').strip() for k in ('other_id', 'other_id_source', 'model_type', 'other_names'))
                        for q in e.get('quadruples', [])]),
            ('improve_drug_mapping.json', 'drugs', lambda e: [e.get('canSMILES')])):
        path = next((p for p in (os.path.join(directory, name), os.path.join(directory, name + '.gz'))
                     if os.path.exists(p)), None)
        if not path:
            report.warn('mapping-file', 'all', f'{name} not found')
            continue
        opener = gzip.open if path.endswith('.gz') else open
        with opener(path, 'rt') as f:
            entries = json.load(f).get(items_key, [])
        ids = [e['stable_id'] for e in entries]
        if len(ids) != len(set(ids)):
            report.error('mapping-file', 'all', f'{name}: stable ids listed more than once')
        owner = defaultdict(set)
        for e in entries:
            for k in key_fn(e):
                owner[k].add(e['stable_id'])
        multi = {k: v for k, v in owner.items() if len(v) > 1}
        if multi:
            report.error('mapping-file', 'all',
                         f'{name}: {len(multi)} identities map to more than one stable id, e.g. '
                         f'{_examples(f"{k}: {sorted(v)}" for k, v in multi.items())}')


# ------------------------------------------------------------------ main runner

def validate(directory, datasets=None, previous=None):
    report = Report()
    requested = datasets
    datasets = datasets or discover_datasets(directory)
    summary = {'datasets': {}}
    sample_owner = defaultdict(set)

    for ds in datasets:
        info = {}
        spath = _find(directory, ds, 'samples')
        if not spath:
            report.error('missing-file', ds, 'no samples file')
            continue
        sids, samples = sample_checks(ds, spath, report)
        info['samples'] = len(sids)
        info['model_types'] = samples.drop_duplicates('improve_sample_id')['model_type'].value_counts().to_dict()
        for sid in sids:
            sample_owner[sid].add(ds)

        dpath = _find(directory, ds, 'drugs')
        dids = drug_checks(ds, dpath, report) if dpath else set()
        info['drugs'] = len(dids)

        info['omics'] = {}
        for kind in OMICS:
            path = _find(directory, ds, kind)
            if not path:
                continue
            used = _id_values(report, ds, 'sample', kind, _unique(path, 'improve_sample_id', keep_blank=True))
            if kind == 'transcriptomics':
                copied_samples(path, ['entrez_id', 'transcriptomics'], report, ds, 'transcriptomics')
            info['omics'][kind] = len(used)
            _orphans(report, ds, 'sample', kind, used - sids, sids)

        epath = _find(directory, ds, 'experiments')
        if epath:
            duplicate_keys(epath, EXPERIMENT_KEY, report, 'duplicate-result', ds)
            copied_responses(epath, report, ds)
            copied_samples(epath, ['improve_drug_id', 'time', 'dose_response_metric', 'dose_response_value'],
                           report, ds, 'drug responses')
            value_ranges(epath, report, ds)
            pairs, esids, edids = set(), set(), set()
            for chunk in _chunks(epath, ['improve_sample_id', 'improve_drug_id']):
                pairs.update(pd.util.hash_pandas_object(chunk[['improve_sample_id', 'improve_drug_id']],
                                                        index=False).to_numpy().tolist())
                esids.update(chunk['improve_sample_id'].unique())
                edids.update(chunk['improve_drug_id'].unique())
            esids = _id_values(report, ds, 'sample', 'experiments', esids)
            edids = _id_values(report, ds, 'drug', 'experiments', edids)
            info['pairs'] = len(pairs)
            info['samples_with_response'] = len(esids)
            info['drugs_with_response'] = len(edids)
            _orphans(report, ds, 'sample', 'experiments', esids - sids, sids)
            _orphans(report, ds, 'drug', 'experiments', edids - dids, dids)
            no_omics = [k for k, n in info['omics'].items() if esids and n == 0]
            if no_omics:
                report.warn('omics-coverage', ds, f'no samples with drug response have {no_omics}')

        cpath = _find(directory, ds, 'combinations')
        if cpath:
            duplicate_keys(cpath, COMBINATION_KEY, report, 'duplicate-result', ds)
            c = _read(cpath, ['improve_sample_id', 'improve_drug_1', 'improve_drug_2'])
            csids = _id_values(report, ds, 'sample', 'combinations', set(c['improve_sample_id']))
            cdids = _id_values(report, ds, 'drug', 'combinations',
                               set(c['improve_drug_1']) | set(c['improve_drug_2']))
            _orphans(report, ds, 'sample', 'combinations', csids - sids, sids)
            _orphans(report, ds, 'drug', 'combinations', cdids - dids, dids)
            info['combination_rows'] = len(c)

        xpath = _find(directory, ds, 'drug_descriptors')
        if xpath and dids:
            xdids = _id_values(report, ds, 'drug', 'drug_descriptors',
                               _unique(xpath, 'improve_drug_id', keep_blank=True))
            _orphans(report, ds, 'drug', 'drug_descriptors', xdids - dids, dids)
        summary['datasets'][ds] = info

    # Cell lines share ids by design; every other dataset owns its ids.
    leaked = {sid: sorted(owners) for sid, owners in sample_owner.items()
              if len(owners) > 1 and not owners <= CELL_LINE_DATASETS}
    if leaked:
        report.error('sample-id-overlap', 'all',
                     f'{len(leaked)} sample ids are used by more than one dataset: '
                     f'{_examples(f"{k}: {v}" for k, v in leaked.items())}')

    mapping_checks(directory, report)
    if previous:
        compare(summary, previous, report, requested)
    return report, summary


def compare(summary, previous, report, datasets=None):
    """Warn on large changes versus a previous build's summary."""
    prev = previous.get('datasets', {})
    # Datasets left out on purpose (--datasets) are not missing.
    expected = set(prev) if datasets is None else set(prev) & set(datasets)
    for ds in sorted(expected - set(summary['datasets'])):
        report.warn('dataset-missing', ds, 'in the previous build but not in this one')
    for ds, cur in summary['datasets'].items():
        old = prev.get(ds)
        if old is None:
            report.warn('release-change', ds, 'new dataset (not in the previous build)')
            continue
        for field in ('samples', 'drugs', 'pairs', 'samples_with_response', 'drugs_with_response'):
            a, b = old.get(field), cur.get(field)
            if a and b is not None and (a - b) / a > DROP_THRESHOLD:
                report.warn('release-drop', ds, f'{field}: {a:,} -> {b:,} ({(b - a) / a:+.0%})')
            elif a and b is not None and abs(b - a) / a > CHANGE_THRESHOLD:
                report.warn('release-change', ds, f'{field}: {a:,} -> {b:,} ({(b - a) / a:+.0%})')
        for kind in set(old.get('omics', {})) | set(cur.get('omics', {})):
            a, b = old.get('omics', {}).get(kind), cur.get('omics', {}).get(kind)
            if a is None or b is None:
                report.warn('release-change', ds, f'{kind}: {"added" if a is None else "removed"}')
            elif a and (a - b) / a > DROP_THRESHOLD:
                report.warn('release-drop', ds, f'{kind} samples: {a:,} -> {b:,} ({(b - a) / a:+.0%})')
            elif a and abs(b - a) / a > CHANGE_THRESHOLD:
                report.warn('release-change', ds, f'{kind} samples: {a:,} -> {b:,} ({(b - a) / a:+.0%})')
        old_mt, new_mt = set(old.get('model_types', {})), set(cur.get('model_types', {}))
        if old_mt != new_mt:
            report.warn('release-change', ds,
                        f'model types changed: removed {sorted(old_mt - new_mt) or "none"}, '
                        f'added {sorted(new_mt - old_mt) or "none"}')


def load_summary(path):
    """Read a summary written by --write-summary (plain or .gz, as published)."""
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as f:
        return json.load(f)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', required=True, help='folder of build output files')
    ap.add_argument('--datasets', nargs='*', default=None, help='datasets to check (default: all found)')
    ap.add_argument('--previous-summary', default=None, help='summary JSON of the build to compare with')
    ap.add_argument('--write-summary', default=None, help='write this build\'s summary JSON here')
    args = ap.parse_args(argv)

    previous = load_summary(args.previous_summary) if args.previous_summary else None
    report, summary = validate(args.dir, args.datasets, previous)
    print(report.markdown())
    if args.write_summary:
        with open(args.write_summary, 'w') as f:
            json.dump(summary, f, indent=1, sort_keys=True)
    return 1 if report.errors else 0


if __name__ == '__main__':
    sys.exit(main())
