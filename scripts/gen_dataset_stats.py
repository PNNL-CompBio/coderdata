#!/usr/bin/env python3
"""Generate the per-dataset tables on the docs Datasets page from a local build.

Reads a build output directory (the ``all_files_dir`` produced by
``coderbuild/build_all.py``) and writes three CSVs to the docs static
directory so the documentation site can render current numbers:

* ``dataset_counts.csv``: unique samples, unique drugs, unique genes,
  per-modality row counts, and the dose-response metrics present. The
  Datasets and Modalities table in ``datasets_included.rst`` is copied
  from this file.
* ``dataset_summary_statistics.csv``: unique sample-drug pairs, overall and
  restricted to samples with each omics combination.
* ``dataset_curve_metrics_wide.csv``: drugs per dataset and which
  dose-response metrics each dataset provides.

This script reads the build files directly (no network, no installed
``coderdata`` package required), so it can be run against any local build.

Usage:
    python scripts/gen_dataset_stats.py \
        --build-dir local/all_files_dir \
        --output docs/source/_static/dataset_counts.csv
"""

import argparse
import csv
import gzip
import io
import os
from collections import OrderedDict


# Modalities whose presence/row-count we report, in display order.
OMICS = ["transcriptomics", "proteomics", "phosphoproteomics",
         "mutations", "copy_number"]
DRUGISH = ["drugs", "drug_descriptors", "experiments", "combinations"]

# Column that identifies the unique key for count-distinct on each file type.
UNIQUE_KEY = {
    "samples": "improve_sample_id",
    "drugs": "improve_drug_id",
}
# Files that use TSV rather than CSV.
TSV_TYPES = {"drugs", "drug_descriptors", "experiments", "combinations"}
# Gene-bearing omics used for the unique-gene tally.
GENE_TYPES = ["transcriptomics", "proteomics", "mutations", "copy_number"]


def _open_any(path):
    """Open a plain or .gz text file for reading."""
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def _find_file(build_dir, dataset, modality):
    """Return the path to <dataset>_<modality>.(csv|tsv)[.gz] if present."""
    for ext in (".csv", ".tsv", ".csv.gz", ".tsv.gz"):
        p = os.path.join(build_dir, f"{dataset}_{modality}{ext}")
        if os.path.exists(p):
            return p
    return None


def _delimiter(path):
    return "\t" if (".tsv" in path) else ","


def _header_index(path, column):
    """Return (delimiter, column_index) for a file, or (delim, None) if absent.

    Reads only the first line. Uses fast raw splitting rather than csv module.
    """
    delim = _delimiter(path)
    with _open_any(path) as f:
        header = f.readline().rstrip("\n").rstrip("\r")
    cols = header.split(delim)
    try:
        return delim, cols.index(column)
    except ValueError:
        return delim, None


def _count_rows(path):
    """Number of data rows (excludes header). Fast raw line count."""
    with _open_any(path) as f:
        n = sum(1 for _ in f)
    return max(n - 1, 0)


def _iter_column(path, column):
    """Yield the value of one column per data row using fast line splitting.

    Lines containing a quote character are parsed with the csv module, since
    a quoted field earlier in the row (e.g. a samples ``other_names`` value
    with embedded commas) would otherwise shift the column index. Skips the
    header.
    """
    delim, idx = _header_index(path, column)
    if idx is None:
        return
    with _open_any(path) as f:
        f.readline()  # skip header
        for line in f:
            line = line.rstrip("\n").rstrip("\r")
            if '"' in line:
                parts = next(csv.reader([line], delimiter=delim))
            else:
                parts = line.split(delim)
            if idx < len(parts):
                v = parts[idx]
                if v and v != "NA":
                    yield v


def _count_unique(path, column):
    """Number of unique non-empty values in a column (fast path)."""
    return len(set(_iter_column(path, column)))


def _unique_genes(build_dir, dataset):
    """Unique entrez_id across all gene-bearing omics for a dataset."""
    genes = set()
    for modality in GENE_TYPES:
        p = _find_file(build_dir, dataset, modality)
        if not p:
            continue
        genes.update(_iter_column(p, "entrez_id"))
    return len(genes)


def _dose_response_metrics(build_dir, dataset):
    """Sorted set of dose_response_metric values in the experiments file."""
    p = _find_file(build_dir, dataset, "experiments")
    if not p:
        return []
    return sorted(set(_iter_column(p, "dose_response_metric")))


def _model_types(build_dir, dataset):
    """Distinct model_type values in the samples file, ordered by frequency."""
    p = _find_file(build_dir, dataset, "samples")
    if not p:
        return []
    counts = {}
    for v in _iter_column(p, "model_type"):
        counts[v] = counts.get(v, 0) + 1
    return [k for k, _ in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))]


def discover_datasets(build_dir):
    """Dataset names inferred from <dataset>_samples.csv files."""
    names = set()
    for fn in os.listdir(build_dir):
        if fn.endswith("_samples.csv"):
            names.add(fn[: -len("_samples.csv")])
    return sorted(names)


def build_rows(build_dir, datasets):
    rows = []
    for d in datasets:
        rec = OrderedDict()
        rec["dataset"] = d
        rec["model_type"] = ";".join(_model_types(build_dir, d))

        samples_path = _find_file(build_dir, d, "samples")
        rec["samples"] = _count_unique(samples_path, "improve_sample_id") if samples_path else 0

        drugs_path = _find_file(build_dir, d, "drugs")
        rec["drugs"] = _count_unique(drugs_path, "improve_drug_id") if drugs_path else 0

        rec["genes"] = _unique_genes(build_dir, d)

        exp_path = _find_file(build_dir, d, "experiments")
        rec["experiments"] = _count_rows(exp_path) if exp_path else 0

        # per-modality presence (X / blank) and row counts
        for modality in OMICS + ["drug_descriptors", "combinations"]:
            p = _find_file(build_dir, d, modality)
            rec[modality] = _count_rows(p) if p else 0

        rec["dose_response_metrics"] = ";".join(_dose_response_metrics(build_dir, d))
        rows.append(rec)
    return rows


def _sample_key(v):
    """Normalize an improve_sample_id so '3524.0' and '3524' match.

    Some experiments files (e.g. novartis) write sample ids as floats while
    the omics files write integers.
    """
    return v[:-2] if v.endswith(".0") else v


def _iter_pairs(path):
    """Yield (improve_sample_id, improve_drug_id) per experiments row."""
    delim, s_idx = _header_index(path, "improve_sample_id")
    _, d_idx = _header_index(path, "improve_drug_id")
    with _open_any(path) as f:
        f.readline()  # skip header
        for line in f:
            parts = line.rstrip("\n").rstrip("\r").split(delim)
            yield _sample_key(parts[s_idx]), parts[d_idx]


def build_summary_rows(build_dir, datasets):
    """Unique sample-drug pairs, overall and restricted to samples with omics.

    Mirrors the definitions in notebooks/dataset_documentation.ipynb: a cell
    is blank when the experiments file or any required omics file is absent.
    """
    combos = [
        ("sample_drug_transcript_pairs", ["transcriptomics"]),
        ("sample_drug_transcriptomics_mutation_pairs", ["transcriptomics", "mutations"]),
        ("sample_drug_transcriptomics_copynumber_pairs", ["transcriptomics", "copy_number"]),
        ("sample_drug_mutation_copynumber_pairs", ["mutations", "copy_number"]),
    ]
    rows = []
    for d in datasets:
        rec = OrderedDict(dataset=d, sample_drug_pairs="")
        for name, _ in combos:
            rec[name] = ""
        exp_path = _find_file(build_dir, d, "experiments")
        if exp_path:
            pairs = set(_iter_pairs(exp_path))
            rec["sample_drug_pairs"] = len(pairs)
            omics_samples = {}
            for modality in ("transcriptomics", "mutations", "copy_number"):
                p = _find_file(build_dir, d, modality)
                if p:
                    omics_samples[modality] = {
                        _sample_key(v) for v in _iter_column(p, "improve_sample_id")}
            for name, needed in combos:
                if all(m in omics_samples for m in needed):
                    keep = set.intersection(*(omics_samples[m] for m in needed))
                    rec[name] = sum(1 for s, _ in pairs if s in keep)
        rows.append(rec)
    return rows


def build_curve_metric_rows(build_dir, datasets):
    """Per-dataset drug count and X markers for each dose_response_metric."""
    per_dataset = OrderedDict()
    for d in datasets:
        p = _find_file(build_dir, d, "experiments")
        if not p:
            continue
        delim, d_idx = _header_index(p, "improve_drug_id")
        _, m_idx = _header_index(p, "dose_response_metric")
        drugs, metrics = set(), set()
        with _open_any(p) as f:
            f.readline()  # skip header
            for line in f:
                parts = line.rstrip("\n").rstrip("\r").split(delim)
                drugs.add(parts[d_idx])
                if parts[m_idx] and parts[m_idx] != "NA":
                    metrics.add(parts[m_idx])
        per_dataset[d] = (len(drugs), metrics)
    all_metrics = sorted(set().union(*(m for _, m in per_dataset.values())),
                         key=str.lower)
    rows = []
    for d, (n_drugs, metrics) in per_dataset.items():
        rec = OrderedDict(dataset=d, num_drugs=n_drugs)
        for m in all_metrics:
            rec[m] = "X" if m in metrics else ""
        rows.append(rec)
    return rows


def _write_csv(path, rows):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {path} with {len(rows)} datasets and {len(fieldnames)} columns.")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-dir", default="local/all_files_dir",
                    help="build output directory (default: local/all_files_dir)")
    ap.add_argument("--output", default="docs/source/_static/dataset_counts.csv",
                    help="output CSV path (default: docs/source/_static/dataset_counts.csv)")
    ap.add_argument("--summary-output",
                    default="docs/source/_static/dataset_summary_statistics.csv",
                    help="sample-drug pair summary CSV path "
                         "(default: docs/source/_static/dataset_summary_statistics.csv)")
    ap.add_argument("--curve-output",
                    default="docs/source/_static/dataset_curve_metrics_wide.csv",
                    help="curve metric presence CSV path "
                         "(default: docs/source/_static/dataset_curve_metrics_wide.csv)")
    ap.add_argument("--datasets", nargs="*", default=None,
                    help="restrict to these dataset names (default: all discovered)")
    args = ap.parse_args()

    if not os.path.isdir(args.build_dir):
        raise SystemExit(f"build dir not found: {args.build_dir}")

    datasets = args.datasets or discover_datasets(args.build_dir)
    if not datasets:
        raise SystemExit(f"no *_samples.csv found in {args.build_dir}")

    _write_csv(args.output, build_rows(args.build_dir, datasets))
    _write_csv(args.summary_output, build_summary_rows(args.build_dir, datasets))
    _write_csv(args.curve_output, build_curve_metric_rows(args.build_dir, datasets))
    print("Datasets:", ", ".join(datasets))


if __name__ == "__main__":
    main()
