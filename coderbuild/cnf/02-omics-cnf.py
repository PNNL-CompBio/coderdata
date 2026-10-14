#!/usr/bin/env python3
"""
02-omics-cnf.py: generate cnf_transcriptomics.csv, cnf_proteomics.csv,
                 and cnf_phosphoproteomics.csv.

For the cNF dataset:

* Transcriptomics: raw gene-level TPM from per-cohort salmon.merged.gene_tpm.tsv
  files (wide matrix: rows = genes, columns = specimens). Default cohorts:
    syn66352931  — cNF_Cohort_1 (NF0017–NF0023)
    syn70765053  — cNF_Cohort_2 (NF0025, NF0027, NF0035)
  Override via env var CNF_RNA_COHORT_SYN_IDS (comma-separated Synapse IDs).

* Proteomics: batch-corrected global proteomics from syn74815895
  (correctedAbundance column). Override via env var CNF_GLOBAL_PROT_SYN_ID.

* Phosphoproteomics: batch-corrected phospho-proteomics from syn70078415
  (correctedAbundance column). Mapped to phosphosite_ids via a phosphosites.csv
  reference file. Override via env var CNF_PHOSPHO_SYN_ID.

* Protocol optimization samples are excluded from RNA via DROP_SAMPLE_SUBSTRINGS.

Specimens are canonicalized via cnf_utils so that values like
"NF0018.T1.organoid" resolve to the canonical IDs in cnf_samples.csv.
Source rows whose specimen doesn't match any entry in cnf_samples.csv are
dropped silently via the inner join to sample_map.

Usage:
    python 02-omics-cnf.py <genes.csv> <cnf_samples.csv> <phosphosites.csv>
                          [--out_rna  /tmp/cnf_transcriptomics.csv]
                          [--out_prot /tmp/cnf_proteomics.csv]
                          [--out_phospho /tmp/cnf_phosphoproteomics.csv]
"""

import argparse
import logging
import os
import re
import sys

import pandas as pd
import synapseclient

from cnf_utils import canonicalize_specimen_column


# ----------------------------------------------------------------------------
# Defaults / Synapse IDs
# ----------------------------------------------------------------------------
_rna_env = os.environ.get("CNF_RNA_COHORT_SYN_IDS", "")
RNA_TPM_COHORT_SYN_IDS: list[str] = (
    [s.strip() for s in _rna_env.split(",") if s.strip()]
    if _rna_env
    else ["syn66352931", "syn70765053"]
)

GLOBAL_PROT_SYN_ID  = os.environ.get("CNF_GLOBAL_PROT_SYN_ID",  "syn74815895")
PHOSPHO_SYN_ID      = os.environ.get("CNF_PHOSPHO_SYN_ID",      "syn70078415")

# Protocol-optimization samples to drop from RNA.
DROP_SAMPLE_SUBSTRINGS: list[str] = [
    "cNF_organoid_DIA_G_02_11Feb25",
    "cNF_organoid_DIA_G_05_11Feb25",
    "cNF_organoid_DIA_G_06_11Feb25",
    "cNF_organoid_DIA_P_02_29Jan25",
    "cNF_organoid_DIA_P_05_11Feb25",
    "cNF_organoid_DIA_P_06_11Feb25",
]

STUDY               = "cnf"
DEFAULT_OUT_RNA     = "/tmp/cnf_transcriptomics.csv"
DEFAULT_OUT_PROT    = "/tmp/cnf_proteomics.csv"
DEFAULT_OUT_PHOSPHO = "/tmp/cnf_phosphoproteomics.csv"


def configure_logging():
    logging.basicConfig(
        format="[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )


# ----------------------------------------------------------------------------
# Reference data loaders
# ----------------------------------------------------------------------------
def load_gene_map(genes_path: str) -> pd.DataFrame:
    """Load genes.csv. Required columns: gene_symbol, entrez_id."""
    df = pd.read_csv(genes_path)
    needed = {"gene_symbol", "entrez_id"}
    missing = needed - set(df.columns)
    if missing:
        raise KeyError(f"genes file missing columns: {missing}")
    df = df.dropna(subset=["entrez_id", "gene_symbol"]).copy()
    df["entrez_id"] = df["entrez_id"].astype(int)
    df = df.drop_duplicates(subset="gene_symbol")
    logging.info("Loaded %d gene_symbol → entrez_id mappings", len(df))
    return df[["gene_symbol", "entrez_id"]]


def load_sample_map(samples_path: str) -> pd.DataFrame:
    """Load cnf_samples.csv; return (other_id, improve_sample_id) frame."""
    df = pd.read_csv(samples_path)
    df = df[["other_id", "improve_sample_id"]].drop_duplicates()
    df["improve_sample_id"] = df["improve_sample_id"].astype(int)
    logging.info("Loaded %d sample mappings", len(df))
    return df


def load_phosphosite_map(phosphosites_path: str) -> pd.DataFrame:
    """Load phosphosites.csv; return (other_id, phosphosite_id) frame."""
    df = pd.read_csv(phosphosites_path)
    needed = {"other_id", "phosphosite_id"}
    missing = needed - set(df.columns)
    if missing:
        raise KeyError(f"phosphosites file missing columns: {missing}")
    df = df[["other_id", "phosphosite_id"]].drop_duplicates()
    df["phosphosite_id"] = df["phosphosite_id"].astype(int)
    logging.info("Loaded %d phosphosite mappings", len(df))
    return df


def fetch_synapse_table(syn, syn_id: str) -> pd.DataFrame:
    """Download a Synapse file and read it as CSV/TSV."""
    if not syn_id:
        raise ValueError("Synapse ID is empty")
    logging.info("Fetching %s from Synapse", syn_id)
    entity = syn.get(syn_id)
    path = entity.path
    if path.endswith(".tsv") or path.endswith(".tsv.gz"):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


# ----------------------------------------------------------------------------
# Shared helpers
# ----------------------------------------------------------------------------
def _col_stem(col: str) -> str:
    """Return the basename without extension for a path-like column name."""
    base = re.sub(r"^.*[/\\]", "", col)
    return re.sub(r"\.[^.]+$", "", base)


def _drop_protocol_cols(cols: list[str]) -> list[str]:
    """Remove columns whose stem matches any DROP_SAMPLE_SUBSTRINGS entry."""
    kept = []
    for c in cols:
        stem = _col_stem(c)
        if any(sub in stem for sub in DROP_SAMPLE_SUBSTRINGS):
            logging.info("  Dropping protocol optimization sample: %s", stem)
        else:
            kept.append(c)
    return kept


def _normalize_col(df: pd.DataFrame, candidates: list[str], target: str) -> pd.DataFrame:
    """Rename the first matching candidate column to target."""
    for c in candidates:
        if c in df.columns:
            return df.rename(columns={c: target})
    raise KeyError(
        f"None of {candidates} found in columns {list(df.columns)[:20]}"
    )


# ----------------------------------------------------------------------------
# Transcriptomics
# ----------------------------------------------------------------------------
def fetch_wide_tpm_cohort(syn, syn_id: str) -> pd.DataFrame:
    """Fetch a wide salmon.merged.gene_tpm.tsv and return long-format DataFrame.

    nf-core salmon output: gene_id (Ensembl), gene_name (symbol), then one
    column per specimen. Returns (gene_symbol, specimen, transcriptomics, source).
    """
    logging.info("Fetching cohort TPM matrix %s", syn_id)
    entity = syn.get(syn_id)
    path = entity.path
    df = pd.read_csv(path, sep="\t")

    gene_col = None
    for cname in ["gene_name", "gene_symbol", "Symbol", "GeneSymbol", "gene"]:
        if cname in df.columns:
            gene_col = cname
            break
    if gene_col is None:
        gene_col = df.columns[0]
        logging.warning(
            "%s: no gene_name column found; using %s — most rows may fail the gene_map join",
            syn_id, gene_col,
        )

    meta_cols = {"gene_id", "gene_name", "gene_symbol", "Symbol",
                 "GeneSymbol", "gene", "transcript_id"}
    sample_cols = [c for c in df.columns if c not in meta_cols]
    sample_cols = _drop_protocol_cols(sample_cols)

    if not sample_cols:
        raise ValueError(
            f"{syn_id}: no specimen columns found after excluding metadata "
            f"and protocol optimization samples"
        )

    long = (
        df[[gene_col] + sample_cols]
        .rename(columns={gene_col: "gene_symbol"})
        .melt(id_vars="gene_symbol", var_name="specimen", value_name="transcriptomics")
    )
    long["source"] = "Synapse"
    logging.info(
        "  %s: %d genes × %d specimens → %d long rows",
        syn_id, df[gene_col].nunique(), len(sample_cols), len(long),
    )
    return long


def build_transcriptomics(
    syn, gene_map: pd.DataFrame, sample_map: pd.DataFrame
) -> pd.DataFrame:
    """Build the long-format transcriptomics table from per-cohort TPM matrices."""
    if not RNA_TPM_COHORT_SYN_IDS:
        logging.error("RNA_TPM_COHORT_SYN_IDS is empty; skipping transcriptomics.")
        return pd.DataFrame(
            columns=["entrez_id", "improve_sample_id", "source", "study", "transcriptomics"]
        )

    cohort_frames = []
    for syn_id in RNA_TPM_COHORT_SYN_IDS:
        try:
            cohort_frames.append(fetch_wide_tpm_cohort(syn, syn_id))
        except Exception as exc:
            logging.warning("Skipping RNA cohort %s: %s", syn_id, exc)

    if not cohort_frames:
        logging.error("All RNA cohort fetches failed; returning empty transcriptomics.")
        return pd.DataFrame(
            columns=["entrez_id", "improve_sample_id", "source", "study", "transcriptomics"]
        )

    raw = pd.concat(cohort_frames, ignore_index=True)
    logging.info("Combined %d RNA cohorts → %d total rows", len(cohort_frames), len(raw))

    n_before = len(raw)
    raw = canonicalize_specimen_column(raw, "specimen")
    logging.info(
        "Transcriptomics: canonicalized %d → %d rows (%d dropped as unrecognized specimens)",
        n_before, len(raw), n_before - len(raw),
    )

    df = (
        raw.merge(gene_map, on="gene_symbol", how="inner")
        .merge(sample_map, left_on="specimen", right_on="other_id", how="inner")
    )
    df["study"] = STUDY
    df["transcriptomics"] = pd.to_numeric(df["transcriptomics"], errors="coerce")
    df = df.dropna(subset=["transcriptomics"])
    df = df[["entrez_id", "improve_sample_id", "source", "study", "transcriptomics"]]

    logging.info(
        "Built transcriptomics table: %d rows (%d samples × %d genes)",
        len(df), df["improve_sample_id"].nunique(), df["entrez_id"].nunique(),
    )
    return df


# ----------------------------------------------------------------------------
# Proteomics
# ----------------------------------------------------------------------------
def build_proteomics(
    syn, gene_map: pd.DataFrame, sample_map: pd.DataFrame
) -> pd.DataFrame:
    """Build the long-format proteomics table from the batch-corrected Synapse file."""
    raw = fetch_synapse_table(syn, GLOBAL_PROT_SYN_ID)
    raw = _normalize_col(raw,
        ["Specimen", "specimenID", "specimen_id", "specimen", "specimen_norm", "sample_id"],
        "specimen")
    raw = _normalize_col(raw,
        ["gene_symbol", "Gene", "GeneSymbol", "Symbol", "gene", "feature_id"],
        "gene_symbol")

    n_before = len(raw)
    raw = canonicalize_specimen_column(raw, "specimen")
    logging.info(
        "Proteomics: canonicalized %d → %d rows (%d dropped as unrecognized specimens)",
        n_before, len(raw), n_before - len(raw),
    )

    raw = _normalize_col(raw,
        ["correctedAbundance", "log_ratio", "logRatio", "proteomics", "value"],
        "proteomics")

    df = (
        raw.merge(gene_map, on="gene_symbol", how="inner")
        .merge(sample_map, left_on="specimen", right_on="other_id", how="inner")
    )
    df["source"] = "Synapse"
    df["study"]  = STUDY
    df["proteomics"] = pd.to_numeric(df["proteomics"], errors="coerce")
    df = df.dropna(subset=["proteomics"])
    df = df[["entrez_id", "improve_sample_id", "source", "study", "proteomics"]]

    logging.info(
        "Built proteomics table: %d rows (%d samples × %d genes)",
        len(df), df["improve_sample_id"].nunique(), df["entrez_id"].nunique(),
    )
    return df


# ----------------------------------------------------------------------------
# Phosphoproteomics
# ----------------------------------------------------------------------------
def build_phosphoproteomics(
    syn, phosphosite_map: pd.DataFrame, sample_map: pd.DataFrame
) -> pd.DataFrame:
    """Build the long-format phosphoproteomics table from the batch-corrected Synapse file."""
    raw = fetch_synapse_table(syn, PHOSPHO_SYN_ID)
    raw = _normalize_col(raw,
        ["Specimen", "specimenID", "specimen_id", "specimen"],
        "specimen")
    raw = _normalize_col(raw,
        ["site", "Site", "phosphosite", "feature_id"],
        "site")

    n_before = len(raw)
    raw = canonicalize_specimen_column(raw, "specimen")
    logging.info(
        "Phosphoproteomics: canonicalized %d → %d rows (%d dropped as unrecognized specimens)",
        n_before, len(raw), n_before - len(raw),
    )

    raw = _normalize_col(raw,
        ["correctedAbundance", "log_ratio", "logRatio", "phosphoproteomics", "value"],
        "phosphoproteomics")

    df = (
        raw.merge(phosphosite_map, left_on="site", right_on="other_id", how="inner")
        .merge(sample_map, left_on="specimen", right_on="other_id", how="inner")
    )
    df["source"] = "Synapse"
    df["study"]  = STUDY
    df["phosphoproteomics"] = pd.to_numeric(df["phosphoproteomics"], errors="coerce")
    df = df.dropna(subset=["phosphoproteomics"])
    df = df[["phosphosite_id", "improve_sample_id", "source", "study", "phosphoproteomics"]]

    logging.info(
        "Built phosphoproteomics table: %d rows (%d samples × %d sites)",
        len(df), df["improve_sample_id"].nunique(), df["phosphosite_id"].nunique(),
    )
    return df


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("genes",       help="Path to genes.csv (gene_symbol → entrez_id)")
    parser.add_argument("samples",     help="Path to cnf_samples.csv from build_samples")
    parser.add_argument("phosphosites", help="Path to phosphosites.csv (other_id → phosphosite_id)")
    parser.add_argument("--out_rna",    default=DEFAULT_OUT_RNA)
    parser.add_argument("--out_prot",   default=DEFAULT_OUT_PROT)
    parser.add_argument("--out_phospho", default=DEFAULT_OUT_PHOSPHO)
    args = parser.parse_args()

    configure_logging()

    gene_map         = load_gene_map(args.genes)
    sample_map       = load_sample_map(args.samples)
    phosphosite_map  = load_phosphosite_map(args.phosphosites)

    syn = synapseclient.Synapse()
    syn.login()

    rna    = build_transcriptomics(syn, gene_map, sample_map)
    prot   = build_proteomics(syn, gene_map, sample_map)
    phospho = build_phosphoproteomics(syn, phosphosite_map, sample_map)

    for path in (args.out_rna, args.out_prot, args.out_phospho):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    rna.to_csv(args.out_rna,     index=False)
    prot.to_csv(args.out_prot,   index=False)
    phospho.to_csv(args.out_phospho, index=False)
    logging.info("Wrote %s (%d rows)", args.out_rna,     len(rna))
    logging.info("Wrote %s (%d rows)", args.out_prot,    len(prot))
    logging.info("Wrote %s (%d rows)", args.out_phospho, len(phospho))


if __name__ == "__main__":
    sys.exit(main() or 0)
