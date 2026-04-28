#!/usr/bin/env python3
"""
02-omics-cnf.py — generate cnf_transcriptomics.csv and cnf_proteomics.csv.

For the cNF dataset:

* Transcriptomics: gene-level TPM (set CNF_RNA_TPM_SYN_ID to the Synapse
  entity holding the gene-level TPM matrix or long table).
* Proteomics: log-ratio global proteomics from syn70078416 (correctedAbundance).
* Phospho-proteomics is omitted by design — the schema has no slot for it.

Genes come in by HGNC symbol and are mapped to entrez_id via the genes file
passed as the first argument. Specimens come in by NFxxxx_Tx and are mapped
to integer improve_sample_id via the cnf_samples.csv passed as the second
argument.

Only specimens with model_type "patient derived organoid" (i.e. the entries
build_samples.py wrote) are kept. Source rows whose specimen ID does not
match any entry in cnf_samples.csv are dropped with a warning.

Usage:
    python build_omics.py <genes.csv> <cnf_samples.csv>
"""

import argparse
import logging
import os
import sys

import pandas as pd
import synapseclient


# ----------------------------------------------------------------------------
# Defaults / Synapse IDs
# ----------------------------------------------------------------------------
RNA_TPM_SYN_ID_DEFAULT  = os.environ.get("CNF_RNA_TPM_SYN_ID", "")
GLOBAL_PROT_SYN_ID      = os.environ.get("CNF_GLOBAL_PROT_SYN_ID", "syn70078416")
STUDY                   = "cnf"


def configure_logging():
    logging.basicConfig(
        format="[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )


def load_gene_map(genes_path: str) -> pd.DataFrame:
    """Load the genes.csv file — expected columns include gene_symbol and entrez_id."""
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
    """Load the cNF samples file; return other_id → improve_sample_id mapping."""
    df = pd.read_csv(samples_path)
    df = df[["other_id", "improve_sample_id"]].drop_duplicates()
    df["improve_sample_id"] = df["improve_sample_id"].astype(int)
    logging.info("Loaded %d sample mappings", len(df))
    return df


def fetch_long_table(syn, syn_id: str) -> pd.DataFrame:
    """Pull a Synapse file that's stored as long-format CSV/TSV."""
    if not syn_id:
        raise ValueError("Synapse ID is empty; cannot fetch table")
    logging.info("Fetching %s from Synapse", syn_id)
    entity = syn.get(syn_id)
    path = entity.path
    if path.endswith(".tsv") or path.endswith(".tsv.gz"):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def filter_to_organoid_specimens(df: pd.DataFrame, specimen_col: str) -> pd.DataFrame:
    """Drop rows whose specimen ID has _tissue / _skin suffix; keep _organoid and bare."""
    if specimen_col not in df.columns:
        return df
    excluded = ["_tissue", "_skin"]
    mask = ~df[specimen_col].astype(str).str.contains("|".join(excluded), na=False)
    n_drop = (~mask).sum()
    if n_drop:
        logging.info("Dropped %d rows with non-organoid specimens (tissue/skin)", n_drop)
    df = df[mask].copy()
    df[specimen_col] = (
        df[specimen_col].astype(str).str.replace("_organoid", "", regex=False)
    )
    return df


def normalize_specimen_col(df: pd.DataFrame) -> pd.DataFrame:
    """Try common specimen column names; standardize to 'specimen'."""
    candidates = ["Specimen", "specimenID", "specimen_id", "specimen", "sample"]
    for c in candidates:
        if c in df.columns:
            df = df.rename(columns={c: "specimen"})
            return df
    raise KeyError(
        f"Could not find a specimen column among {candidates} "
        f"in columns {list(df.columns)[:20]}"
    )


def normalize_gene_col(df: pd.DataFrame) -> pd.DataFrame:
    """Try common gene-symbol column names; standardize to 'gene_symbol'."""
    candidates = ["gene_symbol", "Gene", "GeneSymbol", "Symbol", "gene"]
    for c in candidates:
        if c in df.columns:
            df = df.rename(columns={c: "gene_symbol"})
            return df
    raise KeyError(
        f"Could not find a gene column among {candidates} "
        f"in columns {list(df.columns)[:20]}"
    )


def build_transcriptomics(
    syn, gene_map: pd.DataFrame, sample_map: pd.DataFrame
) -> pd.DataFrame:
    """Build the long-format transcriptomics table."""
    if not RNA_TPM_SYN_ID_DEFAULT:
        logging.error(
            "CNF_RNA_TPM_SYN_ID not set; skipping transcriptomics. "
            "Set this env var to the Synapse ID of the gene-level TPM table."
        )
        return pd.DataFrame(
            columns=["entrez_id", "improve_sample_id", "source",
                     "study", "transcriptomics"]
        )

    raw = fetch_long_table(syn, RNA_TPM_SYN_ID_DEFAULT)
    raw = normalize_specimen_col(raw)
    raw = normalize_gene_col(raw)
    raw = filter_to_organoid_specimens(raw, "specimen")

    # Find the TPM value column
    value_candidates = ["TPM", "tpm", "transcriptomics", "value"]
    val_col = next((c for c in value_candidates if c in raw.columns), None)
    if val_col is None:
        raise KeyError(
            f"No TPM value column among {value_candidates} "
            f"in columns {list(raw.columns)[:20]}"
        )

    df = (
        raw.rename(columns={val_col: "transcriptomics"})
        .merge(gene_map, on="gene_symbol", how="inner")
        .merge(sample_map, left_on="specimen", right_on="other_id", how="inner")
    )

    df["source"] = f"Synapse:{RNA_TPM_SYN_ID_DEFAULT}"
    df["study"] = STUDY
    df = df[["entrez_id", "improve_sample_id", "source",
             "study", "transcriptomics"]]
    df["transcriptomics"] = pd.to_numeric(df["transcriptomics"], errors="coerce")
    df = df.dropna(subset=["transcriptomics"])
    logging.info("Built transcriptomics table: %d rows", len(df))
    return df


def build_proteomics(
    syn, gene_map: pd.DataFrame, sample_map: pd.DataFrame
) -> pd.DataFrame:
    """Build the long-format proteomics table from global proteomics."""
    raw = fetch_long_table(syn, GLOBAL_PROT_SYN_ID)
    raw = normalize_specimen_col(raw)
    raw = normalize_gene_col(raw)
    raw = filter_to_organoid_specimens(raw, "specimen")

    # correctedAbundance is already log ratio; that's what schema expects.
    val_candidates = ["correctedAbundance", "log_ratio", "logRatio",
                      "proteomics", "value"]
    val_col = next((c for c in val_candidates if c in raw.columns), None)
    if val_col is None:
        raise KeyError(
            f"No proteomics value column among {val_candidates} "
            f"in columns {list(raw.columns)[:20]}"
        )

    df = (
        raw.rename(columns={val_col: "proteomics"})
        .merge(gene_map, on="gene_symbol", how="inner")
        .merge(sample_map, left_on="specimen", right_on="other_id", how="inner")
    )

    df["source"] = f"Synapse:{GLOBAL_PROT_SYN_ID}"
    df["study"] = STUDY
    df = df[["entrez_id", "improve_sample_id", "source", "study", "proteomics"]]
    df["proteomics"] = pd.to_numeric(df["proteomics"], errors="coerce")
    df = df.dropna(subset=["proteomics"])
    logging.info("Built proteomics table: %d rows", len(df))
    return df


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("genes", help="Path to genes.csv (gene_symbol → entrez_id)")
    parser.add_argument("samples", help="Path to cnf_samples.csv from build_samples")
    parser.add_argument("--out_rna", default="/tmp/cnf_transcriptomics.csv")
    parser.add_argument("--out_prot", default="/tmp/cnf_proteomics.csv")
    args = parser.parse_args()

    configure_logging()

    gene_map   = load_gene_map(args.genes)
    sample_map = load_sample_map(args.samples)

    syn = synapseclient.Synapse()
    syn.login()

    rna  = build_transcriptomics(syn, gene_map, sample_map)
    prot = build_proteomics(syn, gene_map, sample_map)

    rna.to_csv(args.out_rna,  index=False)
    prot.to_csv(args.out_prot, index=False)
    logging.info("Wrote %s (%d rows)", args.out_rna, len(rna))
    logging.info("Wrote %s (%d rows)", args.out_prot, len(prot))


if __name__ == "__main__":
    sys.exit(main() or 0)