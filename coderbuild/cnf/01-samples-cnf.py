#!/usr/bin/env python3
"""
01-samples-cnf.py: generate cnf_samples.csv for the cNF organoid drug screen.

Specimens are gathered from four sources:

  1. Drug-screen file index on Synapse (syn51301431, dataType='drug screen')
  2. Global proteomics matrix (env CNF_GLOBAL_PROT_SYN_ID, default syn74815895)
  3. RNA discovery file (env CNF_RNA_DISCOVERY_SYN_ID, default syn71333780).
     This is the corrected-abundance file containing every RNA condition
     processed for the project (organoid, tumor tissue, drug-treated
     organoid, skin), used here for sample discovery only — NOT to be
     confused with CNF_RNA_TPM_SYN_ID, which build_omics.py reads for the
     actual TPM data and may be filtered to organoids only.
  4. Normal Skin Synapse folder (env CNF_NORMAL_SKIN_SYN_ID, default syn74284682)

Each specimen is classified as one of:
  - "organoid"         → model_type = "patient derived organoid"
  - "treated_organoid" → model_type = "patient derived organoid"
                          (drug-treated organoid; canonical ID preserves the
                           treatment label, e.g. NF0021_T1_Onalespid_1uM)
  - "tumor_tissue"     → model_type = "tumor"
  - "normal_skin"      → model_type = "normal_tissue"

Tissue and organoid samples for the same tumor are kept as separate rows with
distinct improve_sample_ids — they have different biology (fresh tumor vs.
cultured organoid) and the schema model_types differ. Treated and untreated
organoids of the same tumor are also separate rows because their canonical
IDs differ (treated samples preserve the treatment label).

Cohort assignment is intentionally omitted. The upstream omics data is
batch-corrected and per-modality cohort membership disagrees for some
specimens, so a single per-sample cohort label cannot be assigned cleanly.
The `other_id_source` field carries a uniform project label.

Sample IDs continue from the previous samples file: max(improve_sample_id) + 1.

Usage:
    python 01-samples-cnf.py <previous_samples.csv> [--output /tmp/cnf_samples.csv]
"""

import argparse
import logging
import os
import re
import sys

import pandas as pd
import synapseclient

from cnf_utils import (
    MODEL_TYPE_MAP,
    classify_specimen,
    patient_from_specimen,
)


# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------
DRUG_SCREEN_TABLE       = "syn51301431"
GLOBAL_PROT_SYN_ID      = os.environ.get("CNF_GLOBAL_PROT_SYN_ID", "syn74815895")
# RNA file used for SPECIMEN DISCOVERY only — needs to contain all
# conditions (organoid + tissue + treated + skin). The corrected-abundance
# file at syn71333780 fits this. Also used by build_omics.py as the TPM
# data source (the column is named correctedAbundance but stores TPM).
RNA_DISCOVERY_SYN_ID    = os.environ.get("CNF_RNA_DISCOVERY_SYN_ID", "syn71333780")
NORMAL_SKIN_SYN_ID      = os.environ.get("CNF_NORMAL_SKIN_SYN_ID", "syn74284682")

DEFAULT_OUTPUT = "/tmp/cnf_samples.csv"

CANCER_TYPE     = "Cutaneous Neurofibroma"
SPECIES         = "Homo sapiens (Human)"
OTHER_ID_SOURCE = "NF1_cNF_Project"


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def configure_logging():
    logging.basicConfig(
        format="[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )


def get_starting_id(prev_samples_path: str) -> int:
    """Read previous samples file and return max(improve_sample_id) + 1."""
    if not os.path.exists(prev_samples_path):
        logging.warning("Previous samples file %s not found; starting at 1.",
                         prev_samples_path)
        return 1

    df = pd.read_csv(prev_samples_path)
    if "improve_sample_id" not in df.columns or df.empty:
        logging.warning("No improve_sample_id column or empty file; starting at 1.")
        return 1

    max_id = int(pd.to_numeric(df["improve_sample_id"], errors="coerce").max())
    logging.info("Previous samples max improve_sample_id = %d", max_id)
    return max_id + 1


def merge_specimen_dicts(*dicts: dict[str, str]) -> dict[str, str]:
    """Union per-source dicts; warn on type disagreements (first-seen wins)."""
    merged: dict[str, str] = {}
    for d in dicts:
        for spec, stype in d.items():
            if spec in merged and merged[spec] != stype:
                logging.warning(
                    "Type conflict for %s: %s (kept) vs %s (ignored)",
                    spec, merged[spec], stype,
                )
                continue
            merged.setdefault(spec, stype)
    return merged


# ----------------------------------------------------------------------------
# Source-specific specimen extraction
# ----------------------------------------------------------------------------
def specimens_from_drug_screen(syn) -> dict[str, str]:
    """Pull canonical_id → sample_type for drug-screen specimens."""
    query = (
        "select id, individualID, specimenID "
        f"from {DRUG_SCREEN_TABLE} "
        "where dataType='drug screen'"
    )
    logging.info("Querying %s for drug-screen specimens", DRUG_SCREEN_TABLE)
    table = syn.tableQuery(query).asDataFrame()
    raw = table["specimenID"].dropna().astype(str).unique().tolist()
    out: dict[str, str] = {}
    for r in raw:
        result = classify_specimen(r)
        if result:
            cid, stype = result
            out[cid] = stype
    by_type = _count_by_type(out)
    logging.info("  Drug screen: %s", by_type)
    return out


def _count_by_type(d: dict[str, str]) -> dict[str, int]:
    out: dict[str, int] = {}
    for stype in d.values():
        out[stype] = out.get(stype, 0) + 1
    return out


def _detect_specimen_column(df: pd.DataFrame) -> str:
    """Return the name of the specimen column in an omics file."""
    candidates = ["Specimen", "specimen", "specimen_norm", "sample_id",
                  "specimenID", "Sample"]
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(
        f"No specimen column among {candidates} in columns {list(df.columns)[:15]}"
    )


def _read_synapse_table(syn, syn_id: str) -> pd.DataFrame:
    """Fetch a Synapse file and read it as CSV/TSV."""
    logging.info("  Fetching %s", syn_id)
    entity = syn.get(syn_id)
    path = entity.path
    if path.endswith(".tsv") or path.endswith(".tsv.gz"):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def specimens_from_proteomics(syn) -> dict[str, str]:
    """Pull canonical_id → sample_type from the global proteomics matrix."""
    if not GLOBAL_PROT_SYN_ID:
        logging.warning("CNF_GLOBAL_PROT_SYN_ID not set; skipping proteomics")
        return {}

    logging.info("Reading proteomics specimens")
    try:
        df = _read_synapse_table(syn, GLOBAL_PROT_SYN_ID)
        col = _detect_specimen_column(df)
        raw = df[col].dropna().astype(str).unique().tolist()
        out: dict[str, str] = {}
        for r in raw:
            result = classify_specimen(r)
            if result:
                cid, stype = result
                out[cid] = stype
        logging.info("  Proteomics: %s", _count_by_type(out))
        return out
    except Exception as exc:
        logging.warning("  Proteomics read failed: %s", exc)
        return {}


def specimens_from_rna_discovery(syn) -> dict[str, str]:
    """Discover RNA specimens from a comprehensive long-format RNA file.

    Reads CNF_RNA_DISCOVERY_SYN_ID (default syn71333780) — the corrected-
    abundance RNA file that contains every condition processed for the
    project: organoid + tumor tissue + drug-treated organoid + skin.

    This is intentionally a different file from the one read by
    build_omics.py (CNF_RNA_TPM_SYN_ID, which is the gene-level TPM matrix
    and may be filtered to organoids only). Sample discovery and data
    ingestion have different requirements: discovery wants the maximum
    sample list; ingestion wants TPM values.

    Specimens are passed through classify_specimen, so drug-treated
    samples like NF0021_T1_Onalespid.1uM get caught regardless of their
    condition_norm value.
    """
    if not RNA_DISCOVERY_SYN_ID:
        logging.warning("CNF_RNA_DISCOVERY_SYN_ID not set; skipping RNA discovery")
        return {}

    logging.info("Reading RNA discovery file (%s)", RNA_DISCOVERY_SYN_ID)
    try:
        df = _read_synapse_table(syn, RNA_DISCOVERY_SYN_ID)
        col = _detect_specimen_column(df)

        raw = df[col].dropna().astype(str).unique().tolist()
        out: dict[str, str] = {}
        for r in raw:
            result = classify_specimen(r)
            if result:
                cid, stype = result
                out[cid] = stype
        logging.info("  RNA discovery: %s", _count_by_type(out))
        return out
    except Exception as exc:
        logging.warning("  RNA discovery read failed: %s", exc)
        return {}


def specimens_from_normal_skin_tree(syn) -> dict[str, str]:
    """Discover skin specimens by listing the Normal Skin Synapse folder.

    The "Normal Skin" tree contains one immediate child per patient (folders
    named NF0017, NF0018, ...) plus, in some cases, top-level skin FASTQ
    files like NF0035-SKIN-1-04-16-2025_*.fastq.gz. Both forms encode the
    patient ID at the start, so we match `^NF\\d{4}` against each child name
    and emit one canonical NFxxxx_skin entry per unique patient.

    This is intentionally one API call deep — it avoids walking the whole
    WGS tree (which has hundreds of empty lane-stub folders). For most
    cases this captures every patient with a matched normal.
    """
    if not NORMAL_SKIN_SYN_ID:
        logging.warning("CNF_NORMAL_SKIN_SYN_ID not set; skipping Normal Skin tree")
        return {}

    logging.info("Listing Normal Skin tree (%s)", NORMAL_SKIN_SYN_ID)
    out: dict[str, str] = {}
    try:
        for child in syn.getChildren(NORMAL_SKIN_SYN_ID):
            name = str(child.get("name", ""))
            m = re.match(r"^(NF\d{4})", name)
            if m:
                patient = m.group(1)
                out[f"{patient}_skin"] = "normal_skin"
        logging.info("  Normal Skin tree: %s", _count_by_type(out))
    except Exception as exc:
        logging.warning("  Normal Skin tree listing failed: %s", exc)
    return out


# ----------------------------------------------------------------------------
# Sample table construction
# ----------------------------------------------------------------------------
def build_samples_table(specimens: dict[str, str], start_id: int) -> pd.DataFrame:
    """Build the schema-compliant Sample table from {canonical_id: sample_type}.

    Sort order: organoid → treated_organoid → tumor_tissue → normal_skin
    (each by ID). Drug-modelable untreated organoids get the lowest
    improve_sample_ids — useful for downstream filtering.
    """
    sort_priority = {
        "organoid":         0,
        "treated_organoid": 1,
        "tumor_tissue":     2,
        "normal_skin":      3,
    }
    sorted_specimens = sorted(
        specimens.items(),
        key=lambda kv: (sort_priority.get(kv[1], 99), kv[0]),
    )

    rows = []
    next_id = start_id
    for specimen, stype in sorted_specimens:
        patient = patient_from_specimen(specimen)
        rows.append({
            "improve_sample_id": next_id,
            "other_id":          specimen,
            "other_id_source":   OTHER_ID_SOURCE,
            "common_name":       specimen,
            "cancer_type":       CANCER_TYPE,
            "other_names":       patient,
            "species":           SPECIES,
            "model_type":        MODEL_TYPE_MAP[stype],
        })
        next_id += 1

    df = pd.DataFrame(rows)
    logging.info("Built %d sample rows; ID range %d–%d",
                 len(df), start_id, next_id - 1)
    return df


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "prev_samples",
        help="Path to previous samples CSV (used to find max improve_sample_id)",
    )
    parser.add_argument("--output", default=DEFAULT_OUTPUT,
                         help=f"Output CSV path (default: {DEFAULT_OUTPUT})")
    args = parser.parse_args()

    configure_logging()

    syn = synapseclient.Synapse()
    syn.login()

    # Pull specimens from each source as canonical_id → sample_type dicts
    drug_d = specimens_from_drug_screen(syn)
    prot_d = specimens_from_proteomics(syn)
    rna_d  = specimens_from_rna_discovery(syn)
    skin_d = specimens_from_normal_skin_tree(syn)

    all_specimens = merge_specimen_dicts(drug_d, prot_d, rna_d, skin_d)
    logging.info(
        "Union: %d total specimens — %s",
        len(all_specimens), _count_by_type(all_specimens),
    )

    # Per-source uniqueness diagnostics
    drug_set, prot_set, rna_set, skin_set = (
        set(drug_d), set(prot_d), set(rna_d), set(skin_d),
    )
    only_drug = drug_set - prot_set - rna_set - skin_set
    only_prot = prot_set - drug_set - rna_set - skin_set
    only_rna  = rna_set  - drug_set - prot_set - skin_set
    only_skin = skin_set - drug_set - prot_set - rna_set
    if only_drug:
        logging.info("  Only in drug screen: %s", sorted(only_drug))
    if only_prot:
        logging.info("  Only in proteomics: %s", sorted(only_prot))
    if only_rna:
        logging.info("  Only in RNA: %s", sorted(only_rna))
    if only_skin:
        logging.info("  Only in Normal Skin tree: %s", sorted(only_skin))

    if not all_specimens:
        logging.error("No specimens found from any source; aborting")
        sys.exit(1)

    start_id = get_starting_id(args.prev_samples)
    samples = build_samples_table(all_specimens, start_id)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    samples.to_csv(args.output, index=False)
    logging.info("Wrote %s (%d rows)", args.output, len(samples))


if __name__ == "__main__":
    sys.exit(main() or 0)
