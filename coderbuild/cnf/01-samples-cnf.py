#!/usr/bin/env python3
"""
01-samples-cnf.py — generate cnf_samples.csv for the cNF (cutaneous neurofibroma)
organoid drug screen.

Reads the most recent samples file (passed as the first argument), finds the
maximum existing improve_sample_id, and mints new integer IDs starting at
max+1 for every cNF specimen with a drug screen on Synapse.

Cohort assignment (Cohort 1 vs Cohort 2) is encoded in `other_id_source`.
The `study` column is always "cnf". Patient IDs (e.g. NF0019) are stored in
`other_names`, while the per-tumor specimen ID (e.g. NF0019_T3) is the
`other_id` and `common_name`.

Usage:
    python build_samples.py <previous_samples.csv> [--output cnf_samples.csv]
"""

import argparse
import os
import sys
import logging
import pandas as pd
import synapseclient


# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------
DRUG_SCREEN_TABLE = "syn51301431"   # Synapse table with drug-screen file metadata
DEFAULT_OUTPUT    = "cnf_samples.csv"

# Patient cohort lookup. Encoded in other_id_source; documented in README.
COHORT_LOOKUP = {
    "NF0017": "cNF_Cohort_1",
    "NF0018": "cNF_Cohort_1",
    "NF0019": "cNF_Cohort_1",
    "NF0020": "cNF_Cohort_1",
    "NF0021": "cNF_Cohort_1",
    "NF0022": "cNF_Cohort_1",
    "NF0023": "cNF_Cohort_1",
    "NF0025": "cNF_Cohort_2",
    "NF0027": "cNF_Cohort_2",
    "NF0035": "cNF_Cohort_2",
}

CANCER_TYPE = "Cutaneous Neurofibroma"
SPECIES     = "Homo sapiens"
MODEL_TYPE  = "patient derived organoid"


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


def query_drug_screen_specimens(syn) -> pd.DataFrame:
    """Query Synapse for cNF drug-screen specimens; return unique specimen list."""
    query = (
        "select id, individualID, specimenID "
        f"from {DRUG_SCREEN_TABLE} "
        "where dataType='drug screen'"
    )
    logging.info("Querying %s for drug-screen specimens", DRUG_SCREEN_TABLE)
    table = syn.tableQuery(query).asDataFrame()

    # Each (individualID, specimenID) pair is one biological tumor; one drug-screen
    # file may have multiple plate replicates. Deduplicate.
    specimens = (
        table[["individualID", "specimenID"]]
        .dropna()
        .drop_duplicates()
        .reset_index(drop=True)
    )
    logging.info("Found %d unique cNF specimens", len(specimens))
    return specimens


def assign_cohort(individual_id: str) -> str:
    """Look up cohort from patient ID; warn if missing."""
    cohort = COHORT_LOOKUP.get(individual_id)
    if cohort is None:
        logging.warning("Patient %s not in cohort lookup; using 'cNF_Unknown'",
                        individual_id)
        return "cNF_Unknown"
    return cohort


def build_samples_table(specimens: pd.DataFrame, start_id: int) -> pd.DataFrame:
    """Build the schema-compliant Sample table."""
    rows = []
    next_id = start_id
    for _, row in specimens.iterrows():
        specimen = str(row["specimenID"]).strip()
        patient  = str(row["individualID"]).strip()
        rows.append({
            "improve_sample_id": next_id,
            "other_id":          specimen,
            "other_id_source":   assign_cohort(patient),
            "common_name":       specimen,
            "cancer_type":       CANCER_TYPE,
            "other_names":       patient,
            "species":           SPECIES,
            "model_type":        MODEL_TYPE,
        })
        next_id += 1

    df = pd.DataFrame(rows)
    logging.info("Built %d sample rows; ID range %d–%d",
                 len(df), start_id, next_id - 1)
    return df


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "prev_samples",
        help="Path to previous samples CSV (used to find max improve_sample_id)",
    )
    parser.add_argument(
        "--output", default=DEFAULT_OUTPUT,
        help=f"Output CSV path (default: {DEFAULT_OUTPUT})",
    )
    args = parser.parse_args()

    configure_logging()

    # Synapse auth: requires SYNAPSE_AUTH_TOKEN in env, or a cached config
    syn = synapseclient.Synapse()
    syn.login()  # uses SYNAPSE_AUTH_TOKEN or ~/.synapseConfig

    start_id  = get_starting_id(args.prev_samples)
    specimens = query_drug_screen_specimens(syn)
    samples   = build_samples_table(specimens, start_id)

    samples.to_csv(args.output, index=False)
    logging.info("Wrote %s (%d rows)", args.output, len(samples))


if __name__ == "__main__":
    sys.exit(main() or 0)