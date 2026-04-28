#!/usr/bin/env python3
"""
build_drugs.py — generate cnf_drugs.tsv and cnf_drug_descriptors.tsv.

Pulls the unique drug names from the cNF drug screen on Synapse, then calls
the standard coderdata utilities:

    coderbuild/utils/pubchem_retrieval.py
    coderbuild/utils/build_drug_descriptor_table.py

These are responsible for:
  * Resolving each drug name to canonical SMILES via PubChem
  * Reusing existing improve_drug_id values when the SMILES already exists
    in any of the prior dataset drug files
  * Computing RDKit descriptors

Usage:
    python build_drugs.py --prev_drugs <file1.tsv,file2.tsv,...>
                          [--out_drugs cnf_drugs.tsv]
                          [--out_desc  cnf_drug_descriptors.tsv]
"""

import argparse
import logging
import os
import subprocess
import sys
import tempfile

import pandas as pd
import synapseclient


# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------
DRUG_SCREEN_TABLE   = "syn51301431"
PUBCHEM_SCRIPT      = os.environ.get(
    "CODERDATA_PUBCHEM_SCRIPT", "/coderbuild/utils/pubchem_retrieval.py")
DESCRIPTOR_SCRIPT   = os.environ.get(
    "CODERDATA_DESCRIPTOR_SCRIPT", "/coderbuild/utils/build_drug_descriptor_table.py")


def configure_logging():
    logging.basicConfig(
        format="[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )


def collect_drug_names(syn) -> list[str]:
    """Pull every drug-screen file from Synapse and collect unique drug names."""
    query = (
        "select id, individualID, specimenID "
        f"from {DRUG_SCREEN_TABLE} "
        "where dataType='drug screen'"
    )
    files = syn.tableQuery(query).asDataFrame()
    logging.info("Reading %d drug-screen files for drug names", len(files))

    names: set[str] = set()
    for _, row in files.iterrows():
        try:
            df = pd.read_csv(syn.get(row["id"]).path)
        except Exception as exc:  # pragma: no cover
            logging.warning("Skipping file %s: %s", row["id"], exc)
            continue
        if "Drug" not in df.columns:
            logging.warning("File %s has no 'Drug' column; skipping", row["id"])
            continue
        names.update(df["Drug"].dropna().astype(str).str.strip().tolist())

    names = sorted(n for n in names if n and n.lower() != "nan")
    logging.info("Collected %d unique drug names", len(names))
    return names


def write_name_file(names: list[str], path: str) -> None:
    """Write a one-name-per-line file for pubchem_retrieval.py."""
    with open(path, "w", encoding="utf-8") as fh:
        for n in names:
            fh.write(n + "\n")


def call_pubchem_retrieval(name_file: str, prev_drug_files: list[str],
                            out_drugs: str) -> None:
    """Call coderbuild/utils/pubchem_retrieval.py to resolve names → SMILES → IDs.

    The exact CLI signature of pubchem_retrieval.py varies; the call below
    matches the documented pattern in the coderdata contribution guide.
    Adjust if the script signature changes upstream.
    """
    cmd = ["python", PUBCHEM_SCRIPT,
           "--input", name_file,
           "--output", out_drugs]
    if prev_drug_files:
        cmd += ["--existing", ",".join(prev_drug_files)]
    logging.info("Running pubchem_retrieval: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def call_descriptor_table(drug_file: str, out_desc: str) -> None:
    """Call build_drug_descriptor_table.py to compute RDKit descriptors."""
    cmd = ["python", DESCRIPTOR_SCRIPT,
           "--input", drug_file,
           "--output", out_desc]
    logging.info("Running descriptor builder: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prev_drugs", default="",
        help="Comma-delimited list of existing drug TSV files (for ID reuse)",
    )
    parser.add_argument("--out_drugs", default="/tmp/cnf_drugs.tsv")
    parser.add_argument("--out_desc",  default="/tmp/cnf_drug_descriptors.tsv")
    args = parser.parse_args()

    configure_logging()

    syn = synapseclient.Synapse()
    syn.login()

    names = collect_drug_names(syn)
    if not names:
        logging.error("No drug names found; aborting")
        sys.exit(1)

    prev = [p.strip() for p in args.prev_drugs.split(",") if p.strip()]

    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as tmp:
        name_file = tmp.name
    try:
        write_name_file(names, name_file)
        call_pubchem_retrieval(name_file, prev, args.out_drugs)
        call_descriptor_table(args.out_drugs, args.out_desc)
    finally:
        os.unlink(name_file)

    logging.info("Wrote %s and %s", args.out_drugs, args.out_desc)


if __name__ == "__main__":
    sys.exit(main() or 0)