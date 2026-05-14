#!/usr/bin/env bash
# build_drugs.sh — wraps 03-drugs-cnf.py per coderdata convention.
#
# Usage:
#   build_drugs.sh <existing_drug_files>
#
# <existing_drug_files> is a comma-delimited list of drug TSV files from
# previous datasets in the build sequence. Their improve_drug_id values
# get reused when the canonical SMILES match.

set -euo pipefail

PREV_DRUGS="${1:-}"

python 03-drugs-cnf.py \
    --prev_drugs "$PREV_DRUGS" \
    --out_drugs /tmp/cnf_drugs.tsv \
    --out_desc  /tmp/cnf_drug_descriptors.tsv

python build_drug_desc.py --drugtable /tmp/cnf_drugs.tsv --desctable /tmp/cnf_drug_descriptors.tsv.gz
