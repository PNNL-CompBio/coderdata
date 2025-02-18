#!/bin/bash
#set -euo pipefail

#trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 02_get_drug_data.R with /tmp/mpnst_drugs.tsv and $1."
Rscript 02_get_drug_data.R /tmp/mpnst_drugs.tsv $1

echo "Running build_drug_desc.py."
/opt/venv/bin/python3 build_drug_desc.py --drugtable /tmp/mpnst_drugs.tsv --desctable /tmp/mpnst_drug_descriptors.tsv.gz
