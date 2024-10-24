#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 03a-nci60Drugs.py..."
/opt/venv/bin/python 03a-nci60Drugs.py

echo "Running 03-createDrugFile.R..."
Rscript 03-createDrugFile.R CTRPv2,GDSC,gCSI,PRISM,CCLE,FIMM

echo "Running build_drug_desc.py..."
/opt/venv/bin/python build_drug_desc.py \
  --drugtable /tmp/broad_sanger_drugs.tsv \
  --desctable /tmp/broad_sanger_drug_descriptors.tsv.gz