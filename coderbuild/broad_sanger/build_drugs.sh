#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 03a-nci60Drugs.py..."
/opt/venv/bin/python 03a-nci60Drugs.py --output /tmp/nci60_drugs.tsv

echo "Running 03-createDrugFile.R..."
Rscript 03-createDrugFile.R CTRPv2,GDSC,gCSI,PRISM,CCLE,FIMM /tmp/nci60_drugs.tsv

echo "Merging NCI60 and Sanger/DepMap drug TSVs..."
/opt/venv/bin/python 03_joinDrugFiles.py /tmp/nci60_drugs.tsv /tmp/depmap_sanger_drugs.tsv -o /tmp/broad_sanger_drugs.tsv

echo "Removing temporary NCI60 and depmap drugs file..."
rm /tmp/depmap_sanger_drugs.tsv /tmp/nci60_drugs.tsv

echo "Running build_drug_desc.py..."
/opt/venv/bin/python build_drug_desc.py \
  --drugtable /tmp/broad_sanger_drugs.tsv \
  --desctable /tmp/broad_sanger_drug_descriptors.tsv.gz