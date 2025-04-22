#!/bin/bash
set -euo pipefail
echo "the variable is $1"

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

# running the drug python script
echo "Running 03-drug-crc.py with token and PrevDrugs $1."
python3 03-drug-crc.py --Download --Drug --Token $SYNAPSE_AUTH_TOKEN --PrevDrugs $1

# running the drug descriptor python script
python3 build_drug_desc.py --drugtable /tmp/crc_organoids_drugs.tsv --desctable /tmp/crc_organoids_drug_descriptors.csv.gz