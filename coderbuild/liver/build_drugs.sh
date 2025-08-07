#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

# running the drug python script
echo "Running 03-drug-liver.py with token and PrevDrugs $1."
python3 03-drug-liver.py --Download --Drug --Token $SYNAPSE_AUTH_TOKEN --PrevDrugs $1

# running the drug descriptor python script
python3 build_drug_desc.py --drugtable /tmp/liver_drugs.tsv --desctable /tmp/liver_drug_descriptors.tsv.gz