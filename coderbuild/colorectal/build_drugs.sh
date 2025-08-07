#!/bin/bash
set -euo pipefail
echo "the variable is $1"

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

# running the drug python script
echo "Running 03-drug-colorectal.py with token and PrevDrugs $1."
python3 03-drug-colorectal.py --Download --Drug --Token $SYNAPSE_AUTH_TOKEN --PrevDrugs $1

# running the drug descriptor python script
python3 build_drug_desc.py --drugtable /tmp/colorectal_drugs.tsv --desctable /tmp/colorectal_drug_descriptors.tsv.gz