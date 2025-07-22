#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running script with token and drugFile $1"
# for running locally (from build directory):
python3 -m 03-drugs-novartispdx --token $SYNAPSE_AUTH_TOKEN -d $1 -o /tmp/novartispdx_drugs.tsv
#python3 novar
#python3 03-drugs-novartispdx.py --token $SYNAPSE_AUTH_TOKEN -d $1 -o /tmp/novartispdx_drugs.tsv

echo "Running build_drug_desc.py..."
#for running locally: 
python3 utils/build_drug_desc.py --drugtable /tmp/novartispdx_drugs.tsv --desctable /tmp/novartispdx_drug_descriptors.tsv.gz
#python3 build_drug_desc.py --drugtable /tmp/novartispdx_drugs.tsv --desctable /tmp/novartispdx_drug_descriptors.tsv.gz