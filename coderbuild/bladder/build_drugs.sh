#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running script with token and drugFile $1"
# for running locally (from build directory):
#python3 -m bladderpdo.02_createBladderPDODrugsFile --token $SYNAPSE_AUTH_TOKEN -d $1 -o ./bladderpdo/bladderpdo_drugs.tsv
python3 02_createBladderDrugsFile.py --token $SYNAPSE_AUTH_TOKEN -d $1 -o /tmp/bladder_drugs.tsv

echo "Running build_drug_desc.py..."
#for running locally: 
#python3 utils/build_drug_desc.py --drugtable ./bladderpdo/bladderpdo_drugs.tsv --desctable ./bladderpdo/bladderpdo_drug_descriptors.tsv.gz
python3 build_drug_desc.py --drugtable /tmp/bladder_drugs.tsv --desctable /tmp/bladder_drug_descriptors.tsv.gz