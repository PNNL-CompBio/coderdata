#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running GetBeatAML.py with token and drugFile $1"
python GetBeatAML.py --token $SYNAPSE_AUTH_TOKEN --drugs --drugFile $1

echo "Running build_drug_desc.py..."
python build_drug_desc.py --drugtable /tmp/beataml_drugs.tsv --desctable /tmp/beataml_drug_descriptors.tsv.gz
