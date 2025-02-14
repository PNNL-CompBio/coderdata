#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running script with token and curSamples $1 and drugFile $2."
python3 03_createSarcPDOExperimentFile.py -t $SYNAPSE_AUTH_TOKEN -s sarcpdo_samples.csv -d sarcpdo_drugs.tsv
