#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running script with token, curSamples $2, and genes $1."
# for mutation data (-m)
python3 02-omics-novartis.py --download --transcriptomics --mutations --copy_number --samples $2 --genes $1 --token $SYNAPSE_AUTH_TOKEN