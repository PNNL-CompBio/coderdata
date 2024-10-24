#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running GetBeatAML.py with token, curSamples $2, and genes $1."
python GetBeatAML.py --token $SYNAPSE_AUTH_TOKEN --omics --curSamples $2 --genes $1
