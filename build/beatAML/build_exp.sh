#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running GetBeatAML.py with token and curSamples $1 and drugFile $2."
python GetBeatAML.py --exp --token $SYNAPSE_AUTH_TOKEN --curSamples $1 --drugFile $2
