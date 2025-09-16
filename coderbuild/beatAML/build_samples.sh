#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running GetBeatAML.py with token and prevSamples $1."
python GetBeatAML.py --token $SYNAPSE_AUTH_TOKEN --samples --prevSamples $1
