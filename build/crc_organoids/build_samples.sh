#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 01-createSamples-crc.py with token and prevSamples $1."
#   download the data and then create sample sheet  
python3 01-samples-crc.py  --download --samples --token $SYNAPSE_AUTH_TOKEN --prevSamples $1