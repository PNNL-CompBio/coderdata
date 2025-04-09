#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 03-drug-crc.py with token and prevSamples $1."
python 03-drug-crc.py --Download --Drugs --token $SYNAPSE_AUTH_TOKEN --prevSamples $1

