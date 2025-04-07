#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 00_createBladderPDOSampleFile.py with token and previous sample file $1"
python3 00_createBladderPDOSampleFile.py --token $SYNAPSE_AUTH_TOKEN -p $1