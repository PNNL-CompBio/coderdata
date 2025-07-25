#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 01-samples-novartispdx.py with token and previous sample file $1"
python3 01-samples-novartispdx.py --token $SYNAPSE_AUTH_TOKEN -p $1