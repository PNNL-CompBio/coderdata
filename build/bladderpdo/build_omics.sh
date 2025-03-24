#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running script with token, curSamples $2, and genes $1."
# for mutation data (-m)
#python3 01_createBladderPDOOmicsFiles.py --token $SYNAPSE_AUTH_TOKEN -s $2 -g $1 -m
# for expressiondata (-e)
#python3 01_createBladderPDOOmicsFiles.py --token $SYNAPSE_AUTH_TOKEN -s $2 -g $1 -e
# for copynumber
python3 01_createBladderPDOOmicsFiles.py --token $SYNAPSE_AUTH_TOKEN -s $2 -g $1 -c
