#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 04-drug_dosage_and_curves.py with drugfile $2 and curSampleFile $1"
/opt/venv/bin/python 04-drug_dosage_and_curves.py --drugfile $2 --curSampleFile $1