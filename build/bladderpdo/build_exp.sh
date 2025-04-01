#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 04-drug_dosage_and_curves.py with drugfile $2 and curSampleFile $1"
python 03_createBladderPDOExperimentFile.py --token $SYNAPSE_AUTH_TOKEN --drugfile $2 --curSampleFile $1 --output /tmp/bladderpdo_doserep.tsv

python fit_curve.py --input /tmp/bladderpdo_doserep.tsv --output /tmp/bladderpdo_doserep.tsv