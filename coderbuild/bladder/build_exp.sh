#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 04-drug_dosage_and_curves.py with drugfile $2 and curSampleFile $1"
python3 03_createBladderExperimentFile.py --token $SYNAPSE_AUTH_TOKEN --drugfile $2 --curSampleFile $1 --output /tmp/bladder_doserep.tsv

python3 fit_curve.py --input /tmp/bladder_doserep.tsv --output /tmp/bladder_experiments.tsv
rm /tmp/bladder_doserep.tsv
mv /tmp/bladder_experiments.tsv.0 /tmp/bladder_experiments.tsv