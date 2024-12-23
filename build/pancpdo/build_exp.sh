#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 04-drug_dosage_and_curves.py with drugfile $2 and curSampleFile $1"
/opt/venv/bin/python 04-getPancPDOExperiments.py --pat $SYNAPSE_AUTH_TOKEN --drugs $2 --samples $1 --output /tmp/pancpdo_doserep.tsv
/opt/venv/bin/python fit_curv.py --input /tmp/panpdo_doserep.tsv --output /tmp/pancpdo_experiments.tsv.gz
