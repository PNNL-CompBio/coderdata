#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 04-drug_dosage_and_curves.py with drugfile $2 and curSampleFile $1"
python 04-getPancreaticExperiments.py --pat $SYNAPSE_AUTH_TOKEN --drugs $2 --samples $1 --output /tmp/pancreatic_doserep.tsv
python fit_curve.py --input /tmp/pancreatic_doserep.tsv --output /tmp/pancreatic_doserep.tsv

##now move file and gzip
mv /tmp/pancreatic_doserep.tsv.0 /tmp/pancreatic_experiments.tsv

python 05-addPrecalcAUC.py --samples $1 --drugs $2 --expfile /tmp/pancreatic_experiments.tsv
gzip /tmp/pancreatic_experiments.tsv
