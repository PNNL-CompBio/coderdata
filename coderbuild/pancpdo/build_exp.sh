#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 04-drug_dosage_and_curves.py with drugfile $2 and curSampleFile $1"
python 04-getPancPDOExperiments.py --pat $SYNAPSE_AUTH_TOKEN --drugs $2 --samples $1 --output /tmp/pancpdo_doserep.tsv
python fit_curve.py --input /tmp/pancpdo_doserep.tsv --output /tmp/pancpdo_doserep.tsv

##now move file and gzip
mv /tmp/pancpdo_doserep.tsv.0 /tmp/pancpdo_experiments.tsv

python 05-addPrecalcAUC.py --samples $1 --drugs $2 --expfile /tmp/pancpdo_experiments.tsv
gzip /tmp/pancpdo_experiments.tsv
