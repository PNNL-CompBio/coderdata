#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

# running the drug python script
echo "Running 04-experiments-colorectal.py with token, samples file $1 and drugs file $2."
python3 04-experiments-colorectal.py --Download --Experiment --Token $SYNAPSE_AUTH_TOKEN --Samples $1 --Drugs $2

# running the drug descriptor python script
python3 fit_curve.py --input /tmp/colorectal_experiments_for_curve_fitting.tsv --output /tmp/colorectal_experiments.tsv

# change name of script 
mv /tmp/colorectal_experiments.tsv.0 /tmp/colorectal_experiments.tsv
