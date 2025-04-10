#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

# running the drug python script
echo "Running 04-experiments-crc.py with token, samples file $1 and drugs file $2."
python 04-experiments-crc.py --Download --Experiment --token $SYNAPSE_AUTH_TOKEN --Samples $1 --Drugs $2

# running the drug descriptor python script
python fit_curve.py --input /tmp/crc_experiments_for_curve_fitting.tsv --output /tmp/crc_doserep.tsv