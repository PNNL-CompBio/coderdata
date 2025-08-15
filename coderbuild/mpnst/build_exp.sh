#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 03_get_experiments.R with $SYNAPSE_AUTH_TOKEN, $1, and $2."
Rscript 03_get_experiments.R $SYNAPSE_AUTH_TOKEN $1 $2 mpnst
rm /tmp/mpnst_pdx_experiments.tsv /tmp/mpnst_mt_experiments.tsv /tmp/mpnst_mt_curve_data.tsv /tmp/mpnst_pdx_curve_data.tsv

