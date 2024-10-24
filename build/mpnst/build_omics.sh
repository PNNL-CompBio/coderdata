#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 01_mpnst_get_omics.R with $SYNAPSE_AUTH_TOKEN, $2, and $1."
Rscript 01_mpnst_get_omics.R $SYNAPSE_AUTH_TOKEN $2 $1
