#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 03_get_drug_response_data.R with $SYNAPSE_AUTH_TOKEN, $1, and $2."
Rscript 03_get_drug_response_data.R $SYNAPSE_AUTH_TOKEN $1 $2
