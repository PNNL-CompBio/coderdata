# build_omics.sh
#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 01_combined_omics.R with $SYNAPSE_AUTH_TOKEN, $2, and $1."
Rscript 01_combined_omics.R $SYNAPSE_AUTH_TOKEN $2 $1
echo "Running 01b_combined_omics_treated.R with PAT, samples=$2, genes=$1"
Rscript 01b_combined_omics_treated.R "$SYNAPSE_AUTH_TOKEN" "$2" "$1"