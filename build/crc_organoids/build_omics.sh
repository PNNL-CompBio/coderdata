#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 02-omics-crcPDO.py with token, curSamples $2, and genes $1."
python3 02-omics-crcPDO.py --parse --transcriptomics --mutations --copy_number --ids $2 --genes $1