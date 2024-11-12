#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running getCptacData.py with geneFile $1 and curSampleFile=$2."
python getCptacData.py --geneFile $1 --curSampleFile=$2
