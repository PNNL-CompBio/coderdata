#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 01-createHCMISamplesFile.py with prevSamples $1."
python 01-createHCMISamplesFile.py --prevSamples $1
