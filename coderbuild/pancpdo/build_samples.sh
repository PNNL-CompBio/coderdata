#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 01-createPancPDOSamplesFile.py with prevSamples $1."
python 01-createPancPDOSamplesFile.py --prevSamples $1
