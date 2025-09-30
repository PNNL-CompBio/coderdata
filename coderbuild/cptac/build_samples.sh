#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running getCptacData.py with prevSampleFile=$1."
python3 getCptacData.py --prevSampleFile=$1
