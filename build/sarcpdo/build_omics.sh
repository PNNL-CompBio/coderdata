#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running script with token, curSamples $2, and genes $1."

