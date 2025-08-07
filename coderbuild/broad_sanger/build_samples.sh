#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 01-broadSangerSamples.R."
Rscript 01-broadSangerSamples.R
