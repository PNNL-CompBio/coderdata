#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 00_sample_gen.R with $1."
Rscript 00_sample_gen.R $1
