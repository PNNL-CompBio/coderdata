#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 02a-broad_sanger_proteomics.py with gene file $1 and sample file $2."
/opt/venv/bin/python 02a-broad_sanger_proteomics.py --gene $1 --sample $2

echo "Running 02-broadSangerOmics.R with gene file $1 and sample file $2,"
Rscript 02-broadSangerOmics.R $1 $2
