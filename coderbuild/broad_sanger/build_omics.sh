#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

# DepMap files are no longer downloadable from the DepMap portal (Cloudflare
# challenge). Pull the static copy from Synapse first; the R script reads them
# from $DEPMAP_DIR and deletes each one once consumed.
#
# $DEPMAP_DIR is container-local (NOT under /tmp, which is the bind-mounted
# host local/ directory), so these transient inputs never touch host disk.
# The trap is a backstop in case the R script exits before consuming them.
echo "Fetching DepMap files from Synapse..."
DEPMAP_DIR="$(/opt/venv/bin/python fetch_depmap.py --print-dest)"
export DEPMAP_DIR
trap 'rm -rf "$DEPMAP_DIR"' EXIT
/opt/venv/bin/python fetch_depmap.py --dest "$DEPMAP_DIR"

echo "Running 02a-broad_sanger_proteomics.py with gene file $1 and sample file $2."
/opt/venv/bin/python 02a-broad_sanger_proteomics.py --gene $1 --sample $2

echo "Running 02-broadSangerOmics.R with gene file $1 and sample file $2,"
Rscript 02-broadSangerOmics.R $1 $2
