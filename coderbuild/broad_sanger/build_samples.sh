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

echo "Running 01-broadSangerSamples.R."
Rscript 01-broadSangerSamples.R
