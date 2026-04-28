#!/usr/bin/env bash
# build_omics.sh — wraps 02-omics-cnf.py per coderdata convention.
#
# Usage:
#   build_omics.sh <genes_file> <samples_file>
#
# - <genes_file>: genes.csv with gene_symbol → entrez_id mappings (from prior coderdata builds)
# - <samples_file>: cnf_samples.csv produced by build_samples.sh

set -euo pipefail

GENES="${1:?Usage: build_omics.sh <genes.csv> <cnf_samples.csv>}"
SAMPLES="${2:?Usage: build_omics.sh <genes.csv> <cnf_samples.csv>}"

python 02-omics-cnf.py "$GENES" "$SAMPLES" \
    --out_rna  /tmp/cnf_transcriptomics.csv \
    --out_prot /tmp/cnf_proteomics.csv