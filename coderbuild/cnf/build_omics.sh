#!/usr/bin/env bash
# build_omics.sh - wraps 02-omics-cnf.py per coderdata convention.
#
# Usage:
#   build_omics.sh <genes_file> <samples_file> <phosphosites_file>
#
# - <genes_file>:       genes.csv with gene_symbol → entrez_id mappings
# - <samples_file>:     cnf_samples.csv produced by build_samples.sh
# - <phosphosites_file>: phosphosites.csv (other_id → phosphosite_id)

set -euo pipefail

GENES="${1:?Usage: build_omics.sh <genes.csv> <cnf_samples.csv> <phosphosites.csv>}"
SAMPLES="${2:?Usage: build_omics.sh <genes.csv> <cnf_samples.csv> <phosphosites.csv>}"
PHOSPHOSITES="${3:?Usage: build_omics.sh <genes.csv> <cnf_samples.csv> <phosphosites.csv>}"

python 02-omics-cnf.py "$GENES" "$SAMPLES" "$PHOSPHOSITES" \
    --out_rna     /tmp/cnf_transcriptomics.csv \
    --out_prot    /tmp/cnf_proteomics.csv \
    --out_phospho /tmp/cnf_phosphoproteomics.csv