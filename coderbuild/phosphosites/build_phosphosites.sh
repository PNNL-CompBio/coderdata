#!/usr/bin/env bash
# build_phosphosites.sh: build the global phosphosites reference CSV.
#
# Usage:
#   build_phosphosites.sh <genes.csv> [<prev_phosphosites.csv>]
#
# Sources (in priority order):
#   1. Ochoa et al. (2020) Nat Biotechnol Supplementary Table S3 (~116k sites)
#   2. UniProt PTM annotations (~12k additional sites unique to UniProt)
#   3. Synapse supplement: syn70078415 (cnf raw phospho — captures experiment-
#      specific sites not yet in either database)

set -euo pipefail

GENES="${1:?Usage: build_phosphosites.sh <genes.csv> [<prev_phosphosites.csv>]}"
PREV="${2:-}"

PREV_ARG=""
if [ -n "$PREV" ] && [ -f "$PREV" ]; then
    PREV_ARG="--prev $PREV"
fi

python /coderbuild/phosphosites/00-buildPhosphositeFile.py \
    "$GENES" \
    --out /tmp/phosphosites.csv \
    --synapse_supplements syn70078415 \
    --synapse_site_col site \
    $PREV_ARG

echo "Phosphosites built from Ochoa + UniProt PTM + Synapse supplement."
