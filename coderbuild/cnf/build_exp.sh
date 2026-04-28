#!/usr/bin/env bash
# build_exp.sh — wraps 04-experiments-cnf.py per coderdata convention.
#
# Usage:
#   build_exp.sh <samples_file> <drugs_file>

set -euo pipefail

SAMPLES="${1:?Usage: build_exp.sh <cnf_samples.csv> <cnf_drugs.tsv>}"
DRUGS="${2:?Usage: build_exp.sh <cnf_samples.csv> <cnf_drugs.tsv>}"

python 04-experiments-cnf.py "$SAMPLES" "$DRUGS" --output /tmp/cnf_experiments.tsv