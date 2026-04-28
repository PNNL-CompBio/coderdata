#!/usr/bin/env bash
# build_samples.sh — wraps 01-samples-cnf.py per coderdata convention.
#
# Usage:
#   build_samples.sh <previous_samples_file>
#
# The previous samples file is the latest samples CSV from the coderdata
# build (with the highest improve_sample_id used so far). 01-samples-cnf.py
# reads it, finds the max ID, and increments from there.

set -euo pipefail

PREV_SAMPLES="${1:?Usage: build_samples.sh <previous_samples.csv>}"

python 01-samples-cnf.py "$PREV_SAMPLES" --output /tmp/cnf_samples.csv