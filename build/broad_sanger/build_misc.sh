#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

cp /tmp/broad_sanger* .

echo "Running 05a_remove_problem_drugs.py..."
/opt/venv/bin/python 05a_remove_problem_drugs.py

echo "Running 05b_separate_datasets.py..."
/opt/venv/bin/python 05b_separate_datasets.py

echo "Removing broad_sanger* files..."
rm broad_sanger*
