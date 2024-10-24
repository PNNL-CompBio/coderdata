#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 05_separate_datasets.py..."
/opt/venv/bin/python 05_separate_datasets.py

echo "Removing broad_sanger* files..."
rm broad_sanger*
