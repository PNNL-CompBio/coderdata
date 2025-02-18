#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 02-getHCMIData.py for transcriptomics."
python 02-getHCMIData.py -m full_manifest.txt -t transcriptomics -o /tmp/hcmi_transcriptomics.csv.gz -g $1 -s $2

echo "Running 02-getHCMIData.py for copy_number."
python 02-getHCMIData.py -m full_manifest.txt -t copy_number -o /tmp/hcmi_copy_number.csv.gz -g $1 -s $2

echo "Running 02-getHCMIData.py for mutations."
python 02-getHCMIData.py -m full_manifest.txt -t mutations -o /tmp/hcmi_mutations.csv.gz -g $1 -s $2