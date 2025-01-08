#!/bin/bash
set -euo pipefail

trap 'echo "Error on or near line $LINENO while executing: $BASH_COMMAND"; exit 1' ERR

echo "Running 02-getPancPDOData.py for transcriptomics."
python 02-getPancPDOData.py -m full_manifest.txt -t transcriptomics -o /tmp/pancpdo_transcriptomics.csv.gz -g $1 -s $2

#echo "Running 02-getPancPDOData.py for copy_number."
#python 02-getPancPDOData.py -m full_manifest.txt -t copy_number -o /tmp/pancpdo_copy_number.csv.gz -g $1 -s $2

#echo "Running 02-getPancPDOData.py for mutations."
#python 02-getPancPDOData.py -m full_manifest.txt -t mutations -o /tmp/pancpdo_mutations.csv.gz -g $1 -s $2
