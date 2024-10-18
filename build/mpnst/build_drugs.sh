#!/bin/bash
Rscript 02_get_drug_data.R /tmp/mpnst_drugs.tsv $1
/opt/venv/bin/python3 build_drug_desc.py --drugtable /tmp/mpnst_drugs.tsv --desctable /tmp/mpnst_drug_descriptors.tsv.gz
