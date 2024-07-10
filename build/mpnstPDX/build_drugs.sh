Rscript 02_get_drug_data.R $SYNAPSE_AUTH_TOKEN $1 /tmp/mpnstpdx_drugs.tsv
/opt/venv/bin/python3 build_drug_desc.py --drugtable /tmp/mpnstpdx_drugs.tsv --desctable /tmp/mpnstpdx_drug_descriptors.tsv.gz
