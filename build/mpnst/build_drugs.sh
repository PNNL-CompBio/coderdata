Rscript 02_get_drug_data.R $SYNAPSE_AUTH_TOKEN $1 /tmp/mpnst_drugs.tsv
/opt/venv/bin/python build_drug_desc.py --drugtable /tmp/mpnst_drugs.tsv --desctable /tmp/mpnst_drug_descriptors.tsv.gz
