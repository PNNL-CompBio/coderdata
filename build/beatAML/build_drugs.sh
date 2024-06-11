python GetBeatAML.py --token $SYNAPSE_AUTH_TOKEN --drugs --drugFile $1
python build_drug_desc.py --drugtable /tmp/beataml_drugs.tsv --desctable /tmp/beataml_drug_desciptors.tsv
