
/opt/venv/bin/python3 03-getPancPDODrugs.py --pat $SYNAPSE_AUTH_TOKEN --prevDrugFile=$1 --output=/tmp/pancpdo_drugs.tsv.gz
/opt/venv/bin/python3 build_drug_desc.py --drugtable /tmp/pancpdo_drugs.tsv.gz --desctable /tmp/pancpdo_drug_descriptors.tsv.gz
