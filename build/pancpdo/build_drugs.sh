
python 03-getPancPDODrugs.py --prevDrugFile=$1 --output=/tmp/pancpdo_drugs.tsv
python build_drug_desc.py --drugtable /tmp/pancpdo_drugs.tsv --desctable /tmp/pancpdo_drug_descriptors.tsv.gz

