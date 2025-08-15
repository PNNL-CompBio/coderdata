
python 03-getPancreaticDrugs.py --prevDrugFile=$1 --output=/tmp/pancreatic_drugs.tsv
python build_drug_desc.py --drugtable /tmp/pancreatic_drugs.tsv --desctable /tmp/pancreatic_drug_descriptors.tsv.gz

