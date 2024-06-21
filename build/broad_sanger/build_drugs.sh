/opt/venv/bin/python 03a-nci60Drugs.py 
Rscript 03-createDrugFile.R CTRPv2,GDSC,gCSI,PRISM,CCLE,FIMM
yes
/opt/venv/bin/python build_drug_desc.py --drugtable /tmp/broad_sanger_drugs.tsv --desctable /tmp/broad_sanger_drug_descriptors.tsv.gz
