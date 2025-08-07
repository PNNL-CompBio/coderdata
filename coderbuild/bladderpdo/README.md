# Bladder Patient Derived Organoid Dataset

This dataset comes from a [2018 publication](https://www.sciencedirect.com/science/article/pii/S0092867418302976?via%3Dihub) (Lee et. al., 2018) that created several bladder PDO lines from xenografts in mice and performed sequencing and drug experiments. 

## Samples

Sample information can be found on Synapse [here](https://www.synapse.org/Synapse:syn64765486).

## Omics Data

This study had copy number variation, RNA-seq and mutation data that was included. 

Copy number data is available at [this synapse link](https://www.synapse.org/Synapse:syn64765499). Mutation data is available at [this synapse link](https://www.synapse.org/Synapse:syn64765525). RNA-Seq data was acquired from GEO, we used the 'GSE1039990_Normalized_counts.txt.gz' table available on GEO [GSE103990](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103990). 

## Drug Dose-Response Experiments

For these data, we had complete drug dose-response data to recalculate curves, available on Synapse in the [Lee Bladder PDO Datasets](https://www.synapse.org/Synapse:syn64765430) folder. 

## To Run Code 

### First sample and omics steps are the same, by hand locally or in full coderbuild process
```
python3 00_createBladderPDOSampleFile.py --token $SYNAPSE_AUTH_TOKEN -p prevSamples

# for mutation data (-m)
python3 01_createBladderPDOOmicsFiles.py --token $SYNAPSE_AUTH_TOKEN -s curSamples -g genes.csv -m
# for expressiondata (-e)
python3 01_createBladderPDOOmicsFiles.py --token $SYNAPSE_AUTH_TOKEN -s curSamples -g genes.csv -e
# for copynumber (-c)
python3 01_createBladderPDOOmicsFiles.py --token $SYNAPSE_AUTH_TOKEN -s curSamples -g genes.csv -c
```

### For drug and experiment steps, command depends on location of helper scripts
```
# for running locally (from coderbuild directory):
python3 -m bladderpdo.02_createBladderPDODrugsFile --token $SYNAPSE_AUTH_TOKEN -d prevDrugFilePath -o ./bladderpdo/bladderpdo_drugs.tsv
# for running in Docker as part of full build
python3 02_createBladderPDODrugsFile.py --token $SYNAPSE_AUTH_TOKEN -d prevDrugFilePath -o /tmp/bladderpdo_drugs.tsv

# for running locally (from coderbuild directory): 
python3 utils/build_drug_desc.py --drugtable ./bladderpdo/bladderpdo_drugs.tsv --desctable ./bladderpdo/bladderpdo_drug_descriptors.tsv.gz
# for running in docker as part of full build
python3 build_drug_desc.py --drugtable /tmp/bladderpdo_drugs.tsv --desctable /tmp/bladderpdo_drug_descriptors.tsv.gz


python3 03_createBladderPDOExperimentFile.py --token $SYNAPSE_AUTH_TOKEN --drugfile curDrugFile --curSampleFile curSampleFile --output /tmp/bladderpdo_doserep.tsv

python3 fit_curve.py --input /tmp/bladderpdo_doserep.tsv --output /tmp/bladderpdo_experiments.tsv
```
