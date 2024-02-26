## Building DepMap data
The DepMap data is the first to be built, and requires the
following commands. All scripts write files in to the `/tmp/`
directory, so mounting to that directly will help output the files


### Docker image and gene file
First step is to build the docker file and the genes.csv file. This is
required for all future data files.
```
docker build -f ../../build/docker/Dockerfile.depmap -t depmap ../../
docker run -v $PWD:/tmp/ depmap Rscript 00-buildGeneFile.R

```

### DepMap reference samples and identifiers
Next we retrieve all the standard cell line identifiers we can, from diverse
sources, and map them to IMPROVE sample identifiers for future reference.
```
docker run -v $PWD:/tmp/ depmap Rscript 01-depmapSamples.R

```

### Omics data for DepMap cell lines
Third we collect the omics data for these cell lines, again from
diverse sources. Currently we have a single script for each
source. Each script takes our list of gene and sample identifiers
```
docker run -v $PWD:/tmp/ depmap Rscript 02-pullDepMap.R /tmp/genes.csv /tmp/depmap_samples.csv
docker run -v $PWD:/tmp/ depmap Rscript 02b-pullSanger.R /tmp/genes.csv /tmp/depmap_samples.csv

```

### Drug data for all experiments

Fourth we collect drugs and map them to pubchem, then structure. This
is a slow step as we collect from diverse studies including:
1. CTRPv2
2. GDSCv1
3. GDSCv2
4. gCSI
5. PRISM2020
6. CCLE
7. FIMM
8. NCI60

```
docker run -v $PWD:/tmp depmap Rscript 03-createDrugFile.R CTRPv2,GDSC,gCSI,PRISM,CCLE,FIMM,NCI60

```
### Dose response and curve fitting
This last command will generate the drug file, drugs.tsv.gz, which we
can pass into the next commands. Then we will collect the dose
response data and fit the curves for the following experiments:
1. CTRPv2
2. GDSCv1
3. GDSCv2
4. gCSI
5. PRISM2020
6. CCLE
7. FIMM
8. NCI60

```
docker run -v $PWD:/tmp/ depmap /opt/venv/bin/python 04-drug_dosasge_and_curves.py --drugfile=/tmp/drugs.tsv.gz --curSampleFile=/tmp/cell_line_samples.csv

```

### Cell line perturbation data
```
docker run -v $PWD:/tmp depmap /opt/venv/bin/python 04-cellLineDrugs_LINCS.py --drugFile /tmp/drugs.tsv.gz
docker run -v $PWD:/tmp/ depmap Rscript 05-LINCS_perturbations.R /tmp/genes.csv /tmp/drugs.tsv.gz /tmp/cell_line_samples.csv

```


