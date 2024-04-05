# LINCS Perturbation data

The LINCS perturbation data relies on teh depmap_sanger samples and
	drugs but requires additional tools
	
	
## Build Docker

```
docker build -f ../../build/docker/Dockerfile.lincs -t lincs ../../
```

## Get samples

```
docker run -v $PWD:/tmp/ lincs Rscript 01a-pullSamples_LINCS.R /tmp/depmap_sanger_samples.csv

```

## Get Drugs
```
docker run -v $PWD:/tmp lincs /opt/venv/bin/python 01b-pullDrugs_LINCS.py --drugFile /tmp/depmap_sanger_drugs.tsv

```

## Get perturbations

```
docker run -v $PWD:/tmp/ linkcs Rscript 05-LINCS_perturbations.R /tmp/genes.csv /tmp/lincs_drugs.tsv /tmp/lincs_samples.csv

```
