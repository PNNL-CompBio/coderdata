
# Gene table docker container
This directory contains the data and scripts needed to build a gene
table from Bioconductor


### Build gene table
First we need to build the gene table

1. Build genes docker
```
   docker build -f build/docker/Dockerfile.genes -t genes . --build-arg HTTPS_PROXY=$HTTPS_PROXY 
```

2. Build gene file
```
	docker run -v $PWD:/tmp genes sh build_genes.sh
```
