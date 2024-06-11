## Building Broad and Sanger cell line data
The Broad and Sanger data is the first to be built, and requires the
following commands. All scripts write files in to the `/tmp/`
directory, so mounting to that directly will help output the files. We
broke the Docker image into two to reduce overall size and complexity
of each image. 


### Build gene, sample, and omics data
Below are the steps required to build and test the gene/sample/omics
builds. Commands are designed to be run from the root of the repo. 

1. Build omics docker
```
   docker build -f build/docker/Dockerfile.broad_sanger_omics -t broad_sanger_omics . --build-arg HTTPS_PROXY=$HTTPS_PROXY 
```
2. Build gene file
```
	docker run -v $PWD:/tmp broad_sanger_omics sh build_genes.sh
```

3. Build sample file
```
  docker run -v $PWD:/tmp broad_sanger_omics sh build_samples.sh
```
4. Build omics files
```
  docker run -v $PWD:/tmp broad_sanger_omics sh build_omics.sh
```
This should leave you with the following files for al the cell lines
```
├── broad_sanger_samples.csv.gz
├── broad_sanger_transcriptomics.csv.gz
├── broad_sanger_mutations.csv.gz
├── broad_sanger_copy_number.csv.gz
├── genes.csv

```

### Build out drug files and experiments
Both of these steps can be lengthy - the experiment fitting can be
parallelized but the drug information requires querying PubChem which
can take a while.

1. Build experiment docker fille

```
   docker build -f build/docker/Dockerfile.broad_sanger_exp -t broad_sanger_exp . --build-arg HTTPS_PROXY=$HTTPS_PROXY 


```
2. Build drug files
   ```
   docker run -v $PWD:/tmp broad_sanger_exp sh build_drugs.sh /tmp/build/build_test/test_drugs.tsv
   ```
3. Build experiment files
   ```
   docker run -v $PWD:/tmp   broad_sanger_exp sh build_exp.sh /tmp/broad_sanger_samples.csv /tmp/broad_sanger_drugs.tsv.gz
   ```

### Datasets collected

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


