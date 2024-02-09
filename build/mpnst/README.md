## Build instructions for MPNST dataset

To build the MPNST dataset we have an additional step. 
1. update the syanpse identifier mapping file
   [./synapse_NF-MPNST_sample.csv]()
2. Build the docker image (Run under build directory): 
```
docker build -f dockerfile.mpnst  -t mpnst . --build-arg HTTPS_PROXY=$HTTPS_PROXY 
```
2. Generate new identifers for these samples to create a
   `mpnst_samples.csv` file
```
docker run -v $PWD:/app mpnst Rscript mpnst/sample_gen.R
```
3. Pull the data and map to the samples. 

```
docker run -v $PWD:/app mpnst Rscript mpnst/mpnst_get_rna.R $SYNAPSE_AUTH_TOKEN
docker run -v $PWD:/app mpnst Rscript mpnst/mpnst_get_cnv.R $SYNAPSE_AUTH_TOKEN
docker run -v $PWD:/app mpnst Rscript mpnst/mpnst_get_wes.R $SYNAPSE_AUTH_TOKEN

```
