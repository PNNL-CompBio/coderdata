## Build Instructions for MPNST Dataset

To build the MPNST dataset, follow these steps from the coderdata root
directory. Currently using the test files as input. 

1. Build the Docker image:
   ```
   docker build -f build/docker/Dockerfile.mpnst -t mpnst . --build-arg HTTPS_PROXY=$HTTPS_PROXY
   ```

2. Generate new identifiers for these samples to create a
   `mpnst_samples.csv` file. This pulls from the latest synapse
   project metadata table.
   ```
   docker run -v $PWD:/tmp -e -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN mpnst sh build_samples.sh /tmp/build/build_test/test_samples.csv
   ```

3. Pull the data and map it to the samples. This uses the metadata
   table pulled above.
   ```
   docker run -v $PWD:/tmp -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN mpnst sh build_omics.sh /tmp/build/build_test/test_genes.csv /tmp/mpnst_samples.csv 
   ```

4. Process drug data
   ```
   docker run -v $PWD:/tmp -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN  mpnst sh build_drugs.sh /tmp/build/build_test/test_drugs.tsv
   ```
   
5. Process experiment data. This uses the metadata from above as well as the file directory on synapse:
   ```
   docker run -v $PWD:/tmp -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN mpnst sh build_exp.sh /tmp/mpnst_samples.csv /tmp/mpnst_drugs.tsv.gz
   ```

Please ensure that each step is followed in order for correct dataset compilation.

## MPNST Dataset Structure
The MPNST dataset includes the following output files:
```
├── mpnst_samples.csv
├── mpnst_transcriptomics.csv
├── mpnst_mutations.csv
├── mpnst_copy_number.csv
├── mpnst_drugs.tsv
├── mpnst_drug_descriptors.tsv.gz
├── mpnst_experiments.tsv.gz
```

