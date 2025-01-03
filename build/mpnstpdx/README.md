## Build Instructions for MPNST PDX Dataset

To build the MPNST PDX dataset, follow these steps from the coderdata root
directory. Currently using the test files as input. 

1. Build the Docker image:
   ```
   docker build -f build/docker/Dockerfile.mpnstpdx -t mpnstpdx . --build-arg HTTPS_PROXY=$HTTPS_PROXY
   ```

2. Generate new identifiers for these samples to create a
   `mpnstpdx_samples.csv` file. This pulls from the latest synapse
   project metadata table.
   ```
   docker run -v $PWD:/tmp -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN mpnstpdx sh build_samples.sh /tmp/build/build_test/test_samples.csv
   ```

3. Pull the data and map it to the samples. This uses the metadata
   table pulled above.
   ```
   docker run -v $PWD:/tmp -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN mpnstpdx sh build_omics.sh /tmp/build/build_test/test_genes.csv /tmp/mpnstpdx_samples.csv 
   ```

4. Process drug data
   ```
   docker run -v $PWD:/tmp -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN  mpnstpdx sh build_drugs.sh /tmp/build/build_test/test_drugs.tsv
   ```
   
5. Process experiment data. This uses the metadata from above as well as the file directory on synapse:
   ```
   docker run -v $PWD:/tmp -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN mpnstpdx sh build_exp.sh /tmp/mpnstpdx_samples.csv /tmp/mpnstpdx_drugs.tsv.gz
   ```

Please ensure that each step is followed in order for correct dataset compilation.

## MPNST PDX Dataset Structure
The MPNST dataset includes the following output files:
```
├── mpnstpdx_samples.csv
├── mpnstpdx_transcriptomics.csv
├── mpnstpdx_mutations.csv
├── mpnstpdx_copy_number.csv
├── mpnstpdx_drugs.tsv
├── mpnstpdx_drug_descriptors.tsv.gz
├── mpnstpdx_experiments.tsv.gz
```

