## BeatAML Data generation

This directory builds the data for the BeatAML samples. To build and
test this module, run the following commands from the root directory.

### Build with test data
Build commands should be similar to every other coderdata build module.

1. Build the Docker image:
   ```
   docker build -f build/docker/Dockerfile.beataml -t beataml . --build-arg HTTPS_PROXY=$HTTPS_PROXY 
   ```

2. Generate new identifiers for these samples to create a
   `beataml_samples.csv` file. This pulls from the latest synapse
   project metadata table.
   ```
   docker run -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN -v $PWD:/tmp beataml sh build_samples.sh /tmp/build/build_test/test_samples.csv 
   ```

3. Pull the data and map it to the samples. This uses the metadata
   table pulled above.
   ```
   docker run -v $PWD:/tmp  -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN beataml sh build_omics.sh /tmp/build/build_test/test_genes.csv /tmp/beataml_samples.csv 
   ```

4. Process drug data
   ```
   docker run  -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN -v $PWD:/tmp beataml sh build_drugs.sh /tmp/build/build_test/test_drugs.tsv
   ```
   
5. Process experiment data. This uses the metadata from above as well as the file metadata on synapse:
   ```
   docker run  -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN -v $PWD:/tmp beataml sh build_exp.sh /tmp/beataml_samples.csv /tmp/beataml_drugs.tsv.gz
   ```

Please ensure that each step is followed in order for correct dataset
compilation.


### BeatAML Dataset structure
The build commands above create the following files

```
├── beataml_samples.csv
├── beataml_transcriptomics.csv
├── beataml_mutations.csv
├── beataml_drugs.tsv.gz
├── beataml_drug_descriptors.tsv.gz
├── beataml_experiments.tsv.gz
```
