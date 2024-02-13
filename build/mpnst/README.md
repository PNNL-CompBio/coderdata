## Build Instructions for MPNST Dataset

To build the MPNST dataset, follow these steps:

1. Build the Docker image:
   ```
   docker build -f dockerfile.mpnst -t mpnst . --build-arg HTTPS_PROXY=$HTTPS_PROXY
   ```

2. Generate new identifiers for these samples to create a
   `MPNST_samples.csv` file. This pulls from the latest synapse
   project metadata table.
   ```
   docker run -v $PWD:/tmp mpnst Rscript 00_sample_gen.R /tmp/beatAML/beataml_samples.csv $SYNAPSE_AUTH_TOKEN
   ```

3. Pull the data and map it to the samples. This uses the metadata
   table pulled above.
   ```
   docker run -v $PWD:/tmp mpnst Rscript 01_mpnst_get_omics.R $SYNAPSE_AUTH_TOKEN /tmp/MPNST_samples.csv /tmp/cell_line/genes.csv
   ```

4. Process drug and experiment data. This uses the metadata from above
   as well as the file directory on synapse:
   ```
   docker run -v $PWD:/tmp mpnst Rscript  02_get_drug_data.R $SYNAPSE_AUTH_TOKEN /tmp/MPNST_samples.csv /tmp/cell_line/drugs.tsv.gz
   ```

Please ensure that each step is followed in order for correct dataset compilation.

## MPNST Dataset Structure
The MPNST dataset includes the following output files:
```
├── MPNST_transcriptomics.csv
├── MPNST_mutations.csv
├── MPNST_copy_number.csv
├── MPNST_drugs.tsv.gz
├── MPNST_experiments.tsv.gz
```

