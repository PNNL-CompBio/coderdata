## Build Instructions for MPNST Dataset

To build the MPNST dataset, follow these steps, including an additional step specific to this dataset:

1. Update the Synapse identifier mapping file:
   - [./synapse_NF-MPNST_sample.csv]()
   - [./no_combo_manifest_updated_2_5_24.csv]() (this file only contains single treatment files' Synapse IDs).

2. Build the Docker image:
   ```
   docker build -f dockerfile.mpnst -t mpnst . --build-arg HTTPS_PROXY=$HTTPS_PROXY
   ```

3. Generate new identifiers for these samples to create a `mpnst_samples.csv` file:
   ```
   docker run -v $PWD:/app mpnst Rscript mpnst/00_sample_gen.R
   ```

4. Pull the data and map it to the samples:
   ```
   docker run -v $PWD:/app mpnst Rscript mpnst/01_mpnst_get_rna.R $SYNAPSE_AUTH_TOKEN
   docker run -v $PWD:/app mpnst Rscript mpnst/02_mpnst_get_cnv.R $SYNAPSE_AUTH_TOKEN
   docker run -v $PWD:/app mpnst Rscript mpnst/03_mpnst_get_wes.R $SYNAPSE_AUTH_TOKEN
   ```

5. Process drug and experiment data:
   ```
   docker run -v $PWD:/app mpnst Rscript mpnst/04_mpnst_get_synapse.R $SYNAPSE_AUTH_TOKEN
   docker run -v $PWD:/app mpnst /opt/venv/bin/python utils/fit_curve.py --input tmp_drug/combined_data.tsv --output tmp_drug/zzz.out
   docker run -v $PWD:/app mpnst Rscript mpnst/05_mpnst_get_drug_experiment.R
   ```

Please ensure that each step is followed in order for correct dataset compilation.

## MPNST Dataset Structure
The MPNST dataset includes the following output files:
```
├── MPNST_RNA_seq.csv
├── MPNST_WES_mutation_seq.csv
├── MPNST_cnn_mutation_seq.csv
├── MPNST_drugs.tsv.gz
├── MPNST_experiments.csv.gz
```
