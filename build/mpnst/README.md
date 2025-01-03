## Build Instructions for MPNST Dataset

To build the MPNST dataset, follow these steps from the coderdata root
directory.

### Step 1: Set the SYNAPSE_AUTH_TOKEN Environment Variable.
This is required to download the data.
```
export SYNAPSE_AUTH_TOKEN="Your Synapse Token"
```
### Step 2: Choose an option below depending on your needs.
---
### Option 1: QuickBuild the test dataset using build_dataset.py

This quick build process does not map sample identifers with previous data versions and is only for personal use.
```
python build/build_dataset.py --dataset mpnst --build 
```
---
### Option 2: Build the test dataset using build_dataset.py with a previous dataset.

This build process assumes you already built or have access to a previously built dataset. This previous dataset must be located in `$PWD/local`. The validate argument ensures the output aligns with the schema.
```
python build/build_dataset.py --dataset mpnst --build --validate --use_prev_dataset beataml
```
---
### Option 3: Build each test file one at a time.
This process does not map sample identifers with previous data versions and is only for personal use.

1. Create an empty local directory in the coderdata root directory.
   ```
   mkdir local
   ```
2. Build the Docker image with the optional HTTPS_PROXY argument:
   ```
   docker build -f build/docker/Dockerfile.mpnst -t mpnst . --build-arg HTTPS_PROXY=$HTTPS_PROXY
   ```

3. Generate new identifiers for these samples to create a
   `mpnst_samples.csv` file. This pulls from the latest synapse
   project metadata table.
   ```

   docker run -v "$PWD/local":/tmp -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN mpnst bash build_samples.sh [Previous Samples file or Empty Quotes ("")]


4. Pull the data and map it to the samples. This uses the metadata
   table pulled above.
   ```
   docker run -v "$PWD/local":/tmp -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN mpnst bash build_omics.sh /tmp/genes.csv /tmp/mpnst_samples.csv 
   ```

5. Process drug data
   ```
   docker run -v "$PWD/local":/tmp -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN  mpnst bash build_drugs.sh [Previous Drugs file or Empty Quotes ("")]
   ```
   
6. Process experiment data. This uses the metadata from above as well as the file directory on synapse:
   ```
   docker run -v "$PWD/local":/tmp -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN mpnst bash build_exp.sh /tmp/mpnst_samples.csv /tmp/mpnst_drugs.tsv
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

