## Data harmonization workflow

All data collected for this package has been collated from stable/reproducible sources using the scripts contained here. The figure below shows a brief description of the process, which is designed to be run serially, as new identifiers are generated as data are added.

## build_all.py script

This script initializes all docker containers, builds all datasets, validates them, and uploads them to figshare and pypi.

It requires the following authorization tokens to be set in the local environment depending on the use case:  
`SYNAPSE_AUTH_TOKEN`: Required for beataml and mpnst datasets. Join the [CoderData team](https://www.synapse.org/#!Team:3503472) on Synapse and generate an access token.
`PYPI_TOKEN`: This token is required to upload to PyPI.
`FIGSHARE_TOKEN`: This token is required to upload to Figshare.

Available arguments:

- `--docker`: Initializes and builds all docker containers.
- `--samples`: Processes and builds the sample data files.
- `--omics`: Processes and builds the omics data files.
- `--drugs`: Processes and builds the drug data files.
- `--exp`: Processes and builds the experiment data files.
- `--all`: Executes all available processes above (docker, samples, omics, drugs, exp).
- `--validate`: Validates the generated datasets using the schema check scripts.
- `--figshare`: Uploads the datasets to Figshare.
- `--pypi`: Uploads the package to PyPI.
- `--high_mem`: Utilizes high memory mode for concurrent data processing.
- `--dataset`: Specifies the datasets to process (default='broad_sanger,hcmi,beataml,mpnst,cptac').
- `--version`: Specifies the version number for the package and data upload title. This is required to upload to figshare and PyPI

Example usage:
```bash
python build/build_all.py --all --high_mem --validate --pypi --figshare --version 0.1.29
```

### Directory structure

We have created a separate directory with scripts that collect data from distinct sources as described below.

|Directory | Description | Files generated |
|---- | ---- | --- |
|[`/depmap`](./depmap) | DepMap omics data and multiple drug sensitivity datasets| |
|[`/cptac`](./cptac) | CPTAC omics data measurements| |
|[`/hcmi`](./hcmi)  | HCMI omics data measurements ||
|[`/mpnst`](./mpnst)  | MPNST omics data measurements ||
|[`/beatAML`](./beatAML) | BeatAML Leukemia data | | 
|[`/nfosi`](./nfosi) | Data collated through the NF Open Science Initiative ||

### Data ingest process

Currently the data ingest process requires going through each directory and running the scripts there-in. We are hoping to containerize and automate this process in the future but for now running locally is our most reliable processes.



### Data Storage

Final copies of the data will be put onto Figshare. 
