## Data harmonization workflow

All data collected for this package has been collated from stable/reproducible sources using the scripts contained here. The figure below shows a brief description of the process, which is designed to be run serially, as new identifiers are generated as data are added.

### Directory structure

We have created a separate directory with scripts that collect data from distinct sources as described below.

| Directory | Description | Files generated |
| ---- | ---- | --- |
| [`/cell_line`](./cell_line) | DepMap omics data and multiple drug sensitivity datasets| |
| [`/cptac`](./cptac) | CPTAC omics data measurements| |
| [`/hcmi`](./hcmi)  | HCMI omics data measurements ||
| [`/beatAML`](./beatAML) | BeatAML Leukemia data | | 
| [`/nfosi`](./nfosi) | Data collated through the NF Open Science Initiative ||

### Data ingest process

Currently the data ingest process requires going through each directory and running the scripts there-in. We are hoping to containerize and automate this process in the future but for now running locally is our most reliable processes.



### Data Storage

Final copies of the data will be put onto Figshare. 
