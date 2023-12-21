---
layout: default
title: CoderData
---

<link rel="stylesheet" href="assets/css/style.css">

<!-- # Cancer Omics and Drug Experiment Response Data (`coderdata`) Python Package -->

### Introduction
CoderData is a cancer benchmark data package developed in Python and R. 
There are two aspects of this package, the backend build section and the user facing python package.
The build section is a github workflow that generates four cancer datasets in a format that is easy for users and algorithms to ingest. 
The python package allows users to easily download the data, load it into python and reformat it as desired.

### Installation
To install `coderdata`, simply run the following command in your terminal:

```bash
pip install coderdata
```

### Usage
##### Bash / Command line
To download datasets, simply run the following command in your terminal. Remove the prefix argument if you'd like to install all datasets.

```bash
coderdata download --prefix hcmi
```

##### Python
To download, load, and call datasets in python, simply run the following commands. 

```python
import coderdata
coderdata.download_data_by_prefix('hcmi') # download all hcmi data to local directory.
hcmi_data = coderdata.DatasetLoader('hcmi') # load hcmi data from local directory into the DatasetLoader object.
hcmi_data.transcriptomics # call transcriptomics data from the DatasetLoader object.
```


### Datasets

<div class="dataset-section">

    <div class="dataset-container">
        <a href="#datasets/cell-line" class="dataset-link">Cell Line</a>
        <div class="dataset-blurb">
            <p>Cell Lines: </p>
            <p>Genes: </p>
            <p>Drugs: </p>
        </div>
    </div>

    <div class="dataset-container">
        <a href="#datasets/cptac" class="dataset-link">CPTAC</a>
        <div class="dataset-blurb">
            <p>Cancer Types: </p>
            <p>Genes: </p>
            <p>Drugs: </p>
        </div>
    </div>

    <div class="dataset-container">
        <a href="#datasets/hcmi" class="dataset-link">HCMI</a>
        <div class="dataset-blurb">
            <p>Cancer Types: </p>
            <p>Genes: </p>
            <p>Drugs: </p>
        </div>
    </div>

    <div class="dataset-container">
        <a href="#datasets/beataml" class="dataset-link">BeatAML</a>
        <div class="dataset-blurb">
            <p>Cancer Types: </p>
            <p>Genes: </p>
        </div>
    </div>

</div>
