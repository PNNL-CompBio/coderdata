---
layout: default
title: CoderData
---

<link rel="stylesheet" href="assets/css/style.css">

## Usage of CoderData

### Introduction
CoderData is a comprehensive package designed for handling cancer benchmark data in Python. It offers functionalities to download datasets, load them into Python environments, and reformat them according to user needs.

### Downloading Data
The `download` function in CoderData facilitates the downloading of datasets from Figshare. Users can specify a dataset prefix to filter the required files.

#### Command Line Usage
To download data via the command line, execute the following command:
<div class="code-box">
    <p>coderdata-download --prefix [PREFIX]</p>
</div>
Replace [PREFIX] with the desired dataset prefix (e.g., 'hcmi', 'beataml'). Omitting the prefix or using 'all' downloads all available datasets.

#### Python Usage
In Python, the download process is handled through the download_data_by_prefix function from the downloader module.
<div class="code-box">
    <p>import coderdata as cd</p>
    <p>cd.download_data_by_prefix('beataml') # Download a specific dataset</p>
    <p>cd.download_data_by_prefix() # Download all datasets</p>
</div>

### Loading Data
The DatasetLoader class in CoderData is designed for loading datasets into Python. It automatically initializes attributes for each dataset type like transcriptomics, proteomics, and mutations.
<div class="code-box">
    <p>import coderdata as cd</p>
    <p>hcmi = cd.DatasetLoader('hcmi') # Initialize the loader for a specific dataset type</p>
    <p>transcriptomics_data = hcmi.transcriptomics # Access the loaded datasets</p>
    <p>hcmi.samples # View pandas formatted preview of the data</p>
    <p>hcmi.transcriptomics</p>
</div>

### Reformatting Datasets
You can reformat datasets into long or wide formats using the reformat_dataset method. By default, data is in the long format.
Reformatting from long to wide retains three data types, entrez_id and improve_sample_id, value of interest (such as transcriptomics).
<div class="code-box">
    <p>hcmi.reformat_dataset('transcriptomics', 'wide') # Reformat a specific dataset</p>
    <p>hcmi.reformat_dataset('wide') # Reformat all datasets</p>
    <p>hcmi.reformat_dataset('long') # Reformat all back to long datasets</p>
</div>

### Reloading Datasets
The reload_datasets method is useful for reloading specific datasets or all datasets from local storage, especially if the data files have been updated or altered.
<div class="code-box">
    <p>hcmi.reload_datasets('transcriptomics') # Reload a specific dataset</p>
    <p>hcmi.reload_datasets() # Reload all datasets</p>
</div>

### Information Function (Not yet created)
The information function tells you which datatypes are available and which format they are in for a particular DatasetLoader object.
<div class="code-box">
    <p>hcmi.information() # Get information about the datasets</p>
</div>

### Conclusion
CoderData provides a robust and flexible way to work with cancer benchmark data. By using these functionalities, researchers and data scientists can easily manipulate and analyze complex datasets in their Python environments