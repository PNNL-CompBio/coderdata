---
layout: default
title: CoderData
---

<link rel="stylesheet" href="assets/css/style.css">

<!-- ## Usage of CoderData -->

## Introduction
CoderData is a comprehensive package designed for handling cancer benchmark data in Python. It offers functionalities to download datasets, load them into Python environments, and reformat them according to user needs.

## Downloading Data
The `download` function in CoderData facilitates the downloading of datasets from Figshare. Users can specify a dataset prefix to filter the required files.

### Command Line Usage
To download data via the command line, execute the following command:
<div class="code-box">
    <p>coderdata-download --prefix [PREFIX]</p>
</div>
Replace [PREFIX] with the desired dataset prefix (e.g., 'hcmi', 'beataml'). Omitting the prefix or using 'all' downloads all available datasets.

### Python Usage
In Python, the download process is handled through the download_data_by_prefix function from the downloader module.
<div class="code-box">
    <p>import coderdata as cd</p>
    <p><span class="code-comment"># Download a specific dataset</span></p>
    <p>cd.download_data_by_prefix('beataml')</p>
    <p><span class="code-comment"># Download all datasets</span></p>
    <p>cd.download_data_by_prefix()</p>
</div>

## Loading Data
The DatasetLoader class in CoderData is designed for loading datasets into Python. It automatically initializes attributes for each dataset type like transcriptomics, proteomics, and mutations.
<div class="code-box">
    <p>import coderdata as cd</p>
    <p><span class="code-comment"># Initialize the loader for a specific dataset type</span></p>
    <p>hcmi = cd.DatasetLoader('hcmi')</p>
    <p><span class="code-comment"># Access a datatype of the loaded dataset</span></p>
    <p>transcriptomics_data = hcmi.transcriptomics</p>
    <p><span class="code-comment"># View pandas formatted preview of the samples data</span></p>
    <p>hcmi.samples</p>
    <p><span class="code-comment"># View pandas formatted preview of the transcriptomics data</span></p>
    <p>hcmi.transcriptomics</p>
</div>

## Reformatting Datasets
You can reformat datasets into long or wide formats using the reformat_dataset method. By default, data is in the long format.
Reformatting from long to wide retains three data types, entrez_id and improve_sample_id, value of interest (such as transcriptomics).
<div class="code-box">
    <p><span class="code-comment"># Reformat a specific dataset</span></p>
    <p>hcmi.reformat_dataset('transcriptomics', 'wide') </p>
    <p><span class="code-comment"># Reformat all datasets</span></p>
    <p>hcmi.reformat_dataset('wide')</p>
    <p><span class="code-comment"># Reformat all datatypes back to 'long' datasets</span></p>
    <p>hcmi.reformat_dataset('long') </p>
</div>

## Reloading Datasets
The reload_datasets method is useful for reloading specific datasets or all datasets from local storage, especially if the data files have been updated or altered.
<div class="code-box">
    <p><span class="code-comment"># Reload a specific dataset</span></p>
    <p>hcmi.reload_datasets('transcriptomics')</p>
    <p><span class="code-comment"># Reload all datasets</span></p>
    <p>hcmi.reload_datasets()</p>
</div>

## Information Function (Not yet created)
The information function tells you which datatypes are available and which format they are in for a particular DatasetLoader object.
<div class="code-box">
    <p><span class="code-comment"># Get information about the datasets</span></p>
    <p>hcmi.information()</p>
</div>

## Conclusion
CoderData provides a robust and flexible way to work with cancer benchmark data. By using these functionalities, researchers and data scientists can easily manipulate and analyze complex datasets in their Python environments
