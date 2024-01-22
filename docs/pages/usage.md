---
layout: default
title: CoderData
---

<link rel="stylesheet" href="assets/css/style.css">

<!-- ## Usage of CoderData -->

## Introduction
CoderData is a comprehensive package designed for handling cancer benchmark data in Python.  
It offers functionalities to download datasets, load them into Python environments, and reformat them according to user needs.

## Installation
To install, confirm that you have python avilable and then run the following command in your terminal:

```bash
pip install coderdata
```

## Downloading Data
The `download` function in CoderData facilitates the downloading of datasets from Figshare. Users can specify a dataset prefix to filter the required files.

### Command Line Usage
To download data via the command line, execute the following command:
<div class="code-box">
    <p>coderdata-download --prefix [PREFIX]</p>
</div>
Replace [PREFIX] with the desired dataset prefix (e.g., 'hcmi', 'beataml'). Omitting the prefix or using 'all' downloads all available datasets.

### Python Usage
In Python, the download process is handled through the `download_data_by_prefix` function from the downloader module.
<div class="code-box">
    <p>import coderdata as cd</p>
    <p><span class="code-comment"># Download a specific dataset</span></p>
    <p>cd.download_data_by_prefix('beataml')</p>
    <p><span class="code-comment"># Download all datasets</span></p>
    <p>cd.download_data_by_prefix()</p>
</div>

## Loading Data
The `DatasetLoader` class in CoderData is designed for loading datasets into Python.  
It automatically initializes attributes for each dataset type like transcriptomics, proteomics, and mutations.
<div class="code-box">
    <p>import coderdata as cd</p>
    <p><span class="code-comment"># Initialize the DatasetLoader for a specific dataset type</span></p>
    <p>hcmi = cd.DatasetLoader('hcmi')</p>
    <p><span class="code-comment"># Access a datatype of the loaded dataset</span></p>
    <p>transcriptomics_data = hcmi.transcriptomics</p>
    <p><span class="code-comment"># View pandas formatted preview of the samples data</span></p>
    <p>hcmi.samples</p>
    <p><span class="code-comment"># View pandas formatted preview of the transcriptomics data</span></p>
    <p>hcmi.transcriptomics</p>
</div>

## Joining Datasets
The `join_datasets` function in CoderData is designed for joining and loading datasets in Python with the most flexibility possible.
It is capable of joining initialized, previously joined, or non-initialized datasets. This means you may modify a dataset before joining it with another.
<div class="code-box">
    <p>import coderdata as cd</p>
    <p><span class="code-comment"># Initialize the DatasetLoader for a specific dataset type</span></p>
    <p>hcmi = cd.DatasetLoader('hcmi')</p>
    <p><span class="code-comment"># Access a datatype of the loaded dataset</span></p>
    <p>beataml = cd.DatasetLoader('beataml')</p>
    <p><span class="code-comment"># Join two previously initialized datasets</span></p>
    <p>joined_dataset1 = cd.join_datasets(beataml, hcmi)</p>
    <p><span class="code-comment"># Join a previously joined dataset with a non-initialized dataset</span></p>
    <p><span class="code-comment"># Quotes around a dataset name will load from local files using the DatasetLoader function.</span></p>
    <p>joined_dataset2 = cd.join_datasets(joined_dataset1, "cell_line")</p>
    <p><span class="code-comment"># Join multiple datasets using every method available</span></p>
    <p>joined_dataset3 = cd.join_datasets("cell_line", beataml)</p>
    <p>joined_dataset4 = cd.join_datasets(joined_dataset3, "cptac", hcmi)</p>
</div>

## Reformatting Datasets
You can reformat datasets into long or wide formats using the `reformat_dataset` method. By default, data is in the long format.  
Reformatting from long to wide retains three data types, entrez_id and improve_sample_id, value of interest (such as transcriptomics).  
Datasets cannot be joined while there is a datatype in the wide format.
<div class="code-box">
    <p>import coderdata as cd</p>
    <p><span class="code-comment"># Reformat a specific dataset</span></p>
    <p>hcmi.reformat_dataset('transcriptomics', 'wide') </p>
    <p><span class="code-comment"># Reformat all datasets</span></p>
    <p>hcmi.reformat_dataset('wide')</p>
    <p><span class="code-comment"># Reformat all datatypes back to 'long' datasets</span></p>
    <p>hcmi.reformat_dataset('long') </p>
</div>

## Reloading Datasets
The `reload_datasets` method is useful for reloading specific datasets or all datasets from local storage, especially if the data files have been updated or altered.
<div class="code-box">
    <p>import coderdata as cd</p>
    <p><span class="code-comment"># Reload a specific dataset</span></p>
    <p>hcmi.reload_datasets('transcriptomics')</p>
    <p><span class="code-comment"># Reload all datasets</span></p>
    <p>hcmi.reload_datasets()</p>
</div>

## Info Function 
The `info` method tells you which datatypes are available, their long/wide format, and which datasets they came from.
<div class="code-box">
    <span class="code-comment"># Get information about the joined datasets</span><br>
    joined_dataset4.information()<br>
    <span class="code-comment"># The output is as follows - </span><br>
    <span class="code-comment">
    This is a joined dataset comprising of:<br>
    - beataml: Beat acute myeloid leukemia (BeatAML) data was collected though GitHub and Synapse.<br>
    - hcmi: Human Cancer Models Initiative (HCMI) data was collected though the National Cancer Institute (NCI) Genomic Data Commons (GDC) Data Portal.<br>
    - cell_line: The cell line datasets were collected from numerous resources such as the LINCS project, DepMap, and the Sanger Institute.<br>
    - cptac: The Clinical Proteomic Tumor Analysis Consortium (CPTAC) project is a collaborative network funded by the National Cancer Institute (NCI).<br>

    Available Datatypes and Their Formats<br>
    - copy_number: long format<br>
    - mutations: long format<br>
    - proteomics: long format<br>
    - samples: long format<br>
    - transcriptomics: long format<br>

    Datatype Origins:<br>
    - proteomics: Data from beataml, cell_line, cptac<br>
    - transcriptomics: Data from beataml, cell_line, hcmi, cptac<br>
    - copy_number: Data from cell_line, hcmi, cptac<br>
    - mutations: Data from beataml, cell_line, hcmi, cptac<br>
    - samples: Data from beataml, cell_line, hcmi, cptac 
    </span>
</div>

## Conclusion
CoderData provides a robust and flexible way to work with cancer benchmark data.   
By using these functionalities, researchers and data scientists can easily manipulate and analyze complex datasets in their Python environments
