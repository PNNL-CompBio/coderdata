---
layout: default
title: CoderData
cell_line_cell_lines: 0
cell_line_genes: 0
cell_line_drugs: 0
cptac_cancer_types: 0
cptac_genes: 0
cptac_drugs: 0
hcmi_cancer_types: 0
hcmi_genes: 0
beataml_cancer_types: 0
beataml_genes: 0
beataml_drugs: 0
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
<!-- 
```python
import coderdata
coderdata.download_data_by_prefix('hcmi') # download all hcmi data to local directory.
hcmi_data = coderdata.DatasetLoader('hcmi') # load hcmi data from local directory into the DatasetLoader object.
hcmi_data.transcriptomics # call transcriptomics data from the DatasetLoader object.
``` -->
<div class="code-box">
    <p>import coderdata</p>
    <p>coderdata.download_data_by_prefix('hcmi') # download all hcmi data to local directory.</p>
    <p>hcmi_data = coderdata.DatasetLoader('hcmi') # load hcmi data from local directory into the DatasetLoader object.</p>
    <p>hcmi_data.transcriptomics # call transcriptomics data from the DatasetLoader object.</p>
</div>


### Datasets

<div class="legend">
    <p>Transcriptomics<span class="dot dot_transcriptomics"></span></p>
    <p>Proteomics<span class="dot dot_proteomics"></span></p>
    <p>Mutations<span class="dot dot_mutations"></span></p>
    <p>Copy Number<span class="dot dot_copy_number"></span></p>
</div>

<div class="dataset-section">
    {% assign datasets = 'cell_line,cptac,hcmi,beataml' | split: ',' %}
    {% for dataset in datasets %}
    <p>Dataset: {{ dataset }}</p>
    <p>Data file: {{ dataset | append: '_table' }}</p>
    {% for row in site.data[dataset | append: '_table'] %}
        <p>{{ row | jsonify }}</p>
    {% endfor %}
    <div class="dataset-container">
        <a href="datasets/{{ dataset }}" class="dataset-link">{{ dataset }}</a>
        <div class="dataset-blurb">
            {% for row in site.data[dataset | append: '_table'] %}
            <!-- {% unless forloop.first %}  -->
                <span class="dot dot_{{ row[0] | downcase }}"></span> 
                <p>{{ row | jsonify }}</p>
            {% endunless %}
            {% endfor %}
            {% case dataset %}
            {% when 'cell_line' %}
                <p>Cell Lines: {{ page.cell_line_cell_lines }} </p>
                <p>Genes: {{ page.cell_line_genes }} </p>
                <p>Drugs: {{ page.cell_line_drugs }} </p>
            {% when 'cptac' %}
                <p>Cancer Types: {{ page.cptac.cancer_types }} </p>
                <p>Genes: {{ page.cptac_genes }} </p>
                <p>Drugs: {{ page.cptac_drugs }} </p>
            {% when 'hcmi' %}
                <p>Cancer Types: {{ page.hcmi_cancer_types }} </p>
                <p>Genes: {{ page.hcmi_genes }} </p>
            {% when 'beataml' %}
                <p>Cancer Types: {{ page.beataml_cancer_types }}</p>
                <p>Genes: {{ page.beataml_genes }}</p>
                <p>Drugs: {{ page.beataml_drugs }}</p>
            {% endcase %}
        </div>
    </div>
    {% endfor %}

</div>


<!-- 
<div class="dataset-section">

    <div class="dataset-container">
        <a href="datasets/cell-line" class="dataset-link">Cell Line</a>
        <div class="dataset-blurb">
            <p>Cell Lines: </p>
            <p>Genes: </p>
            <p>Drugs: </p>
        </div>
    </div>

    <div class="dataset-container">
        <a href="datasets/cptac" class="dataset-link">CPTAC</a>
        <div class="dataset-blurb">
            <p>Cancer Types: </p>
            <p>Genes: </p>
            <p>Drugs: </p>
        </div>
    </div>

    <div class="dataset-container">
        <a href="datasets/hcmi" class="dataset-link">HCMI</a>
        <div class="dataset-blurb">
            <p>Cancer Types: </p>
            <p>Genes: </p>
            <p>Drugs: </p>
        </div>
    </div>

    <div class="dataset-container">
        <a href="datasets/beataml" class="dataset-link">BeatAML</a>
        <div class="dataset-blurb">
            <p>Cancer Types: </p>
            <p>Genes: </p>
        </div>
    </div>

</div> -->

### Data Overview

<div class="flex-container"> 
    <div class="flex-item">
        <embed src="{{ 'assets/stats/Fig0_Overview.pdf' | relative_url }}" type="application/pdf" />
    </div>
    <div class="flex-item">
        <embed src="{{ 'assets/stats/Fig5_Sample_Summary.pdf' | relative_url }}" type="application/pdf" />
    </div>
</div>



### To do - add data types for each daaset like transcriptomics ,proteomics, etc Maybe have a little legend and different colored dots for each in the Datasets section.
### These little dots could just be in the in a line in the box above cancer /cell line types.