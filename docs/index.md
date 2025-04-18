---
layout: default
title: CoderData
---

<link rel="stylesheet" href="assets/css/style.css">

<!-- # Cancer Omics and Drug Experiment Response Data (`coderdata`) Python Package -->

## Introduction
CoderData is a cancer benchmark data package developed in Python and R. 
There are two aspects of this package, the backend build section and the user facing python package.
The build section is a github workflow that generates four cancer datasets in a format that is easy for users and algorithms to ingest. 
The python package allows users to easily download the data, load it into python and reformat it as desired.

## Installation and Usage
### Installation

Assuming `python>=3.9` is installed on the system, simply run the following command in the terminal to install the most recent release of the coderdata API:

```bash
$ pip install coderdata
```

### Bash / Command line
A full list of available datasets can be retrieved via:
```sh
$ coderdata --list
```

To download datasets, simply run the following command in your terminal substituting `<DATASET>` with the desired dataset (e.g. `beataml`). To download all datasets use `--name all`.

```bash
$ coderdata download --name <DATASET>
```

### Python

To download, load, and call datasets in python, simply run the following commands.

```python
>>> import coderdata as cd
>>> cd.download(name='beataml')
>>> beataml = cd.load('beataml')
>>> beataml.experiments
         source  improve_sample_id improve_drug_id    study  time time_unit dose_response_metric  dose_response_value
0       synapse               3907       SMI_11123  BeatAML    72       hrs              fit_auc               0.0564
1       synapse               3907       SMI_11211  BeatAML    72       hrs              fit_auc               0.9621
2       synapse               3907       SMI_12192  BeatAML    72       hrs              fit_auc               0.1691
3       synapse               3907       SMI_12254  BeatAML    72       hrs              fit_auc               0.4245
4       synapse               3907       SMI_12469  BeatAML    72       hrs              fit_auc               0.7397
...         ...                ...             ...      ...   ...       ...                  ...                  ...
233775  synapse               3626        SMI_7110  BeatAML    72       hrs                  dss               0.0000
233776  synapse               3626        SMI_7590  BeatAML    72       hrs                  dss               0.0000
233777  synapse               3626        SMI_8159  BeatAML    72       hrs                  dss               0.1946
233778  synapse               3626        SMI_8724  BeatAML    72       hrs                  dss               0.0000
233779  synapse               3626         SMI_987  BeatAML    72       hrs                  dss               0.7165

[233780 rows x 8 columns]
```

For more indepth instructions view our [Usage](pages/usage.md) page.


## Datasets

<table>
  <thead>
    <tr>
      <th>Dataset</th>
      <th>Cancer Types</th>
      <th>Samples</th>
      <th>Drugs</th>
      <th>Transcriptomics</th>
      <th>Proteomics</th>
      <th>Mutations</th>
      <th>Copy Number</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="datasets/broad_sanger">Broad Sanger</a></td>
      <td>{{site.data.stats.broad_sanger.cancer_types}}</td>
      <td>{{site.data.stats.broad_sanger.samples}}</td>
      <td>{{site.data.stats.broad_sanger.drugs}}</td>
      <td>{{site.data.stats.broad_sanger.transcriptomics}}</td>
      <td>{{site.data.stats.broad_sanger.proteomics}}</td>
      <td>{{site.data.stats.broad_sanger.mutations}}</td>
      <td>{{site.data.stats.broad_sanger.copy_number}}</td>
    </tr>
    <tr>
      <td><a href="datasets/cptac">CPTAC</a></td>
      <td>{{site.data.stats.cptac.cancer_types}}</td>
      <td>{{site.data.stats.cptac.samples}}</td>
      <td>{{site.data.stats.cptac.drugs}}</td>
      <td>{{site.data.stats.cptac.transcriptomics}}</td>
      <td>{{site.data.stats.cptac.proteomics}}</td>
      <td>{{site.data.stats.cptac.mutations}}</td>
      <td>{{site.data.stats.cptac.copy_number}}</td>
    </tr>
    <tr>
      <td><a href="datasets/hcmi">HCMI</a></td>
      <td>{{site.data.stats.hcmi.cancer_types}}</td>
      <td>{{site.data.stats.hcmi.samples}}</td>
      <td>{{site.data.stats.hcmi.drugs}}</td>
      <td>{{site.data.stats.hcmi.transcriptomics}}</td>
      <td>{{site.data.stats.hcmi.proteomics}}</td>
      <td>{{site.data.stats.hcmi.mutations}}</td>
      <td>{{site.data.stats.hcmi.copy_number}}</td>
    </tr>
    <tr>
      <td><a href="datasets/beataml">BeatAML</a></td>
      <td>{{site.data.stats.beataml.cancer_types}}</td>
      <td>{{site.data.stats.beataml.samples}}</td>
      <td>{{site.data.stats.beataml.drugs}}</td>
      <td>{{site.data.stats.beataml.transcriptomics}}</td>
      <td>{{site.data.stats.beataml.proteomics}}</td>
      <td>{{site.data.stats.beataml.mutations}}</td>
      <td>{{site.data.stats.beataml.copy_number}}</td>
    </tr>
    <tr>
      <td><a href="datasets/mpnst">MPNST</a></td>
      <td>{{site.data.stats.mpnst.cancer_types}}</td>
      <td>{{site.data.stats.mpnst.samples}}</td>
      <td>{{site.data.stats.mpnst.drugs}}</td>
      <td>{{site.data.stats.mpnst.transcriptomics}}</td>
      <td>{{site.data.stats.mpnst.proteomics}}</td>
      <td>{{site.data.stats.mpnst.mutations}}</td>
      <td>{{site.data.stats.mpnst.copy_number}}</td>
    </tr>
  </tbody>
</table>




<!-- <div class="dataset-section">
    {% assign datasets = 'depmap,cptac,hcmi,beataml' | split: ',' %}
    {% for dataset in datasets %}
        <div class="dataset-container">
            <a href="datasets/{{ dataset }}" class="dataset-link">{{ dataset | capitalize }}</a>
            <div class="dataset-blurb">
                {% case dataset %}
                    {% when 'depmap' %}
                        <p>Cancer Types: {{ site.data.stats.depmap.depmaps }} </p>
                        <p>Samples: {{ site.data.stats.depmap.samples }} </p>
                        <p>Genes: {{ site.data.stats.depmap.genes }} </p>
                        <p>Drugs: {{ site.data.stats.depmap.drugs }} </p>
                        <span class="dot dot_transcriptomics"></span> 
                        <span class="dot dot_proteomics"></span> 
                        <span class="dot dot_mutations"></span> 
                        <span class="dot dot_copy_number"></span> 
                    {% when 'cptac' %}
                        <p>Cancer Types: {{ site.data.stats.cptac.cancer_types }} </p>
                        <p>Samples: {{ site.data.stats.cptac.samples }} </p>
                        <p>Genes: {{ site.data.stats.cptac.genes }} </p>
                        <p>Drugs: {{ site.data.stats.cptac.drugs }} </p>
                        <span class="dot dot_transcriptomics"></span> 
                        <span class="dot dot_proteomics"></span> 
                        <span class="dot dot_mutations"></span> 
                        <span class="dot dot_copy_number"></span> 
                    {% when 'hcmi' %}
                        <p>Cancer Types: {{ site.data.stats.hcmi.cancer_types }} </p>
                        <p>Samples: {{ site.data.stats.hcmi.samples }} </p>
                        <p>Genes: {{ site.data.stats.hcmi.genes }} </p>
                        <p>Drugs: {{ site.data.stats.hcmi.drugs }} </p>
                        <span class="dot dot_transcriptomics"></span> 
                        <span class="dot dot_proteomics"></span> 
                        <span class="dot dot_mutations"></span> 
                        <span class="dot dot_copy_number"></span> 
                    {% when 'beataml' %}
                        <p>Cancer Types: {{ site.data.stats.beataml.cancer_types }}</p>
                        <p>Samples: {{ site.data.stats.beataml.samples }} </p>
                        <p>Genes: {{ site.data.stats.beataml.genes }}</p>
                        <p>Drugs: {{ site.data.stats.beataml.drugs }}</p>
                        <span class="dot dot_transcriptomics"></span> 
                        <span class="dot dot_proteomics"></span> 
                {% endcase %}
            </div>
        </div>
    {% endfor %}

</div> -->


<!-- <div class="legend">
    <p>Transcriptomics<span class="dot dot_transcriptomics"></span></p>
    <p>Proteomics<span class="dot dot_proteomics"></span></p>
    <p>Mutations<span class="dot dot_mutations"></span></p>
    <p>Copy Number<span class="dot dot_copy_number"></span></p>
</div> -->


## Data Overview

<div class="flex-container"> 
    <div class="flex-item">
        <img src="{{ 'assets/stats/Fig0_Overview.png' | relative_url }}" alt="Summary 1" />
    </div>
    <div class="flex-item">
        <img src="{{ 'assets/stats/Fig6_Sample_Summary.png' | relative_url }}" alt="Summary 2" />
    </div>
</div>
