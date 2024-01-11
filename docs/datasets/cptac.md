---
layout: default
title: CoderData
---

<link rel="stylesheet" href="assets/css/style.css">


### Introduction
CoderData is a cancer benchmark data package developed in Python and R. 
There are two aspects of this package, the backend build section and the user facing python package.
The build section is a github workflow that generates four cancer datasets in a format that is easy for users and algorithms to ingest. 
The python package allows users to easily download the data, load it into python and reformat it as desired.

### CPTAC Summary

The Clinical Proteomic Tumor Analysis Consortium ([CPTAC](https://gdc.cancer.gov/about-gdc/contributed-genomic-data-cancer-research/clinical-proteomic-tumor-analysis-consortium-cptac#:~:text=The%20National%20Cancer%20Institute's%20Clinical,and%20genome%20analysis%2C%20or%20proteogenomics.
)) project is a collaborative network funded by the National Cancer Institute (NCI) focused on improving our understanding of cancer biology through the integration of transcriptomic, proteomic, and genomic data. 

{% if site.data.cptac_table %}
<table>
  {% for row in site.data.cptac_table %}
    {% if forloop.first %}
    <tr>
      {% for pair in row %}
        <th>{{ pair[0] }}</th>
      {% endfor %}
    </tr>
    {% endif %}

    {% tablerow pair in row %}
      {{ pair[1] }}
    {% endtablerow %}
  {% endfor %}
</table>
{% else %}
<p>CPTAC table is not available.</p>
{% endif %}

### Visualization

<div class="flex-container"> 
    <div class="flex-item">
        <img src="{{ 'assets/stats/Fig3_CPTAC.png' | relative_url }}" alt="CPTAC Figure" />

    </div>
    <div class="flex-item">
        <img src="{{ 'assets/stats/cptac_circos.png' | relative_url }}" alt="Cell Line Circos" />
    </div>
</div>
