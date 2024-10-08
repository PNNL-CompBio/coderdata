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

### HCMI Summary
Human Cancer Models Initiative (HCMI) data was collected though the National Cancer Institute (NCI) Genomic Data Commons (GDC) <a href="https://portal.gdc.cancer.gov/projects/HCMI-CMDC" target="_blank">Data Portal</a>.
This data encompasses numerous cancer types and includes cell line, organoid, and tumor data. Data includes the transcriptomics, somatic mutation, and copy number datasets.

{% if site.data.hcmi_table %}
<table>
  {% for row in site.data.hcmi_table %}
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
<p>HCMI table is not available.</p>
{% endif %}


### Visualization

<div class="flex-container"> 
    <div class="flex-item">
        <img src="{{ 'assets/stats/Fig1_HCMI.png' | relative_url }}" alt="HCMI Figure" />
    </div>
    <div class="flex-item">
        <img src="{{ 'assets/stats/hcmi_circos.png' | relative_url }}" alt="HCMI Circos" />
    </div>
</div>


