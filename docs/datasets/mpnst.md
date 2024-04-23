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

### MPNST Summary
Malignant Peripheral Nerve Sheath Tumor is a rare, agressive sarcoma that affects peripheral nerves throughout the body. Data was collected from the <a href="https://nf.synapse.org/" target="_blank">NF Data Portal</a> and includes transcriptomics, mutation, copy number, and drug response data.



{% if site.data.mpnst_table %}
<table>
  {% for row in site.data.mpnst_table %}
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
<p>MPNST table is not available.</p>
{% endif %}


### Visualization

<div class="flex-container"> 
    <div class="flex-item">
        <img src="{{ 'assets/stats/Fig5_MPNST.png' | relative_url }}" alt="MPNST Figure" />
    </div>
    <div class="flex-item">
        <img src="{{ 'assets/stats/mpnst_circos.png' | relative_url }}" alt="MPNST Circos" />
    </div>
</div>


