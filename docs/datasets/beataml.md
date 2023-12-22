---
layout: default
title: CoderData
pdf_exists: true
png_exists: false
---

<link rel="stylesheet" href="assets/css/style.css">


### Introduction
CoderData is a cancer benchmark data package developed in Python and R. 
There are two aspects of this package, the backend build section and the user facing python package.
The build section is a github workflow that generates four cancer datasets in a format that is easy for users and algorithms to ingest. 
The python package allows users to easily download the data, load it into python and reformat it as desired.

### BeatAML Summary

Beat acute myeloid leukemia (BeatAML) data was collected though the [GitHub](https://biodev.github.io/BeatAML2/) and [Synapse](https://www.synapse.org/#!Synapse:syn24171150).
This data focuses on acute myeloid leukemia tumor data. Data includes drug response, proteomics, and transcriptomics datasets.


{% if site.data.beataml_table %}
<table>
  {% for row in site.data.beataml_table %}
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
<p>Cell line table is not available.</p>
{% endif %}
{% else %}
<p>BeatAML table is not available.</p>
{% endif %}


### Visualization

<div class="flex-container"> 
    {% if page.pdf_exists %}
    <div class="flex-item">
        <embed src="{{ 'assets/stats/Fig2_BeatAML.pdf' | relative_url }}" type="application/pdf" />
    </div>
    {% endif %}
    {% if page.png_exists %}
    <div class="flex-item">
        <img src="{{ 'assets/stats/beataml_circos.png' | relative_url }}" alt="BeatAML Circos" />
    </div>
    {% endif %}
</div>