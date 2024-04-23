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

### Broad Sanger Summary

The Broad Sanger datasets were collected from numerous resources such as the <a href="https://lincsproject.org/" target="_blank">LINCS project</a>, <a href="https://depmap.org/portal/" target="_blank">DepMap</a>, and the <a href="https://www.sanger.ac.uk/" target="_blank">Sanger Institute</a>.
This data will allow scientists to explore the drugs response for thousands of drugs across hundreds of cell lines.



{% if site.data.broad_sanger_table %}
<table>
  {% for row in site.data.broad_sanger_table %}
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
<p>Borad Sanger table is not available.</p>
{% endif %}




### Visualization

<div class="flex-container"> 
    <div class="flex-item">
        <img src="{{ 'assets/stats/Fig4_Broad_Sanger.png' | relative_url }}" alt="Broad Sanger Figure" />
    </div>
    <div class="flex-item">
        <img src="{{ 'assets/stats/broad_sanger_circos.png' | relative_url }}" alt="Broad Sanger Circos" />
    </div>
</div>