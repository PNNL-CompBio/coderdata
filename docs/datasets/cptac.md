---
layout: default
title: CoderData
pdf_exists: true
png_exists: true
---

<link rel="stylesheet" href="assets/css/style.css">


### Introduction
CoderData is a cancer benchmark data package developed in Python and R. 
There are two aspects of this package, the backend build section and the user facing python package.
The build section is a github workflow that generates four cancer datasets in a format that is easy for users and algorithms to ingest. 
The python package allows users to easily download the data, load it into python and reformat it as desired.

### CPTAC Summary


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
    {% if page.pdf_exists %}
    <div class="flex-item">
        <embed src="{{ 'assets/stats/Fig3_CPTAC.pdf' | relative_url }}" type="application/pdf" />
    </div>
    {% endif %}
    {% if page.png_exists %}
    <div class="flex-item">
        <img src="{{ 'assets/stats/cptac_circos.png' | relative_url }}" alt="Cell Line Circos" />
    </div>
    {% endif %}
</div>
