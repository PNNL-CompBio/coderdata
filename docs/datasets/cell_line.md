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

### Cell Line Summary

{% if site.data.cell_line_table %}
<table>
  <thead>
    <tr>
      {% assign first_row = site.data.cell_line_table[0] %}
      {% for cell in first_row %}
      <th>{{ cell }}</th> 
      {% endfor %}
    </tr>
  </thead>
  <tbody>
    {% for row in site.data.cell_line_table offset:1 %} 
    <tr>
      {% for cell in row %}
      <td>{{ cell }}</td> 
      {% endfor %}
    </tr>
    {% endfor %}
  </tbody>
</table>
{% else %}
<p>Cell line table is not available.</p>
{% endif %}


### Visualization

<div class="flex-container"> 
    {% if page.pdf_exists %}
    <div class="flex-item">
        <embed src="{{ 'assets/stats/Fig4_CCLE.pdf' | relative_url }}" type="application/pdf" />
    </div>
    {% endif %}
    {% if page.png_exists %}
    <div class="flex-item">
        <img src="{{ 'assets/stats/cell_line_circos.png' | relative_url }}" alt="Cell Line Circos" />
    </div>
    {% endif %}
</div>