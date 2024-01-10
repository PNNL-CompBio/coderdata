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

### Cell Line Summary

The cell line datasets were collected from numerous resources such as the [LINCS project](https://lincsproject.org/), [DepMap](https://depmap.org/portal/), and the [Sanger Institute](https://www.sanger.ac.uk/).
This data will allow scientists to explore the drugs response for thousands of drugs across hundreds of cell lines.


{% if site.data.cell_line_table %}
<table>
  {% for row in site.data.cell_line_table %}
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




### Visualization

<div class="flex-container"> 
    <div class="flex-item">
        <img src="{{ 'assets/stats/Fig4_CCLE.png' | relative_url }}" alt="Cell Line Figure" />
    </div>
    {% endif %}
    <div class="flex-item">
        <img src="{{ 'assets/stats/cell_line_circos.png' | relative_url }}" alt="Cell Line Circos" />
    </div>
    {% endif %}
</div>