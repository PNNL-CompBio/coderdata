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

<div class="flex-container"> 
    <div class="flex-item">
        <embed src="{{ 'assets/stats/Fig4_CCLE.pdf' | relative_url }}" type="application/pdf" />
    </div>
    <div class="flex-item">
        <img src="{{ 'assets/stats/cell_line_circos.png' | relative_url }}" alt="Cell Line Circos" />
    </div>
</div>