[![Continuous Integration](https://github.com/PNNL-CompBio/coderdata/actions/workflows/main.yml/badge.svg?branch=builder_branch_JJ&event=push)](https://github.com/PNNL-CompBio/coderdata/actions/workflows/main.yml)

## Cancer Omics Drug Experiment Response Dataset 

There is a recent explosion of deep learning algorithms that to tackle the computational problem of predicting drug treatment outcome from baseline molecular measurements. To support this,we have built a benchmark dataset that harmonizes diverse datasets to better assess algorithm performance.

This package collects diverse sets of paired molecular datasets with corresponding drug sensitivity data. All data here is reprocessed and standardized so it can be easily used as a benchmark dataset for the 
This repository leverages existing datasets to collect the data
required for deep learning model development. Since each deep learning model
requires distinct data capabilities, the goal of this repository is to
collect and format all data into a schema that can be leveraged for
existing models.

The goal of this repository is two-fold: First, it aims to collate and
standardize the data for the broader community. This requires
running a series of scripts to build and append to a standardized data
model. Second, it has a series of scripts that pull from the data
model to create model-specific data files that can be run by the data
infrastructure. 

# coderdata Data Model

The coderdata schema is maintained in [LinkML](schema/coderdata.yaml)
and can be udpated via a commit to the repository. For more details,
please see the [schema description](schema/README.md).

## Building the data package

The data package is currently assembled via continuous automation,
