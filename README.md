## Cancer Omics Drug Experiment Response Dataset 

There is a recent explosion of deep learning algorithms that to tackle
the computational problem of predicting drug treatment outcome from
baseline molecular measurements. To support this,we have built a
Python package that enables access to and facile usage of cancer drug
sensitivity datsets for AI applications. 

This package collects diverse sets of paired molecular datasets with corresponding drug sensitivity data. All data here is reprocessed and standardized so it can be easily used as a benchmark dataset for the 
This repository leverages existing datasets to collect the data
required for deep learning model development. Since each deep learning model
requires distinct data capabilities, the goal of this repository is to
collect and format all data into a schema that can be leveraged for
existing models.

![Coderdata Motivation](coderdata_overview.jpg?raw=true "Motivation behind
coderdata develompent")

## Installation
To install the CoderData Python package:
```
pip install coderdata
```

## Usage
The Python package is designed to facilitate the training and
validating of computational models that predict drug
response. Currently the package supports the following commands:

1. `list`: Lists names of datasets to download. This depends on what
   datasets have been in the main build. 
2. `download`: Downloads the dataset by name. Case insensitive, but
   should contain full name of dataset returned from the `list` command. 
2. `load`: This return a `dataset` object that houses all the
   data. This object has the following functions:
   1. `train_test_validate`: Splits the object into 
   2. `types`: Each datasaet has different types of data included in
      it. For all possible data types see the
      [schema](schema/README.md). These can include:
      - transcriptomics
      - mutations
      - copy_numbers
      - proteomics
      - experiments
      - combinations
      - drugs
      - genes
      - samples
   3. `format`: Format performs data type-specific formatting. First
      argument is name of data type, next arguments are
      data-type-specific, last argument is `use_polars` to return a
      polars instead of a pandas data frame.
      - transcriptomics: `ds.format('transcriptomics')` returns a
        pandas or polars data frame with each row representing a gene
        and each column representing a sample.
      - mutations:
        `ds.format('mutations',['Frame_Shift_Del','Frame_Shift_Ins','Missense_Muation','Start_Codon_SNP'])`
        will return a binary  matrix with rows representing genes and
        columns representing samples, and a `1` value if there there
        is a mutation in given gene/sample that falls into the class
        provided by the second argument. 
      - copy_numbers: `ds.format('copy_number','copy_number')` returns a
        pandas or polars data frame with each row representing a gene
        and each column representing the average copy_number value for
        each gene. If the second argument is `copy_call` the data
        frame values are a discrete measurement of copy number
        provided by the schema. 
      - proteomics: `ds.format('proteomics')` returns a
        pandas or polars data frame with each row representing a gene
        and each column representing a sample.
      - experiments: `ds.format('experiments', 'fit_auc')` returns a
      matrix with drugs represented by rows and samples represneted by
      columns the numeric values represent the measurement provided by
      the second argument. 
      - combinations: `ds.format('combinations')` returns ???
      - drugs: Not sure what to return here - just table? how about descriptors?
      - samples: just return table?
   4. `save`: saves the object to a file. 

## Additional documentation
For the access to the latest version of CoderData, please visit our
[documentation site](https://pnnl-compbio.github.io/coderdata/) which provides access to Figshare and
instructions for using the Python package to download the data.


