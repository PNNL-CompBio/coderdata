---
layout: default
title: CoderData
---

<link rel="stylesheet" href="assets/css/style.css">

<!-- ## Usage of CoderData -->

## Introduction
CoderData is a comprehensive package designed for handling cancer benchmark data in Python.  
It offers functionalities to download datasets, load them into Python environments, and reformat them according to user needs.

## Installation
`coderdata` requires `python>=3.9` to be installed. The installed version can be checked via
```shell
$ python --version
Python 3.13.1
```
If a Python version older that 3.9 is installed please referr to the instruction at [python.org](https://www.python.org/about/gettingstarted/#installing) on how to install / update Python.

The preferred way to install `coderdata` is via `pip`. Executing the command below will install the most recent published version of `coderdata` including all required dependencies.
```shell
$ pip install coderdata
```

To check if the package has been sucessfully installed open an interactive python termial and import the package. See an example of what to expect below.
```python
>>> import coderdata as cd
>>> cd.__version__
'0.1.40'
```

## Usage
The primary way to interact with coderdata is through the `coderdata` API. Additionally a command line interface with limited functionality (primarily to download data) is also available.

### CLI
Invoking `coderdata` from the command line will by default print a help / usage message and exit (see below):
```sh
$ coderdata
usage: coderdata [-h] [-l | -v] {download} ...

options:
  -h, --help     show this help message and exit
  -l, --list     prints list of available datasets and exits program.
  -v, --version  prints the versions of the coderdata API and dataset and exits the program

commands:
  {download}
    download     subroutine to download datasets. See "coderdata download -h" for more options.
```

The primary use case of the CLI is to retrieve dataset from the repository. This can be done by invoking the `download` routine of `coderdata`. Without defining a specific dataset the whole repository will be downloaded:
```sh
$ coderdata download
Downloaded 'https://ndownloader.figshare.com/files/48032953' to '/tmp/beataml_drugs.tsv.gz'
Downloaded 'https://ndownloader.figshare.com/files/48032962' to '/tmp/mpnst_drugs.tsv.gz'
...
```

Downloading a specific dataset can be achieve by passing the `-n/--name` argument to the `download` routine:
```sh
$ coderdata download --name beataml
Downloaded 'https://ndownloader.figshare.com/files/48032953' to 'beataml_drugs.tsv.gz'
Downloaded 'https://ndownloader.figshare.com/files/48032959' to 'beataml_samples.csv'
...
```

A full list of available arguments of the `download` function including a short explanation can retrieved via the command shown below:
```sh
$ coderdata download -h
usage: coderdata download [-h] [-n DATASET_NAME] [-p LOCAL_PATH] [-o]

options:
  -h, --help            show this help message and exit
  -n, --name DATASET_NAME
                        name of the dataset to download (e.g., "beataml"). Alternatively, "all" will download the full repository of coderdata datasets. See "coderdata --list" for a
                        complete list of available datasets. Defaults to "all"
  -p, --local_path LOCAL_PATH
                        defines the folder the datasets should be stored in. Defaults to the current working directory if omitted.
  -o, --overwrite       allow dataset files to be overwritten if they already exist.
```

Additionally to the `download` functionality, the CLI currently supports displaying basic information such as the version numbers of the package and the dataset (see example call below)
```sh
$ coderdata --version
package version: 0.1.40
dataset version: 0.1.4
```
as well as listing the dataset that are available for download (example output below)
```sh
$ coderdata --list

Available datasets
------------------

beataml: Beat acute myeloid leukemia (BeatAML) focuses on acute myeloid leukemia tumor data. Data includes drug response, proteomics, and transcriptomics datasets.
cptac: The Clinical Proteomic Tumor Analysis Consortium (CPTAC) project is a collaborative network funded by the National Cancer Institute (NCI) focused on improving our understanding of cancer biology through the integration of transcriptomic, proteomic, and genomic data.
hcmi: Human Cancer Models Initiative (HCMI) encompasses numerous cancer types and includes cell line, organoid, and tumor data. Data includes the transcriptomics, somatic mutation, and copy number datasets.
mpnst: Malignant Peripheral Nerve Sheath Tumor is a rare, agressive sarcoma that affects peripheral nerves throughout the body.

------------------

To download individual datasets run "coderdata download -name DATASET_NAME" where "DATASET_NAME" is for example "beataml".
```

### API

#### Downloading data
Using the `coderdata` API, the download process is handled through the `download` function in the downloader module.
```python
>>> import coderdata as cd
>>> cd.download(name='beataml')
Downloaded 'https://ndownloader.figshare.com/files/48032953' to 'beataml_drugs.tsv.gz'
Downloaded 'https://ndownloader.figshare.com/files/48032959' to 'beataml_samples.csv'
Downloaded 'https://ndownloader.figshare.com/files/48032965' to 'beataml_mutations.csv.gz'
Downloaded 'https://ndownloader.figshare.com/files/48032968' to 'beataml_proteomics.csv.gz'
Downloaded 'https://ndownloader.figshare.com/files/48032974' to 'beataml_experiments.tsv.gz'
Downloaded 'https://ndownloader.figshare.com/files/48033052' to 'beataml_transcriptomics.csv.gz'
Downloaded 'https://ndownloader.figshare.com/files/48033058' to 'beataml_drug_descriptors.tsv.gz'
```
As with the CLI download functionality, the local path where to store the downloaded files, as well as a flag the defines whether existing files should be overwritten can be defined in the `download()` function. For example the function call below will download all 'BeatAML' related datasets to the local path `/tmp/coderdata/` and will overwrite files if they already exist.
```python
>>> cd.download(name='beataml', local_path='/tmp/coderdata/', exist_ok=True)
```
Note that if `exist_ok==False` (the default if omitted) and a downloaded file already exists a warning will be given and the file won't be stored. Finally, if all datasets should be downloaded the `name` argument can manually set to `name='all'` or omitted all together as the `name` defaults to `'all'`.

#### The `Dataset` object

The `Dataset` object is the central data structure in CoderData. It automatically initializes attributes for each dataset type like tumor samples, drug response data, as well as associated omics data like proteomics. Each datatype in a `Dataset` is internally stored in a [`pandas.DataFrame`](https://pandas.pydata.org/docs/reference/frame.html).

##### Loading data into a `Dataset` object
The code snippet will load the [previously downloaded](#downloading-data) 'BeatAML' dataset into a `Dataset` object called `beataml`.

```python
>>> beataml = cd.load(name='beataml', local_path='/tmp/coderdata')
Importing raw data ...
Importing 'transcriptomics' from /tmp/coderdata/beataml_transcriptomics.csv.gz ... DONE
Importing 'drugs' from /tmp/coderdata/beataml_drugs.tsv.gz ... DONE
Importing 'proteomics' from /tmp/coderdata/beataml_proteomics.csv.gz ... DONE
Importing 'drug_descriptors' from /tmp/coderdata/beataml_drug_descriptors.tsv.gz ... DONE
Importing 'mutations' from /tmp/coderdata/beataml_mutations.csv.gz ... DONE
Importing 'samples' from /tmp/coderdata/beataml_samples.csv ... DONE
Importing 'experiments' from /tmp/coderdata/beataml_experiments.tsv.gz ... DONE
Importing raw data ... DONE
```

Additionally, the `load()` function also allows for loading data from a previously pickled `Dataset` object (see [Saving manipulated `Dataset` objects](#saving-manipulated-dataset-objects)).

##### Displaying the datatypes in a `Dataset` object

The data types associated with a dataset can be displayed via the `Dataset.types()` function. The function will return a simple list of available datatypes.
```python
>>> beataml.types()
['transcriptomics', 'proteomics', 'mutations', 'samples', 'drugs', 'experiments']
```
Individual datatypes can be adressed and manipulated by subscripting the dataset. For example extracting the underlying `pandas.DataFrame` that contains drug response values for 'BeatAML' can be done via the command below:
```python
>>> beataml.experiments
         source  improve_sample_id improve_drug_id    study  time time_unit dose_response_metric  dose_response_value
0       synapse               3907       SMI_11123  BeatAML    72       hrs              fit_auc               0.0564
1       synapse               3907       SMI_11211  BeatAML    72       hrs              fit_auc               0.9621
2       synapse               3907       SMI_12192  BeatAML    72       hrs              fit_auc               0.1691
3       synapse               3907       SMI_12254  BeatAML    72       hrs              fit_auc               0.4245
4       synapse               3907       SMI_12469  BeatAML    72       hrs              fit_auc               0.7397
...         ...                ...             ...      ...   ...       ...                  ...                  ...
233775  synapse               3626        SMI_7110  BeatAML    72       hrs                  dss               0.0000
233776  synapse               3626        SMI_7590  BeatAML    72       hrs                  dss               0.0000
233777  synapse               3626        SMI_8159  BeatAML    72       hrs                  dss               0.1946
233778  synapse               3626        SMI_8724  BeatAML    72       hrs                  dss               0.0000
233779  synapse               3626         SMI_987  BeatAML    72       hrs                  dss               0.7165

[233780 rows x 8 columns]
```

##### Reformatting and exporting datatypes

Internally all data is stored in long format. If different formats are needed for further analysis or as input for the training of machine learning models, the `Dataset.format(data_type, **kwargs)` function is able to return individual data types in altered formats.

For example the drug response data can be reformatted into wide format via the following command:
```python
>>> beataml.format(data_type='experiments', shape='wide', metrics=['fit_auc', 'dss'])
        source  improve_sample_id improve_drug_id    study  time time_unit     dss  fit_auc
0      synapse               3190       SMI_11123  BeatAML    72       hrs  0.4244   0.5447
1      synapse               3190       SMI_12192  BeatAML    72       hrs  0.2782   0.4848
2      synapse               3190       SMI_12254  BeatAML    72       hrs  0.0000   0.5872
3      synapse               3190       SMI_12469  BeatAML    72       hrs  0.2973   0.4435
4      synapse               3190       SMI_12953  BeatAML    72       hrs  0.0000   0.5566
...        ...                ...             ...      ...   ...       ...     ...      ...
23373  synapse               3916        SMI_7590  BeatAML    72       hrs  0.4537   0.5689
23374  synapse               3916        SMI_8063  BeatAML    72       hrs  0.0000   0.5640
23375  synapse               3916        SMI_8159  BeatAML    72       hrs  0.0000   0.5340
23376  synapse               3916        SMI_8724  BeatAML    72       hrs  0.7033   0.7172
23377  synapse               3916         SMI_987  BeatAML    72       hrs  0.0000   0.4842
```
Note that the `Dataset.format(data_type, **kwargs)` function behaves slightly different for different `data_type` values. For example for `data_type='experiments'` accepted keyword arguments are `shape` & `metrics`. `shape` defines which format the resulting `pandas.DataFrame` should be in (e.g. `long`, `wide` or `matrix`). `metrics` defines the drug response metrics that should be filtered for.

A full list of parameters for the individual data types can be found below:
- `Dataset.format(data_type='transcriptomics')` returns a `matrix` like `pandas.DataFrame` where each cell contains the measured transcriptomics value for a gene (row - `entrez_id`) in a specific cancer sample (column - `improve_sample_id`).
- `Dataset.format(data_type='mutations', mutation_type=...)` will return a binary `matrix` like `pandas.DataFrame` with rows representing genes and columns representing samples. `mutation_type` can be any of the recoreded mutation types available (e.g. `'Frame_Shift_Del'`,`'Frame_Shift_Ins'`,`'Missense_Muation'` or `'Start_Codon_SNP'` among others). Cells contain the value of `1` if a mutation in given gene/sample falls into the category defined by `mutation_type`.
- `Dataset.format(data_type='copy_number', copy_call=False)` returns a `matrix` like `pandas.DataFrame` where cells report the `mean` copy number value for each combination of gene (row - `entrez_id`) and cancer sample (column - `improve_sample_id`). If `copy_call=True` cells report the discretized measurement ('deep del', 'het loss', 'diploid', 'gain', 'amp') of copy number provided by the schema.
- `Dataset.format('data_type=proteomics')` returns a `matrix` like `pandas.DataFrame` where each cell contains the measured proteomics value for a gene (row - `entrez_id`) in a specific cancer sample (column - `improve_sample_id`).
- `Dataset.format(data_type='experiments', shape=..., metrics=...)`, returns a formatted `pandas.DataFrame` according to defined `shape` (`shape` can be of values `'long'`, `'wide'` and `'matrix'`). `metrics` further defines which drug response metrics the resulting output `DataFrame` should be filtederd for. Examples are `'fit_auc'`, `'fit_ec50` or `'dss'`. If `shape=wide`, a list can be passed to `metric` containing more than one value.
- `Dataset.format(data_type='drug_descriptor', shape=..., drug_descriptor_type=...)` returns a `pandas.DataFrame` formatted either in `long` or `wide` (depending on the `shape` argument). `drug_descriptor_type` can be defined as a list of desired `structural_descriptors` in conjunction with `shape=wide`, to limit the resulting `DataFrame` to only list the desired `structual_descriptors` as columns.
- `Dataset.format(data_type='drugs')` is equal to `Dataset.drugs`. It returns the underlying `pandas.DataFrame` containing the drug information.
- `Dataset.format(data_type='genes')` is equal to `Dataset.genes`. It returns the underlying `pandas.DataFrame` containing the gene information.
- `Dataset.format(data_type='samples')` is equal to `Dataset.samples`. It returns the underlying `pandas.DataFrame` containing the cancer sample data information.

##### Creating training / testing and validation splits with `coderdata`

Using the `Dataset.train_test_validate()` function the dataset can be split into trining, testing and validation sets. The function will return a `Split` object (a python `@dataclass`) that contains three `Dataset` objects that can be adressed and retrieved by subscripting with eiter `Split.train`, `Split.test` or `Split.validate`. 

```python
>>> split = beataml.train_test_validate()
>>> split.train.experiments.shape
(187020, 8)
>>> split.test.experiments.shape
(23380, 8)
>>> split.validate.experiments.shape
(23380, 8)
```

By default the returned splits will be `mixed-set` (drugs and cancer samples can appear in all three folds), with a ratio of 8:1:1, no stratification and no set random state (seed). This behaviour can be changed by passing `split_type`, `ratio`, `stratified_by` and `random_state` to the function. 

`split_type` can be either `'mixed-set'`, `'drug-blind'` or `'drug-blind'`:
- `mixed-set`: Splits randomly independent of drug / cancer association of the samples. Individual drugs or cancer types can appear in all three splits
- `drug-blind`: Splits according to drug association. Any sample associated with a drug will be unique to one of the splits. For example samples with association to drug A will only be present in the train split, but never in test or validate.
- `cancer-blind`: Splits according to cancer association. Equivalent to drug-blind, except cancer types will be unique to splits.

`ratio` can be used to adjust the split ratios using a 3 item tuple containing integers. For example `ratio=(5:3:2)` would result in a split where train, test and validate contain roughly 50%, 30% and 20% of the original data respectively.

`random_state` defines a seed values for the random number generator. Defining a `random_state` will guarantee reproducability as two runs with the same `random_state` will result in the same splits.

`stratify_by` Defines if the training, testing and validation sets should be stratified. Stratification tries to maintain a similar distribution of feature classes across different splits. For example assuming a drug respones value threshold that defines positive and negative classes (e.g. reduced vs. no change in cancer cell viability) the splitting algorithm could attempt to assign the same amount of positive class instances as negative class instances to each split. Stratification is performed by `drug_response_value`. Any value other than `None` indicates stratification and defines which `drug_response_value` should be used as basis for the stratification. `None` indicates that no stratfication should be performed. Which type of stratification should be performe can further be customized with keyword arguments (`thresh`, `num_classes`, `quantiles`).

An example call to create a 70/20/10 drug-blind split that is stratified by `fit_auc` could look like this:
```python
>>> split = beataml.train_test_validate(
...     split_type='drug-blind',
...     ratio=[7,2,1],
...     random_state=42,
...     stratify_by='fit_auc',
...     thresh=0.8
...     )
>>> split.train.experiments.shape
(154840, 8)
>>> split.test.experiments.shape
(65750, 8)
>>> split.validate.experiments.shape
(13190, 8)
```

##### Saving manipulated `Dataset` objects (e.g. saving splits)
In order to save a `Dataset` for later use, the `Dataset.save()` function can be used.

```python
>>> split.train.save(path='/tmp/coderdata/beataml_train.pickle')
>>> split.test.save(path='/tmp/coderdata/beataml_test.pickle')
>>> split.validate.save(path='/tmp/coderdata/beataml_validate.pickle')
```

This function can be used to either save the individual splits (as demonstrated above), or raw `Dataset` that was the basis for the splits for example if any modifications of the dataset were performed.

To reload the splits (or the full dataset) the `coderdata.load()` function (see also [Loading data into a `Dataset` object](#loading-data-into-a-dataset-object)) can be used. To load a pickled `Dataset`, the argument `from_pickle=True` must be passed to the function:

```python
>>> beataml_train = cd.load('beataml_train', local_path='/tmp/coderdata/', from_pickle=True)
Importing pickled data ... DONE
```
Note that only individual splits (e.g. only train) can be saved and loaded and not the full `Split` object.

## Conclusion
CoderData provides a robust and flexible way to work with cancer benchmark data.   
By using these functionalities, researchers and data scientists can easily manipulate and analyze complex datasets in their Python environments
