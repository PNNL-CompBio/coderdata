---
layout: default
title: CoderData
---

<link rel="stylesheet" href="assets/css/style.css">

## Contribute to CoderData

CoderData is a data assembly pipeline that pulls from
original data sources of drug sensitivity and omics datasets and
assembles them so they can be integrated into a python package for
AI/ML integration. 

CoderData is indeed a work in progress. If you have specific requests
or bugs, please file an issue on our [GitHub
page](https://github.com/PNNL-CompBio/coderdata/issues) and we will
begin a conversation about details and how to fix the issue. If you
would like to create a new feature to address the issue, you are welcome to fork the
repository and create a pull request to discuss it in more
detail. These will be triaged by the CoderData team as they are received.

The rest of this document is focused on how to contribute to and
augment CoderData, either for use by the community or your own
purpopses. 

### CoderData build process

To build your own internal Coderdata dataset, or to augment it, it is
important to understand how the package is built. 

The build process is managed in the [build
directory](https://github.com/PNNL-CompBio/coderdata/tree/main/build)
primarily by the [`build_all.py` script](https://github.com/PNNL-CompBio/coderdata/blob/main/build/build_all.py). This script calls the
[`build_dataset.py` script](https://github.com/PNNL-CompBio/coderdata/blob/main/build/build_dataset.py) for each dataset in CoderData in
order. Because our sample and drug identifiers are unique, we must
finish the generation of one dataset before we move to the next. This
process is depicted below. 


![Coderdata Build](coderDataBuild.jpg?raw=true "Modular build
process")

The build process is slow, partially due to our querying of PubChem,
and also because of our extensive curve fitting. However, it can be
run locally so that you can still leverage the Python package
functionality with your own datasets.

If you want to add a new dataset, you must create a docker image that
contains all the scripts to pull the data and reformat it into our
[LinkML Schema](https://github.com/PNNL-CompBio/coderdata/blob/main/schema/coderdata.yaml). Once complete, you can modify `build_dataset.py` to
call your Docker image and associated scripts. More details are below.

## Adding your own dataset

To add your own data, you must add a Docker image with the following
constraints:

1. Be named `Dockerfile.[dataset_name]` and reside in the
   `/build/docker` directory
2. Possess scripts called `build_omics.sh`, `build_samples.sh`,
   `build_drugs.sh`,  `build_exp.sh` , and if needed, a
   `build_misc.sh`. These will all be called directly by
   `build_dataset.py`. 
3. Create tables that mirror the schema described by the [LinkML YAML
   file](https://github.com/PNNL-CompBio/coderdata/blob/main/schema/coderdata.yaml).
   
Files are generated in the following order as described below.


### Sample generation

The first step of any dataset build is to create a unique set of
sample identifies and store them in a `[dataset_name]_samples.csv`
file. We recommend following these steps:

1. Build a python script that pulls the sample identifier information
   from a stable repository and generates Improve identiefiers for
   each sample while also ensuring that no sample identifiers are
   clashing with prior samples. Examples can be found [here](https://github.com/PNNL-CompBio/coderdata/blob/main/build/mpnst/00_sample_gen.R) and [here](https://github.com/PNNL-CompBio/coderdata/blob/main/build/broad_sanger/01-broadSangerSamples.R). If
   you are using the Genomic Data Commons, you can leverage our
   existing scripts [here](https://github.com/PNNL-CompBio/coderdata/blob/main/build/hcmi/01-createHCMISamplesFile.py). 
2. Create a `build_samples.sh` script that calls your script with an
   existing sample file as the first argument. 
3. Test the `build_samples.sh` script with a [test sample
   file](https://github.com/PNNL-CompBio/coderdata/blob/main/build/build_test/test_samples.csv). 
4. Validate the file with the [linkML validation tool](https://linkml.io/linkml/cli/validate.html) and our
   [schema file](https://github.com/PNNL-CompBio/coderdata/blob/main/schema/coderdata.yaml). 

### Omics data generation

The overall omics generation process is the same as the samples, with
a few caveats. 

1. Build a python script that maps the omics data and gene data to the
   standardized identifiers and aligns them to the schema. 
   pulls the sample identifier information
   from a stable repository and generates Improve identiefiers for
   each sample while also ensuring that no sample identifiers are
   clashing with prior samples. Examples can be found here and here. If
   you are using the Genomic Data Commons, you can leverage our
   existing scripts here. For each type of omics data (see below), a
   single file is created.It might take more than one script, but you
   can combine those in step 2. 
2. Create a `build_omics.sh` script that calls your script with the
   `genes.csv` file as the first argument and `[dataset_name]_samples.csv` file as second
   argument. 
3. Test the `build_omics.sh` script with your sample file and [test genes
   file](https://github.com/PNNL-CompBio/coderdata/blob/main/build/build_test/test_genes.csv). 
4. Validate the files generated with the [linkML validation tool](https://linkml.io/linkml/cli/validate.html) and our
   [schema file](https://github.com/PNNL-CompBio/coderdata/blob/main/schema/coderdata.yaml). 

The precise data files have varying standards, as described below:

- *Mutation data:* In addition to matching gene identifiers each gene
mutation should be mapped to a specific schema of variations.  The
list of allowed variations can be found [in our linkML
file](https://github.com/PNNL-CompBio/coderdata/blob/8000968dc5f19fbb986a700862c5035a0230b656/schema/coderdata.yaml#L247). 
- *Transcriptomic data:* Transcript data is mapped to the same gene
  identifiers and sample sbut is convered to transcripts per million,
  or TPM. 
- *Copy number data:*  Copy number is assumed to be a value
  represneting the number of copies of that gene in a particular
  sample. 2 is assumed to be diploid. 
- *Proteomic data:* Proteomic measurements are generally log ratio
  values of the abundance measurements normalized to an internal
  control. 
  
The resulting files are then stored as [dataset_name]_[datatype].csv. 

### Drug data generation

The drug generation process can be slow depending on how many drugs
require querying from PubChem. However, with the use of an existing
drug data file, it's possible to shorten this process.

1. Build a python script that maps the drug information to SMILES
   String and IMPROVE identifier. All drugs are given an Improve
   identifier based on the canonical SMILES string to ensure that 
   each drug has a unique structure to be used in the modeling
   process. To standardize this we encourage using
our [standard drug lookup
  script](http://github.com/pnnl-compbio/coderdata/tree/main/build/utils/pubchem_retrieval.py)
  that retrieves drug structure and information by name or
   identifier. [This file of NCI60
   drugs](https://github.com/PNNL-CompBio/coderdata/blob/main/build/broad_sanger/03a-nci60Drugs.py)
   is our most comprehensive script as it pulls over 50k drugs
1a. In cases where the dose and response values are not available, you
   can use the published AUC values instead, and use the
   `published_auc` as the `drug_response_metric` value in the table. 
2. Create a  `build_drugs.sh` script that takes as its first argument
an existing drug file and calls the script created in step 1
above. Once the drugs for a dataset are retrieved, we have a second utility
script that [builds the drug descriptor table](https://github.com/PNNL-CompBio/coderdata/blob/cbf017326b83771c55f12317189f4b2dbd9d900a/schema/coderdata.yaml#L94). Add this to the
shell script to generate the drug descriptor file.
3. Test the `build_drugs.sh` script with the [test drugs
   file] (TBD). 
4. Validate the files generated with the [linkML validation tool](https://linkml.io/linkml/cli/validate.html) and our
   [schema file](https://github.com/PNNL-CompBio/coderdata/blob/main/schema/coderdata.yaml). 

 The resulting files should be `[dataset_name]_drugs.tsv` and
 `[dataset_name]_drug_descriptors.tsv`. 
 
### Experiment data generation

The experiment file maps the sample information to the drugs of
interest with various drug response metrics. The experiment data
varies based on the type of system:
- Cell line and organoid data use the [drug curve fitting
  tool](http://github.com/pnnl-compbio/coderdata/tree/main/build/utils/fit_curve.py)
  that maps doses of drugs (in Moles) to drug response measurements
  (in percent) to a variety of curve fitting metrics described in our
  [schema file](https://github.com/PNNL-CompBio/coderdata/blob/8000968dc5f19fbb986a700862c5035a0230b656/schema/coderdata.yaml#L200). 
- Patient derived xenografts require an alternate script that [creates
  PDX-speciic metrics](https://github.com/PNNL-CompBio/coderdata/blob/main/build/utils/calc_pdx_metrics.py). 
  
Otherwise the steps for building an experiment file are similar:
1. Build a python script that maps the drug information and sample
information to the DOSE and GROWTH values, then calls the appropriate
curve fitting tool described above. 
2. Create a  `build_exp.sh` script that takes as its first argument
the samples file and the second argument the drug file. 
3. Test the `build_exp.sh` script with the drug and samples files. 
4. Validate the files generated with the [linkML validation tool](https://linkml.io/linkml/cli/validate.html) and our
   [schema file](https://github.com/PNNL-CompBio/coderdata/blob/main/schema/coderdata.yaml). 

  
### Dockerize and test

All scripts described above go into a single directory with the name
of the dataset under the [build](http://github.com/pnnl-compbio/coderdata/tree/main/build) directory, with instructions to add everything in the [docker](http://github.com/pnnl-compbio/coderdata/tree/main/build/docker)
directory. Make sure to include any requirements for building in the
folder and docker image as well. 

Once the Dockerfile builds and runs, you can modify the
`build_dataset.py` script so that it runs and validates. 

Check out examples! We have numerous Docker files in our
[Dockerfile
directory](http://github.com/pnnl-compbio/coderdata/tree/main/build/docker),
and multiple datasets in our [build
directory](http://github.com/pnnl-compbio/coderdata/tree/main/build). 


---  
Your contributions are essential to the growth and improvement of CoderData. We look forward to collaborating with you!  
