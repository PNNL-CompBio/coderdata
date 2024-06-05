---
layout: default
title: CoderData
---

<link rel="stylesheet" href="assets/css/style.css">

## Contribute to CoderData

CoderData is designed to be a customizable resources that can be
altered and appended for your own needs. 

## Issues with current version

CoderData is indeed a work in progress. If you have specific requests
or bugs, please file an issue on our [GitHub
page](https://github.com/PNNL-CompBio/coderdata/issues) and we will
begin a conversation about details and how to fix the issue. If you
would like to create a new feature to address the issue, you are welcome to fork the
repository and create a pull request to discuss it in more
detail. These will be triaged by the CoderData team as they are received.

## Add your own data

CoderData is designed to be federated and therefore you can build your
own dataset that can be accessed locally. 

###  Documentation and steps
This process is documented
on our [GitHub site](http://github.com/pnnl-compbio/coderdata) under
'Adding a new dataset'.

### Considerations for building a dataset

Considerations to include are:
- Do you have unique sample identifiers and metadata for your samples
  needed to populate the `Sample` table of the [linkML
  schema](http://github.com/pnnl-compbio/coderdata/schema)?
- Do you have drug structure or names that can be used to carry out
  the lookup process using our [standard drug lookup
  script](http://github.com/pnnl-compbio/coderdata/build/utils/pubchem_retrieval.py)?
- Have you standardized your omics values to conform with the schema? 
- Do you have dose and response information for your experiment files
  so that the curves can be refit using the [drug curve fitting
  tool](http://github.com/pnnl-compbio/coderdata/build/utils/fit_curve.py)?
- Are your omics datasets in the correct format for our schema? 

Lastly, check out examples! We have numerous Docker files in our
[Dockerfile
directory](http://github.com/pnnl-compbio/coderdata/build/docker),
and multiple datasets in our [build
directory](http://github.com/pnnl-compbio/coderdata/build). 

---  
Your contributions are essential to the growth and improvement of CoderData. We look forward to collaborating with you!  
