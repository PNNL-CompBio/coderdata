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


## Data Source Reference List

| Dataset | Data Source | Resource | Authors | AACR Reference Number |
|---------|-------------|----------|---------|-----------------------|
| DepMap / Sanger | PharmacoGx - CCLE | The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity | Jordi Barretina et al. | 1
| DepMap / Sanger | PharmacoGx - gCSI | Reproducible pharmacogenomic profiling of cancer cell line panels | Peter M Haverty et al. | 2
| DepMap / Sanger | PharmacoGx - gCSI | A comprehensive transcriptional portrait of human cancer cell lines | Christiaan Klijn et al. | 3
| DepMap / Sanger | PharmacoGx - GDSC | Systematic identification of genomic markers of drug sensitivity in cancer cells | Mathew J Garnett et al. | 4
| DepMap / Sanger | PharmacoGx - GDSC | A Landscape of Pharmacogenomic Interactions in Cancer | Francesco Iorio et al. | 5
| DepMap / Sanger | PharmacoGx - GDSC | Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells | Wanjuan Yang et al. | 6
| DepMap / Sanger | PharmacoGx - PRISM | Discovering the anti-cancer potential of non-oncology drugs by systematic viability profiling | Steven M Corsello et al. | 7
| DepMap / Sanger | PharmacoGx - PRISM | High-throughput identification of genotype-specific cancer vulnerabilities in mixtures of barcoded tumor cell lines | Channing Yu et al. | 8
| DepMap / Sanger | PharmacoGx - CTRP | Harnessing Connectivity in a Large-Scale Small-Molecule Sensitivity Dataset | Brinton Seashore-Ludlow et al. | 9
| DepMap / Sanger | PharmacoGx - FIMM | Consistency in drug response profiling | John Patrick Mpindi et al. | 10
| DepMap / Sanger | PharmacoGx - FIMM | Individualized Systems Medicine Strategy to Tailor Treatments for Patients with Chemorefractory Acute Myeloid Leukemia | Tea Pemovska et al. | 11
| DepMap / Sanger | PharmacoGx - NCI60 | The NCI60 human tumour cell line anticancer drug screen | Robert H. Shoemaker | 12
| DepMap / Sanger | PharmacoGx - CTRPv2 | Correlating chemical sensitivity and basal gene expression reveals mechanism of action | Rees et al. | 13
| DepMap / Sanger | PharmacoGx - CTRPv2 | Harnessing Connectivity in a Large-Scale Small-Molecule Sensitivity Dataset | Seashore-Ludlow et al. | 14
| DepMap / Sanger | PharmacoGx - CTRPv2 | An Interactive Resource to Identify Cancer Genetic and Lineage Dependencies Targeted by Small Molecules | Bodycombe Basu et al. | 15
| DepMap / Sanger | COSMIC | COSMIC: a curated database of somatic variants and clinical data for cancer | Zbyslaw Sondka et al. | 16
| DepMap / Sanger | Cellosaurus | The Cellosaurus, a Cell-Line Knowledge Resource | Amos Bairoch | 17
| DepMap / Sanger | Cancer Cell Line Encyclopedia | Quantitative Proteomics of the Cancer Cell Line Encyclopedia | David P Nusinow et al | 18
| DepMap / Sanger | Cancer Cell Line Encyclopedia | The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity | Jordi Barretina | 19
| CPTAC | Clinical Proteomic Tumor Analysis Consortium | Simplified and Unified Access to Cancer Proteogenomic Data | Caleb M Lindgren et al. | 20
| HCMI | NCI Genomic Data Commons | Human Cancer Models Initiative | - | 21
| BeatAML | NCI Genomic Data Commons | Integrative analysis of drug response and clinical outcome in acute myeloid leukemia | Daniel Bottomly et al. | 22
| BeatAML | NCI Proteomic Data Commons | Mapping the proteogenomic landscape enables prediction of drug response in acute myeloid leukemia | James Pino et al. | 23
| MPNST | NF Data Portal | Chromosome 8 gain is associated with high-grade transformation in MPNST | David P Nusinow et al. | 24
