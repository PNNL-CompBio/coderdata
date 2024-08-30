## Coderdata build process

All data collected for this package has been collated from
stable/reproducible sources using the scripts contained here. The
figure below shows a brief description of the process, which is
designed to be run serially, as new identifiers are generated as data
are added.

![Build process](coderDataBuild.jpg?raw=true "Build process")

## build_all.py script

This script initializes all docker containers, builds all datasets, validates them, and uploads them to figshare and pypi.

It requires the following authorization tokens to be set in the local environment depending on the use case:  
`SYNAPSE_AUTH_TOKEN`: Required for beataml and mpnst datasets. Join the [CoderData team](https://www.synapse.org/#!Team:3503472) on Synapse and generate an access token.
`PYPI_TOKEN`: This token is required to upload to PyPI.
`FIGSHARE_TOKEN`: This token is required to upload to Figshare.

Available arguments:

- `--docker`: Initializes and builds all docker containers.
- `--samples`: Processes and builds the sample data files.
- `--omics`: Processes and builds the omics data files.
- `--drugs`: Processes and builds the drug data files.
- `--exp`: Processes and builds the experiment data files.
- `--all`: Executes all available processes above (docker, samples, omics, drugs, exp).
- `--validate`: Validates the generated datasets using the schema check scripts.
- `--figshare`: Uploads the datasets to Figshare.
- `--pypi`: Uploads the package to PyPI.
- `--high_mem`: Utilizes high memory mode for concurrent data processing. This has been successfully tested using 32 or more vCPUs.
- `--dataset`: Specifies the datasets to process (default='broad_sanger,hcmi,beataml,mpnst,cptac').
- `--version`: Specifies the version number for the package and data upload title. This is required to upload to figshare and PyPI

Example usage:
Build all datasets and upload to Figshare and PyPI. (Ensure tokens are set locally first)
```bash
python build/build_all.py --all --high_mem --validate --pypi --figshare --version 0.1.29
```

To only build the experiment files use the following. **Note**: Preceding steps will not automatically be run. This assumes that docker images, samples, omics, and drugs were all previously built. 
```bash
python build/build_all.py --exp
```

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


