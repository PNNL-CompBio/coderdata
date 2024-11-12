## Coderdata build process

All data collected for this package has been collated from
stable/reproducible sources using the scripts contained here. The
figure below shows a brief description of the process, which is
designed to be run serially, as new identifiers are generated as data
are added.

![Build process](coderDataBuild.jpg?raw=true "Build process")

## Build a local version using the `build_all.py` script

This script initializes all docker containers, builds all datasets, validates them, and uploads them to figshare and pypi.

It requires the following authorization tokens to be set in the local environment depending on the use case:   
`SYNAPSE_AUTH_TOKEN`: Required for beataml and mpnst datasets. Join the [CoderData team](https://www.synapse.org/#!Team:3503472) on Synapse and generate an access token.
`PYPI_TOKEN`: This token is required to upload to PyPI.
`FIGSHARE_TOKEN`: This token is required to upload to Figshare.
`GITHUB_TOKEN`: This token is required to upload to GitHub.
**Available arguments**:

- `--docker`: Initializes and builds all docker containers.
- `--samples`: Processes and builds the sample data files.
- `--omics`: Processes and builds the omics data files.
- `--drugs`: Processes and builds the drug data files.
- `--exp`: Processes and builds the experiment data files.
- `--all`: Executes all available processes above (docker, samples, omics, drugs, exp). This does not run the validate, figshare, or pypi commands.
- `--validate`: Validates the generated datasets using the schema check scripts. This is automatically included if data upload occurs.
- `--figshare`: Uploads the datasets to Figshare. FIGSHARE_TOKEN must be set in local environment.
- `--pypi`: Uploads the package to PyPI. PYPI_TOKEN must be set in local environment.
- `--high_mem`: Utilizes high memory mode for concurrent data processing. This has been successfully tested using 32 or more vCPUs. 
- `--dataset`: Specifies the datasets to process (default='broad_sanger,hcmi,beataml,mpnst,cptac').
- `--version`: Specifies the version number for the PyPI package and Figshare upload title (e.g., "0.1.29"). This is required for figshare and PyPI upload steps. This must be a higher version than previously published versions.
- `--github-username`: GitHub username matching the GITHUB_TOKEN. Required to push the new Tag to the GitHub Repository.
- `--github-email`: GitHub email matching the GITHUB_TOKEN. Required to push the new Tag to the GitHub Repository.

**Example usage**:
- Build all datasets and upload to Figshare and PyPI and GitHub.
Required tokens for the following command: `SYNAPSE_AUTH_TOKEN`, `PYPI_TOKEN`, `FIGSHARE_TOKEN`, `GITHUB_TOKEN`.
```bash
python build/build_all.py --all --high_mem --validate --pypi --figshare --version 0.1.41 --github-username jjacobson95 --github-email jeremy.jacobson3402@gmail.com
```
  
- Build only the experiment files.
**Note**: Preceding steps will not automatically be run. This assumes that docker images, samples, omics, and drugs were all previously built. Ensure all required tokens are set.
```bash
python build/build_all.py --exp
```

## Build/test individual datset using the `build_dataset.py` script
This script builds a single dataset for **debugging purposes only**. It can help determine if a dataset will build correctly in isolation. Note that the sample and drug identifiers generated may not align with those from other datasets, so this script is not suitable for building production datasets.

It requires the following authorization tokens to be set in the local environment depending on the dataset:

`SYNAPSE_AUTH_TOKEN`: Required for beataml and mpnst datasets. Follow the directions above to use gain access.

Available arguments:
- `--dataset`: Required. Name of the dataset to build.
- `--use_prev_dataset`: Optional. Prefix of the previous dataset for sample and drug ID continuation. The previous dataset files must be in the "local" directory.
- `--validate`: Optional. Runs the schema checker on the built files.
- `--continue`: Optional. Continues from where the build left off by skipping existing files in "local" directory.
Example usage:

Build the broad_sanger dataset:
```bash
python build/build_dataset.py --dataset broad_sanger
```
Build the mpnst dataset continuing from broad_sanger sample and drug IDs:
```bash
python build/build_dataset.py --dataset mpnst --use_prev_dataset broad_sanger
```
Build the hcmi dataset and run validation:
```bash
python build/build_dataset.py --dataset hcmi --validate
```
Build the broad_sanger dataset but skip previously built files in "local" directory:
```bash
python build/build_dataset.py --dataset broad_sanger --continue
```

## Adding a new dataset

We have standardized the build process so an additional dataset can be
built locally or as part of the next version of coder. Here are the
steps to follow:

1. First visit the [build
directory](build/README.md) and ensure you can build a local copy of
CoderData. 

2. Checkout this repository and  create a subdirectory of the
[build directory](build) with your own build files. 

3. Develop your scripts to build the data files according to our
[LinkML Schema](schema/coderdata.yaml]). This will require collecting
the following metadata:
- entrez gene identifiers (or you can use the `genes.csv` file
- sample information such as species and model system type
- drug name that can be searched on PubChem

You can validate each file by
using the [linkML
validator](https://linkml.io/linkml/data/validating-data) together
with our schema file. 

You can use the following scripts as part of your build process:
- [build/utils/fit_curve.py](build/utils/fit_curve.py): This script
  takes dose-response data and generates the dose-response statistics
  required by CoderData/
- [build/utils/pubchem_retrieval.py](build/utils/pubchem_retreival.py):
  This script retreives structure and drug synonym information
  required to populate the `Drug` table. 

4. Wrap your scripts in standard shell scripts with the following names
and arguments:

| shell script     | arguments                | description         |
|------------------|--------------------------|---------------------|
| `build_samples.sh` | [latest_samples] | Latest version of samples generated by coderdata build |
| `build_omics.sh` | [gene file] [samplefile] | This includes the `genes.csv` that was generated in the original build as well as the sample file generated above. |
| `build_drugs.sh` | [drugfile1,drugfile2,...]       | This includes a comma-delimited list of all drugs files generated from previous build  |
| `build_exp.sh`| [samplfile ] [drugfile] | sample file and drug file generated by previous scripts |

5. Put the Docker container file inside the [Docker
directory](./build/docker) with the name
`Dockerfile.[datasetname]`. 

6. Run `build_all.py` from the root directory, which should now add in
your Dockerfile in the mix and call the scripts in your Docker
container to build the files.






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

