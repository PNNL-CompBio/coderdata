## Coderdata build process

All data collected for this package has been collated from
stable/reproducible sources using the scripts contained here. The
figure below shows a brief description of the process, which is
designed to be run serially, as new identifiers are generated as data
are added.

![CoderBuild process](coderDataBuild.jpg?raw=true "CoderBuild process")

## build_all.py script

This script initializes all docker containers, builds all datasets, validates them, and uploads them to figshare.

It requires the following authorization tokens to be set in the local environment depending on the use case:   
`SYNAPSE_AUTH_TOKEN`: Required for the beataml, mpnst, bladder, pancreatic, liver, sarcoma, novartis, colorectal and broad_sanger (DepMap files) datasets. Join the [CoderData team](https://www.synapse.org/#!Team:3503472) on Synapse and generate an access token.  
`FIGSHARE_TOKEN`: This token is required to upload to Figshare.  
`GITHUB_TOKEN`: This token is required to upload to GitHub.  

**Available arguments**:

- `--docker`: Initializes and builds all docker containers.
- `--samples`: Processes and builds the sample data files.
- `--omics`: Processes and builds the omics data files.
- `--drugs`: Processes and builds the drug data files.
- `--exp`: Processes and builds the experiment data files.
- `--misc`: Runs the final post-build step (e.g. splitting broad_sanger into its datasets).
- `--prev_drugs`: Comma-separated list of previously built drug files. They seed the PubChem lookup so only new drugs are queried.
- `--all`: Executes all available processes above (docker, samples, omics, drugs, exp, misc). This does not run the validate or figshare commands.
- `--validate`: Validates the generated datasets with the build data checks (see below) and the schema check scripts. This is automatically included if data upload occurs.
- `--previous-summary`: Optional. `validation_summary.json(.gz)` from the previous release, compared against by the build data checks.
- `--skip-build-checks`: Skips the build data checks during `--validate` (not recommended).
- `--figshare`: Uploads the datasets to Figshare. FIGSHARE_TOKEN must be set in local environment.
- `--high_mem`: Utilizes high memory mode for concurrent data processing. This has been successfully tested using 32 or more vCPUs. 
- `--dataset`: Specifies the datasets to process (default='broad_sanger,beataml,pancreatic,bladder,sarcoma,liver,novartis,colorectal,mpnst,cptac,hcmi').
- `--depmap-ready`: Confirms that the DepMap files on Synapse are up to date for the current DepMap release. Required for any build that includes broad_sanger.
- `--continue`: Resumes a partial build, skipping steps whose output already exists in `local/`. Progress is logged to `local/build_progress.log`.
- `--version`: Specifies the version number for the Figshare upload title (e.g., "0.1.29"). This must be a higher version than previously published versions.
- `--github-username`: GitHub username matching the GITHUB_TOKEN. Required to push the new Tag to the GitHub Repository.
- `--github-email`: GitHub email matching the GITHUB_TOKEN. Required to push the new Tag to the GitHub Repository.

**Example usage**:  
- Build all datasets and upload to Figshare and GitHub.  
Required tokens for the following command: `SYNAPSE_AUTH_TOKEN`, `FIGSHARE_TOKEN`, `GITHUB_TOKEN`.  
```bash
python coderbuild/build_all.py --all --high_mem --validate --depmap-ready --figshare --version 0.1.41 --github-username jjacobson95 --github-email jeremy.jacobson3402@gmail.com
```
  
- Build only the experiment files.  
**Note**: Preceding steps will not automatically be run. This assumes that docker images, samples, omics, and drugs were all previously built. Ensure all required tokens are set.   
```bash
python coderbuild/build_all.py --exp
```

## build_dataset.py script
This script builds a single dataset for **debugging purposes only**. It can help determine if a dataset will build correctly in isolation. Note that the sample and drug identifiers generated may not align with those from other datasets, so this script is not suitable for building production datasets.

It requires the following authorization tokens to be set in the local environment depending on the dataset:

`SYNAPSE_AUTH_TOKEN`: Required for most datasets (see the list above). Follow the directions above to gain access.

Available arguments:
- `--dataset`: Required. Name of the dataset to build. At a minimum, this will build the docker images.
- `--use_prev_dataset`: Optional. Prefix of the previous dataset for sample and drug ID continuation. The previous dataset files must be in the "local" directory.
- `--build`: Optional. Build the desired Dataset.
- `--validate`: Optional. Run the schema checker on the built files.
- `--continue`: Optional. Continues from where the build left off by skipping existing files in "local" directory.
- `--depmap-ready`: Required when building broad_sanger. Confirms the DepMap files on Synapse are up to date.

Example usage:

Build the broad_sanger dataset:
```bash
python coderbuild/build_dataset.py --build --dataset broad_sanger
```
Build the mpnst dataset continuing from broad_sanger sample and drug IDs:
```bash
python coderbuild/build_dataset.py --build --dataset mpnst --use_prev_dataset broad_sanger
```
Build run schema validation on hcmi dataset:
```bash
python coderbuild/build_dataset.py --dataset hcmi --validate
```
Build the broad_sanger dataset but skip previously built files in "local" directory:
```bash
python coderbuild/build_dataset.py --dataset broad_sanger --continue
```

## Build data checks (validate_build.py)

These checks look for problems the schema check cannot catch: duplicate results, ids that do not match between files, sample ids shared between datasets, copied responses, and large changes from the previous release. `build_all.py --validate` runs them automatically. They can also be run on their own against any folder of CoderData files, either a local build (`local/all_files_dir`) or a release downloaded from Figshare. Files may be plain or `.gz`. They need `pandas` and `numpy`; the pytest form also needs `pytest`.

As a command:
```bash
python coderbuild/utils/validate_build.py --dir local/all_files_dir
```
- `--dir`: Required. Folder of files to check.
- `--datasets`: Optional. Check only these datasets (default: all found in `--dir`).
- `--previous-summary`: Optional. Summary of an earlier build to compare counts with. Each release publishes one as `validation_summary.json.gz`.
- `--write-summary`: Optional. Save this build's counts, to use as the next baseline.

Results are printed to the terminal. The command exits with code 1 if there are errors; warnings do not change the exit code.

As tests:
```bash
pytest tests/test_build_data.py --build-dir local/all_files_dir
```
Add `--build-datasets beataml,hcmi` to check only some datasets, or `--previous-summary <file>` to compare with a previous build. Each dataset is one test. Errors fail the test, and warnings show as pytest warnings. Without `--build-dir` these tests are skipped.

In `build_all.py` the checks run after the id mapping and before the schema check, and write `validation_summary.json` into `local/all_files_dir` so it is published with the release. If they find errors, the build stops after the schema check.




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

