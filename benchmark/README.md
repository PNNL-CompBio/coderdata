### Benchmarking Models


This repository contains scripts and utilities for benchmarking deep learning models with drug response prediction. The models were developed by the [CDRP](https://github.com/PNNL-CompBio/cdrp) Repository maintainers. We borrow the following scripts from them data_utils.py, gnn_utils.py, deeptta_rna_model.py, transformer.py, and the shared_input resources.

This directory integrates key components from [CoderData](https://pypi.org/project/coderdata/) for data loading, handling, and splitting. The benchmarks are designed to compare different model architectures and data splits.


Goals:
1. **Balanced Splits**: Evaluate the role of balanced splits in self and cross-dataset performance.
2. **Cross-System Prediction**: Assess the feasibility of training on cell lines and testing on ex vivo systems.
3. **Omics Comparison**: Compare proteomics to transcriptomics as predictors of drug response.
4. **Drug Representation**: Compare various drug representation models across different datasets.



## Data Download & Environment Setup

### Step 1: Download Data

Open a terminal and execute the following commands to create the `data/` directory and download the required datasets using the `coderdata` command-line tool:

```bash
# Create the data directory and navigate into it
mkdir data
cd data

# Download a specific dataset (replace {dataset} with the desired dataset name, e.g., "ccle")
coderdata download --name {dataset}

# Alternatively, download all available datasets
coderdata download --name all
```

## Step 2: Set Up Python Environment

Activate your preferred Python environment (using Conda, virtualenv, etc.) and install the required dependencies:

```bash
# Install the required packages
pip install -r requirements.txt
```

# How to Run the Scripts

This repository includes two main scripts for benchmarking drug response prediction models using CoderData: `train.py` and `test.py`.

---

## Training with `train.py`

The `train.py` script is used for training, validating, and optionally self-testing your deep learning model. It performs the following functions:

- **Data Loading & Preprocessing**: Loads the selected dataset via CoderData, filters and processes transcriptomics and experiments data (including log transformation, scaling, and handling missing values).
- **Data Splitting**: Creates train/validation/test splits when in self-test mode, or train/validation splits when preparing for external testing.
- **Model Training**: Sets up the model (supporting various encoder types such as `transformer`, `gnn`, `morganfp`, and `descriptor`), optimizer, and loss function. Training is conducted with early stopping to avoid overfitting.
- **Checkpointing**: Saves the best model checkpoint during training.
- **Evaluation & Plotting**: In self-test mode, the script evaluates the trained model on the test set and generates a loss plot.
- **Results Logging**: Appends training and evaluation metrics to a user-specified output file.

### Command-Line Options for `train.py`

- **`--dataset`**: Dataset to use for training. Options include: `ccle`, `prism`, `beataml`, `mpnst`, `gdscv1`, `gdscv2`, `gcsi`, `ctrpv2`, `fimm`, `nci60`.  
  *Default*: `ccle`

- **`--epochs`**: Number of training epochs.  
  *Default*: `100`

- **`--encoder`**: Encoder/model type. Options: `transformer`, `gnn`, `morganfp`, `descriptor`.  
  *Default*: `transformer`

- **`--test_type`**:  
  - `self`: Performs a full train/validation/test split and evaluates the model after training.  
  - `external`: Trains and validates only, so that testing can be performed later with `test.py`.  
  *Default*: `self`

- **`--output`**: Path to the file where results will be appended.  
  *Default*: `results/results.txt`

#### Example

To train a model using the GNN encoder on the `ccle` dataset for 100 epochs while preparing for external testing, run:

```bash
python train.py --dataset ccle --epochs 100 --encoder gnn --test_type external --output external_results.txt
```

## Testing with `test.py`

The `test.py` script is used to evaluate a previously trained model on an external test dataset. It performs the following tasks:

- **Data Loading & Preprocessing**: Loads and processes the external test dataset using the same steps as in `train.py` (e.g., log transformation, scaling, and merging of drug descriptor data).
- **Model Loading**: Loads a saved model checkpoint.
- **Evaluation**: Computes evaluation metrics such as RMSE, Pearson correlation, and Spearman correlation.
- **Results Logging**: Appends the evaluation metrics to a user-specified output file.

### Command-Line Options for `test.py`

- **`--dataset`**: External test dataset to use. Options include: `ccle`, `prism`, `beataml`, `mpnst`, `gdscv1`, `gdscv2`, `gcsi`, `ctrpv2`, `fimm`, `nci60`.  
  *Required*

- **`--encoder`**: Encoder/model type used.  
  *Default*: `transformer`

- **`--ckpt`**: Path to the saved model checkpoint file to load.  
  *Required*

- **`--output`**: Path to the file where test results will be appended.  
  *Default*: `results/results.txt`

### Example

To evaluate a model trained on `ccle` with the `tranformer` encoder on the `gcsi` dataset using a saved checkpoint, run:

```bash
python test.py --dataset gcsi --encoder transformer --ckpt tmp/best_ccle_transformer_external.pt --output test_results.txt
```


## Additional Notes

- **Directory Creation**: Both `train.py` and `test.py` automatically create the necessary directories (`models`, `results`, `tmp`, `shared_input`, `plots`, and `data`) if they do not exist.

- **Data Processing**: The scripts filter and process transcriptomics and experiments data to ensure consistency (e.g., matching sample IDs, log transformation, scaling).

- **Early Stopping**: Training employs early stopping to prevent overfitting. The best model checkpoint is saved and later reloaded for evaluation.


