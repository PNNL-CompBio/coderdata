#!/usr/bin/env python
"""
test.py

This script loads a saved model checkpoint and applies it to an external test dataset.
It implements the same data formatting steps as in train.py—including processing
of transcriptomics/proteomics and experiments data—and supports all four encoder/model types 
(transformer, gnn, morganfp, descriptor). Instead of performing its own gene selection,
this version reads in the final gene list produced during training to ensure consistent
gene ordering and prevent size mismatches.
"""

import os
import argparse
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from torch_geometric.loader import DataLoader
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
import coderdata as cd
from sklearn.decomposition import PCA

# Local project imports
from data_utils import DataProcessor, add_smiles, average_dose_response_value, filter_exp_data
from gnn_utils import CreateData, EarlyStopping, test_fn
from deeptta_rna_model import Model

def main():
    # --------------------------
    # Parse command-line arguments
    # --------------------------
    parser = argparse.ArgumentParser(
        description="Test deep learning model on an external dataset"
    )
    parser.add_argument(
        "--dataset",
        type=str,
        choices=["ccle", "prism", "beataml", "mpnst", "gdscv1", "gdscv2", "gcsi", "ctrpv2", "fimm", "nci60"],
        required=True,
        help="External test dataset to use"
    )
    parser.add_argument(
        "--encoder",
        type=str,
        choices=["transformer", "gnn", "morganfp", "descriptor"],
        default="transformer",
        help="Encoder type (default: transformer)"
    )
    parser.add_argument(
        "--ckpt",
        type=str,
        required=True,
        help="Path to the model checkpoint file to load"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="results/results.txt",
        help="Path to results output file (results will be appended; default: results/results.txt)"
    )
    parser.add_argument(
        "--gene_list",
        type=str,
        required=True,
        help="Path to gene list file."
    )
    # New argument for selecting the omics data type (transcriptomics or proteomics)
    parser.add_argument(
        "--omics",
        type=str,
        choices=["transcriptomics", "proteomics"],
        default="transcriptomics",
        help="Omics data type to use (default: transcriptomics)"
    )
    args = parser.parse_args()

    # --------------------------
    # Create necessary directories
    # --------------------------
    for folder in ["models", "results", "tmp", "shared_input", "plots", "data"]:
        if not os.path.exists(folder):
            os.makedirs(folder)
            print(f"'{folder}' directory created.")
    
    # --------------------------
    # Load external test dataset using coderdata
    # --------------------------
    print(f"Loading external test dataset: {args.dataset}")
    cd_test = cd.load(args.dataset, "data")
    
    # --------------------------
    # Ensure omics and experiments have the same samples
    # --------------------------
    # Use the selected omics type instead of hardcoding transcriptomics
    omics_field = args.omics
    omics_data = getattr(cd_test, omics_field)
    common_ids = set(cd_test.experiments['improve_sample_id'].unique()).intersection(
        set(omics_data['improve_sample_id'].unique())
    )
    cd_test.experiments = cd_test.experiments[cd_test.experiments['improve_sample_id'].isin(common_ids)].reset_index(drop=True)
    omics_data = omics_data[omics_data['improve_sample_id'].isin(common_ids)]
    setattr(cd_test, omics_field, omics_data)
    
    # Final Data Pre-Processing - this should eventually be handled upstream.
    cd_test.experiments.dropna(inplace=True)
    cd_test.experiments = cd_test.experiments[cd_test.experiments.dose_response_metric == "fit_auc"]

    # --------------------------
    # Process Omics Data
    # --------------------------
    print(f"Formatting {omics_field} data for {cd_test}...")
    omics_formatted = cd_test.format(data_type=omics_field, shape='wide', inplace=True)
    # Apply log1p transformation only for transcriptomics
    if args.omics == "transcriptomics":
        omics_formatted = np.log1p(omics_formatted)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(omics_formatted)
    scaled_df = pd.DataFrame(scaled_data, index=omics_formatted.index, columns=omics_formatted.columns)
    omics_formatted = scaled_df.dropna(axis=1)
    setattr(cd_test, omics_field, omics_formatted)
    
    # --------------------------
    # Instead of performing gene selection, load the final gene list produced during training.
    # This ensures the same genes (in the same order) are used during testing.
    # --------------------------
    gene_list_file = args.gene_list
    if os.path.exists(gene_list_file):
        with open(gene_list_file, 'r') as f:
            final_gene_list = [line.strip() for line in f if line.strip()]
        print(f"Loaded gene list from {gene_list_file} with {len(final_gene_list)} genes.")
        current_genes = getattr(cd_test, omics_field).columns.tolist()
        print(f"final_gene_list: {final_gene_list}")
        print(f"current_genes: {current_genes}")
        # Assuming the gene identifiers are integers as in the training script
        ordered_genes = [int(gene) for gene in final_gene_list if int(gene) in current_genes]
        print(f"ordered_genes: {ordered_genes}")
        setattr(cd_test, omics_field, getattr(cd_test, omics_field)[ordered_genes])
        print(f"{omics_field} data after gene filtering: {getattr(cd_test, omics_field)}")
    else:
        print("Gene list file not found, proceeding without gene filtering.")
    
    # --------------------------
    # Process Experiments Data
    # --------------------------
    cd_test.experiments = cd_test.format(data_type='experiments', shape='wide', metrics=['fit_auc'])
    cd_test.experiments = pd.merge(cd_test.experiments, cd_test.drugs, on="improve_drug_id", how='left')
    cd_test.experiments = cd_test.experiments[["improve_sample_id", "improve_drug_id", "fit_auc", "canSMILES"]]
    cd_test.experiments.drop_duplicates(ignore_index=True, inplace=True)
    if args.encoder in ["gnn", "morganfp", "descriptor"]:
        cd_test.experiments.rename(columns={"canSMILES": "smiles"}, inplace=True)
    
    n_descriptors = None
    feature_names = None
    if args.encoder == "descriptor":
        cd_test.drug_descriptors = pd.merge(cd_test.drug_descriptors, cd_test.drugs, on="improve_drug_id", how='left')
        cd_test.drug_descriptors = cd_test.drug_descriptors[["improve_drug_id", "canSMILES", "structural_descriptor", "descriptor_value"]]
        cd_test.drug_descriptors.rename(columns={"canSMILES": "smiles"}, inplace=True)
        cd_test.drug_descriptors = cd_test.drug_descriptors[cd_test.drug_descriptors["structural_descriptor"].str.startswith("n")]
        cd_test.drug_descriptors = cd_test.drug_descriptors[~cd_test.drug_descriptors["structural_descriptor"].str.contains("Ring")]
        cd_test.drug_descriptors.drop(columns="improve_sample_id", inplace=True)
        cd_test.drug_descriptors = cd_test.drug_descriptors.pivot_table(
            index=["smiles"],
            columns="structural_descriptor",
            values="descriptor_value",
            aggfunc='first'
        ).reset_index()
        features = cd_test.drug_descriptors
        n_descriptors = features.shape[1] - 1
        feature_names = features.drop(['smiles'], axis=1).columns.tolist()
        cd_test.experiments = pd.merge(cd_test.experiments, features, on='smiles', how='left')
        n_columns = [col for col in cd_test.experiments.columns if col.startswith("n")]
        cd_test.experiments[n_columns] = cd_test.experiments[n_columns].astype(float)
    
    # --------------------------
    # Create DataLoader for test data using the selected omics data
    # --------------------------
    data_creater = CreateData(gexp=getattr(cd_test, omics_field),
                              encoder_type=args.encoder,
                              metric="fit_auc",
                              data_path="shared_input/",
                              feature_names=feature_names)
    print("Creating test dataset...")
    test_ds = data_creater.create_data(cd_test.experiments)
    test_loader = DataLoader(test_ds, batch_size=64, shuffle=False, drop_last=False)
    
    # --------------------------
    # Set up and load the model
    # --------------------------
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    n_genes = len(getattr(cd_test, omics_field).columns)
    if args.encoder == "gnn":
        model = Model(gnn_features=65, n_descriptors=n_descriptors, encoder_type=args.encoder, n_genes=n_genes).to(device)
    else:
        model = Model(gnn_features=None, n_descriptors=n_descriptors, encoder_type=args.encoder, n_genes=n_genes).to(device)
    
    model.load_state_dict(torch.load(args.ckpt, map_location=device))
    print(f"Loaded checkpoint from {args.ckpt}")
    
    # --------------------------
    # Evaluate the model on the external test set
    # --------------------------
    mse, mae, r2, pearson_corr, spearman_corr, target_list, predicted_list = test_fn(test_loader, model, device)
    print(f"\nExternal Test Results -- RMSE: {mse:.3f}, Pearson: {pearson_corr:.3f}, Spearman: {spearman_corr:.3f}")
    
    # --------------------------
    # Append test results to output file
    # --------------------------
    results_str = (
        f"External Test Results for model with encoder {args.encoder} on dataset {args.dataset}\n"
        f"Gene list used from file: {gene_list_file}\n"
        "Test RMSE\tPearson Correlation\tSpearman Correlation\n"
        f"{mse:.3f}\t{pearson_corr:.3f}\t{spearman_corr:.3f}\n\n"
    )
    with open(args.output, 'a') as f:
        f.write(results_str)
    print(f"Results appended to {args.output}")

if __name__ == '__main__':
    main()
