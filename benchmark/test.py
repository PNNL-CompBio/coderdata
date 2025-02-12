#!/usr/bin/env python
"""
test.py

This script loads a saved model checkpoint and applies it to an external test dataset.
It implements the same data formatting steps as in your original code—including processing
of transcriptomics and experiments data—and supports all four encoder/model types 
(transformer, gnn, morganfp, descriptor). Test results are appended to the specified output file.
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
        choices=["ccle", "prism", "beataml", "mpnst", "gdscv1", "gdscv2", "gcsi", "ctrpv2", "fimm","nci60"],
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
    
    # If a gene mapping file exists, filter transcriptomics data accordingly.
    # (The original code filtered by intersecting with the selected genes.)
    selected_gene_path = "./shared_input/graphDRP_landmark_genes_map.txt"
    if os.path.exists(selected_gene_path):
        selected_gene_df = pd.read_csv(selected_gene_path, sep='\t')
        selected_genes = set(selected_gene_df['To'])
        cd_test.transcriptomics = cd_test.transcriptomics[cd_test.transcriptomics['entrez_id'].isin(selected_genes)]

    # Ensure omics and experiments have the same samples
    common_ids = set(cd_test.experiments['improve_sample_id'].unique()).intersection(set(cd_test.transcriptomics['improve_sample_id'].unique()))
    cd_test.experiments = cd_test.experiments[cd_test.experiments['improve_sample_id'].isin(common_ids)].reset_index(drop=True)
    cd_test.transcriptomics = cd_test.transcriptomics[cd_test.transcriptomics['improve_sample_id'].isin(common_ids)]

    # Final Data Pre-Processing - this should be fixed on the data build side eventually.
    cd_test.experiments.dropna(inplace = True)
    cd_test.experiments = cd_test.experiments[cd_test.experiments.dose_response_metric == "fit_auc"]

    # --------------------------
    # Process Transcriptomics Data
    # --------------------------
    print(f"{cd_test} in progress")
    cd_test.transcriptomics = cd_test.format(data_type='transcriptomics', shape='wide', inplace=True)
    # Log transform the transcriptomics data
    cd_test.transcriptomics = np.log1p(cd_test.transcriptomics)
    # Scale the data
    scaler = StandardScaler()
    train_gene_exp_scaled = scaler.fit_transform(cd_test.transcriptomics)
    train_gene_exp_scaled = pd.DataFrame(train_gene_exp_scaled, 
                                         index=cd_test.transcriptomics.index, 
                                         columns=cd_test.transcriptomics.columns)
    # Remove columns with NaN values
    cd_test.transcriptomics = train_gene_exp_scaled.dropna(axis=1)
    
    
    # --------------------------
    # Process Experiments Data
    # --------------------------
    cd_test.experiments = cd_test.format(data_type='experiments', shape='wide', metrics=['fit_auc'])
    cd_test.experiments = pd.merge(cd_test.experiments, cd_test.drugs, on="improve_drug_id", how='left')
    cd_test.experiments = cd_test.experiments[["improve_sample_id", "improve_drug_id", "fit_auc", "canSMILES"]]
    cd_test.experiments.drop_duplicates(ignore_index=True, inplace=True)
    
    # For certain encoder types, rename the "canSMILES" column to "smiles"
    if args.encoder in ["gnn", "morganfp", "descriptor"]:
        cd_test.experiments.rename(columns={"canSMILES": "smiles"}, inplace=True)

    n_descriptors = None
    feature_names = None
    if args.encoder in ["descriptor"]:
        # Process drug descriptors for descriptor encoding.
        cd_test.drug_descriptors = pd.merge(cd_test.drug_descriptors, cd_test.drugs, on="improve_drug_id", how='left')
        cd_test.drug_descriptors = cd_test.drug_descriptors[["improve_drug_id", "canSMILES", "structural_descriptor", "descriptor_value"]]
        cd_test.drug_descriptors.rename(columns={"canSMILES": "smiles"}, inplace=True)
        cd_test.drug_descriptors = cd_test.drug_descriptors[cd_test.drug_descriptors["structural_descriptor"].str.startswith("n")]
        cd_test.drug_descriptors = cd_test.drug_descriptors[~cd_test.drug_descriptors["structural_descriptor"].str.contains("Ring")]
        cd_test.drug_descriptors.drop(columns="improve_drug_id", inplace=True)
        
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
    # Create DataLoader for test data using CreateData function
    # --------------------------
    data_creater = CreateData(gexp=cd_test.transcriptomics,
                              encoder_type=args.encoder,
                              metric="fit_auc",
                              data_path="shared_input/",
                              feature_names=feature_names)
    print("Creating test_ds...")
    test_ds = data_creater.create_data(cd_test.experiments)
    test_loader = DataLoader(test_ds, batch_size=64, shuffle=False, drop_last=False)
    
    # --------------------------
    # Set up and load the model
    # --------------------------
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    n_genes = len(cd_test.transcriptomics.columns)
    if args.encoder == "gnn":
        model = Model(gnn_features=65, n_descriptors=n_descriptors, encoder_type=args.encoder, n_genes=n_genes).to(device)
    else:
        model = Model(gnn_features=None, n_descriptors=n_descriptors, encoder_type=args.encoder, n_genes=n_genes).to(device)
    
    model.load_state_dict(torch.load(args.ckpt, map_location=device))
    print(f"Loaded checkpoint from {args.ckpt}")
    
    # --------------------------
    # Evaluate the model on the external test set
    # --------------------------
    test_rmse, pearson_corr, spearman_corr, _, _ = test_fn(test_loader, model, device)
    print(f"\nExternal Test Results -- RMSE: {test_rmse:.3f}, Pearson: {pearson_corr:.3f}, Spearman: {spearman_corr:.3f}")
    
    # --------------------------
    # Append test results to output file
    # --------------------------
    results_str = (
        f"External Test Results for model with encoder {args.encoder} on dataset {args.dataset}:\n"
        "Test RMSE\tPearson Correlation\tSpearman Correlation\n"
        f"{test_rmse:.3f}\t{pearson_corr:.3f}\t{spearman_corr:.3f}\n\n"
    )
    with open(args.output, 'a') as f:
        f.write(results_str)
    print(f"Results appended to {args.output}")

if __name__ == '__main__':
    main()
