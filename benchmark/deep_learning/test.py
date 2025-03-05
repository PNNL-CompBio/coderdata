#!/usr/bin/env python
"""
test.py

This script loads a saved model checkpoint and applies it to an external test dataset.
It implements the same data formatting steps as in train.py—including processing
of transcriptomics and experiments data—and supports all four encoder/model types 
(transformer, gnn, morganfp, descriptor). Test results are appended to the specified output file.

New arguments --gene_selection and --gene_number have been added so that the test data
are filtered using the same gene selection method as in training.
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
        "--gene_selection",
        type=str,
        choices=["landmark", "pca", "mutation", "variance"],
        default="landmark",
        help=("Gene selection method: "
              "'landmark' = use graphDRP_landmark_genes_map file (default), "
              "'pca' = PCA-based, "
              "'mutation' = Mutation counts-based, "
              "'variance' = Transcriptomic variance-based")
    )
    parser.add_argument(
        "--gene_number",
        type=int,
        default=500,
        help="Number of top genes to retain for ranking methods (default: 500)"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="results/results.txt",
        help="Path to results output file (results will be appended; default: results/results.txt)"
    )
    args = parser.parse_args()
    gene_number = args.gene_number

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
    
    # Apply gene selection to transcriptomics data based on the chosen method.
    if args.gene_selection == "landmark":
        selected_gene_path = "./shared_input/graphDRP_landmark_genes_map.txt"
        if os.path.exists(selected_gene_path):
            try:
                df_landmark = pd.read_csv(selected_gene_path, sep='\t')
                if 'To' not in df_landmark.columns:
                    raise ValueError("Column 'To' not found in landmark gene file.")
                selected_genes = set(df_landmark['To'].dropna().unique())
                cd_test.transcriptomics = cd_test.transcriptomics[cd_test.transcriptomics['entrez_id'].isin(selected_genes)]
                print(f"Landmark gene selection applied: {len(selected_genes)} genes.")
            except Exception as e:
                print(f"Error reading landmark gene file: {e}. Proceeding without gene filtering.")
        else:
            print(f"Warning: {selected_gene_path} not found. Proceeding without gene filtering.")
    
    elif args.gene_selection == "pca":
        transcriptomics_df = cd_test.transcriptomics.copy()
        expr_cols = [col for col in transcriptomics_df.columns if col not in ['improve_sample_id', 'entrez_id']]
        if not expr_cols:
            print("Warning: No expression column found for PCA gene selection. Skipping gene filtering.")
        else:
            try:
                transcriptomics_wide = cd.format(cd_test, "transcriptomics", "wide")
            except Exception as e:
                print(f"Error formatting transcriptomics data for PCA: {e}. Skipping PCA gene selection.")
                transcriptomics_wide = pd.DataFrame()
            if transcriptomics_wide.empty:
                print("Warning: Formatted transcriptomics data is empty. Skipping PCA gene selection.")
            else:
                transcriptomics_wide_filled = transcriptomics_wide.fillna(0)
                scaler = StandardScaler()
                scaled_data = scaler.fit_transform(transcriptomics_wide_filled)
                from sklearn.decomposition import PCA
                pca = PCA()
                pca.fit(scaled_data)
                loadings = pd.DataFrame(
                    pca.components_.T,
                    index=transcriptomics_wide_filled.columns,
                    columns=[f'PC{i+1}' for i in range(pca.n_components_)]
                )
                cumulative_variance = np.cumsum(pca.explained_variance_ratio_)
                # Use a threshold similar to training (e.g., 95%)
                num_pcs = np.where(cumulative_variance >= 0.95)[0][0] + 1
                top_n_genes_per_pc = 10
                selected_pcs = [f'PC{i+1}' for i in range(num_pcs)]
                top_genes_set = set()
                for pc in selected_pcs:
                    abs_loadings = loadings[pc].abs()
                    top_genes_pc = abs_loadings.sort_values(ascending=False).head(top_n_genes_per_pc).index.tolist()
                    top_genes_set.update(top_genes_pc)
                if len(top_genes_set) > gene_number:
                    avg_loadings = loadings.loc[list(top_genes_set), selected_pcs].abs().mean(axis=1)
                    pca_top_genes = avg_loadings.sort_values(ascending=False).head(gene_number).index.tolist()
                else:
                    pca_top_genes = list(top_genes_set)
                print(f"Top {gene_number} genes based on PCA loadings (selected {len(pca_top_genes)} genes):")
                print(pca_top_genes)
                cd_test.transcriptomics = cd_test.transcriptomics[cd_test.transcriptomics['entrez_id'].isin(pca_top_genes)]
    
    elif args.gene_selection == "mutation":
        if hasattr(cd_test, 'mutations'):
            mutations_df = cd_test.mutations.copy()
            non_silent = mutations_df[mutations_df['variant_classification'] != 'Silent']
            mutation_counts = non_silent['entrez_id'].value_counts()
            mutation_counts_df = mutation_counts.reset_index()
            mutation_counts_df.columns = ['entrez_id', 'mutation_count']
            mutation_counts_df_sorted = mutation_counts_df.sort_values(by='mutation_count', ascending=False)
            mutation_top_genes = mutation_counts_df_sorted.head(gene_number)['entrez_id'].tolist()
            print(f"Top {gene_number} genes with highest mutation frequency:")
            print(mutation_top_genes)
            cd_test.transcriptomics = cd_test.transcriptomics[cd_test.transcriptomics['entrez_id'].isin(mutation_top_genes)]
        else:
            print("Warning: Mutations data not available. Proceeding without mutation-based gene selection.")
    
    elif args.gene_selection == "variance":
        transcriptomics_df = cd_test.transcriptomics.copy()
        expr_cols = [col for col in transcriptomics_df.columns if col not in ['improve_sample_id', 'entrez_id']]
        if not expr_cols:
            print("Warning: No expression column found for variance gene selection. Skipping gene filtering.")
        else:
            try:
                transcriptomics_wide = cd.format(cd_test, "transcriptomics", "wide")
            except Exception as e:
                print(f"Error formatting transcriptomics data for variance selection: {e}. Skipping variance-based gene selection.")
                transcriptomics_wide = pd.DataFrame()
            if transcriptomics_wide.empty:
                print("Warning: Formatted transcriptomics data is empty. Skipping variance-based gene selection.")
            else:
                gene_variances = transcriptomics_wide.var(axis=0)
                top_genes_by_variance = gene_variances.sort_values(ascending=False).head(gene_number).index.tolist()
                print(f"Top {gene_number} genes with highest expression variability:")
                print(top_genes_by_variance)
                cd_test.transcriptomics = cd_test.transcriptomics[cd_test.transcriptomics['entrez_id'].isin(top_genes_by_variance)]
    
    # --------------------------
    # Ensure omics and experiments have the same samples
    # --------------------------
    common_ids = set(cd_test.experiments['improve_sample_id'].unique()).intersection(
        set(cd_test.transcriptomics['improve_sample_id'].unique())
    )
    cd_test.experiments = cd_test.experiments[cd_test.experiments['improve_sample_id'].isin(common_ids)].reset_index(drop=True)
    cd_test.transcriptomics = cd_test.transcriptomics[cd_test.transcriptomics['improve_sample_id'].isin(common_ids)]
    
    # Final Data Pre-Processing - this should eventually be handled upstream.
    cd_test.experiments.dropna(inplace=True)
    cd_test.experiments = cd_test.experiments[cd_test.experiments.dose_response_metric == "fit_auc"]

    # --------------------------
    # Process Transcriptomics Data
    # --------------------------
    print(f"Formatting transcriptomics data for {cd_test}...")
    cd_test.transcriptomics = cd_test.format(data_type='transcriptomics', shape='wide', inplace=True)
    cd_test.transcriptomics = np.log1p(cd_test.transcriptomics)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(cd_test.transcriptomics)
    scaled_df = pd.DataFrame(scaled_data, index=cd_test.transcriptomics.index, columns=cd_test.transcriptomics.columns)
    cd_test.transcriptomics = scaled_df.dropna(axis=1)
    
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
    # Create DataLoader for test data
    # --------------------------
    data_creater = CreateData(gexp=cd_test.transcriptomics,
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
        f"External Test Results for model with encoder {args.encoder} on dataset {args.dataset}\n"
        f"Gene Selection: {args.gene_selection}, Top Genes: {gene_number}\n"
        "Test RMSE\tPearson Correlation\tSpearman Correlation\n"
        f"{test_rmse:.3f}\t{pearson_corr:.3f}\t{spearman_corr:.3f}\n\n"
    )
    with open(args.output, 'a') as f:
        f.write(results_str)
    print(f"Results appended to {args.output}")

if __name__ == '__main__':
    main()