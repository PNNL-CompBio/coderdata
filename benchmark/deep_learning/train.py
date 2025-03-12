#!/usr/bin/env python
"""
train.py

This script trains, validates, and (in self-test mode) tests a deep learning model
for drug response prediction. It implements all of the original data formatting and
model-type logic (for transformer, gnn, morganfp, and descriptor encoders).
Use the command-line options to select the dataset, number of epochs,
encoder type, whether to perform self-testing or prepare for external testing, and
an output file (to which results are appended).

A new argument, --gene_selection, allows choosing among alternative gene-selection methods:
  "landmark"  : Use graphDRP_landmark_genes_map file (default)
  "pca"       : PCA-based gene selection
  "mutation"  : Mutation counts-based selection
  "variance"  : Omics variance-based selection
  "random"    : A random selection - this uses seed for selection.

A new argument, --gene_number, specifies the top X genes to retain for ranking methods (default: 1000).

A new argument, --omics allows switching the omics data type (default: transcriptomics).
"""

import os
import argparse
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from torch_geometric.loader import DataLoader
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
import coderdata as cd
import random

# Local project imports
from data_utils import DataProcessor
from gnn_utils import CreateData, EarlyStopping, test_fn
from deeptta_rna_model import Model

def main():
    # --------------------------
    # Parse command-line arguments
    # --------------------------
    parser = argparse.ArgumentParser(
        description="Train deep learning model for drug response prediction"
    )
    parser.add_argument(
        "--dataset",
        type=str,
        choices=["ccle", "prism", "beataml", "mpnst", "gdscv1", "gdscv2", "gcsi", "ctrpv2", "fimm", "nci60"],
        default="ccle",
        help="Dataset to use for training (default: ccle)"
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=100,
        help="Number of training epochs (default: 100)"
    )
    parser.add_argument(
        "--encoder",
        type=str,
        choices=["transformer", "gnn", "morganfp", "descriptor"],
        default="transformer",
        help="Encoder type to use (default: transformer)"
    )
    parser.add_argument(
        "--test_type",
        type=str,
        choices=["self", "external"],
        default="self",
        help="If 'self', the script performs train/val/test splitting on the chosen dataset; if 'external', only training (and validation) is done so that testing may be done later via test.py (default: self)"
    )
    parser.add_argument(
        "--split_type",
        type=str,
        choices=["mixed-set", "cancer-blind", "drug-blind"],
        default="mixed-set",
        help="Split type options."
    )
    parser.add_argument(
        "--output",
        type=str,
        default="results/results.txt",
        help="Path to results output file (results will be appended; default: results/results.txt)"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Select a random seed to use for the split function."
    )
    parser.add_argument(
        "--gene_selection",
        type=str,
        choices=["landmark", "pca", "mutation", "variance", "random"],
        default="landmark",
        help=("Gene selection method: "
              "'landmark' = graphDRP_landmark_genes_map file (default), "
              "'pca' = PCA-based, "
              "'mutation' = Mutation counts-based, "
              "'variance' = Omics variance-based, "
              "'random' = random selection")
    )
    parser.add_argument(
        "--gene_number",
        type=int,
        default=1000,
        help="Number of top genes to retain for ranking methods (default: 1000)"
    )
    parser.add_argument(
        "--omics",
        type=str,
        choices=["transcriptomics", "proteomics"],
        default="transcriptomics",
        help="Omics data type to use (default: transcriptomics)"
    )
    
    args = parser.parse_args()
    gene_number = args.gene_number

    # --------------------------
    # Hyperparameters and directories
    # --------------------------
    data_split_seed = args.seed
    split_method = args.split_type
    bs = 32
    lr = 1e-4
    n_epochs = args.epochs
    encoder = args.encoder
    dose_response_metric = "fit_auc"
    
    # Set output prefix and checkpoint path based on test_type.
    output_prefix = f"{args.dataset}_{split_method}_{encoder}_{args.test_type}_{args.gene_selection}_{gene_number}_{args.omics}_{n_epochs}_{data_split_seed}"
    ckpt_path = f"./models/best_{output_prefix}.pt"
    
    # Create necessary directories
    for folder in ["models", "results", "shared_input", "plots"]:
        if not os.path.exists(folder):
            os.makedirs(folder)
            print(f"'{folder}' directory created.")

    # --------------------------
    # Load dataset using coderdata and extract the chosen omics data
    # --------------------------
    print(f"Loading dataset: {args.dataset}")
    cd_data = cd.load(args.dataset, "data")
    
    # Get the specified omics data (default is transcriptomics)
    omics_data = getattr(cd_data, args.omics, None)
    if omics_data is None:
        print(f"Warning: {args.omics} data not found in the dataset. Exiting.")
        exit(1)
    
    # --------------------------
    # Perform gene selection on the full omics data (but delay final gene-number filtering)
    # Compute candidate gene set and ranking (if applicable)
    # --------------------------
    candidate_genes = set()
    candidate_ranking = None  # Will hold an ordered list if available

    if args.gene_selection == "landmark":
        selected_gene_path = "./shared_input/graphDRP_landmark_genes_map.txt"
        if os.path.exists(selected_gene_path):
            try:
                selected_gene_df = pd.read_csv(selected_gene_path, sep='\t')
                if 'To' not in selected_gene_df.columns:
                    raise ValueError("Column 'To' not found in landmark gene file.")
                candidate_genes = set(selected_gene_df['To'])
                # For landmark selection, use alphabetical order for ranking
                candidate_ranking = sorted(candidate_genes)
                omics_data = omics_data[omics_data['entrez_id'].isin(candidate_genes)]
            except Exception as e:
                print(f"Error reading landmark gene file: {e}. Proceeding without gene filtering.")
        else:
            print(f"Warning: {selected_gene_path} not found. Proceeding without gene filtering.")

    elif args.gene_selection == "random":
        random.seed(args.seed)
        # For random selection, take all available genes as candidates (final sampling later)
        candidate_genes = set(omics_data['entrez_id'].unique())
        # Establish a consistent order (e.g., sorted) for reproducibility
        candidate_ranking = sorted(candidate_genes)

    elif args.gene_selection == "pca":
        try:
            omics_wide = cd.format(cd_data, args.omics, "wide")
        except Exception as e:
            print(f"Error formatting {args.omics} data: {e}. Skipping PCA gene selection.")
            omics_wide = pd.DataFrame()
        if omics_wide.empty:
            print("Warning: Formatted omics data is empty. Skipping PCA gene selection.")
            candidate_genes = set(omics_data['entrez_id'].unique())
            candidate_ranking = sorted(candidate_genes)
        else:
            omics_wide_filled = omics_wide.fillna(0)
            scaler = StandardScaler()
            omics_scaled = scaler.fit_transform(omics_wide_filled)
            pca = PCA()
            pca.fit(omics_scaled)
            loadings = pd.DataFrame(
                pca.components_.T,
                index=omics_wide_filled.columns,
                columns=[f'PC{i+1}' for i in range(pca.n_components_)]
            )
            explained_variance_ratio_cumulative = np.cumsum(pca.explained_variance_ratio_)
            num_pcs = np.where(explained_variance_ratio_cumulative >= 0.95)[0][0] + 1
            selected_pcs = [f'PC{i+1}' for i in range(num_pcs)]
            top_genes_set = set()
            top_n_genes_per_pc = 10
            for pc in selected_pcs:
                abs_loadings = loadings[pc].abs()
                top_genes_pc = abs_loadings.sort_values(ascending=False).head(top_n_genes_per_pc).index.tolist()
                top_genes_set.update(top_genes_pc)
            # Do not truncate yet; compute full ranking based on average absolute loadings over selected PCs.
            avg_loadings = loadings.loc[list(top_genes_set), selected_pcs].abs().mean(axis=1)
            candidate_ranking = avg_loadings.sort_values(ascending=False).index.tolist()
            candidate_genes = set(candidate_ranking)
            omics_data = omics_data[omics_data['entrez_id'].isin(candidate_genes)]
    
    elif args.gene_selection == "mutation":
        if hasattr(cd_data, 'mutations'):
            mutations_df = cd_data.mutations.copy()
            non_silent_mutations = mutations_df[mutations_df['variant_classification'] != 'Silent']
            mutation_counts = non_silent_mutations['entrez_id'].value_counts()
            candidate_ranking = mutation_counts.sort_values(ascending=False).index.tolist()
            candidate_genes = set(candidate_ranking)
            omics_data = omics_data[omics_data['entrez_id'].isin(candidate_genes)]
        else:
            print("Warning: Mutations data not found. Proceeding without mutation-based gene selection.")
            candidate_genes = set(omics_data['entrez_id'].unique())
            candidate_ranking = sorted(candidate_genes)
    
    elif args.gene_selection == "variance":
        try:
            omics_wide = cd.format(cd_data, args.omics, "wide")
        except Exception as e:
            print(f"Error formatting {args.omics} data for variance selection: {e}. Skipping variance-based gene selection.")
            omics_wide = pd.DataFrame()
        if omics_wide.empty:
            print("Warning: Formatted omics data is empty. Skipping variance-based gene selection.")
            candidate_genes = set(omics_data['entrez_id'].unique())
            candidate_ranking = sorted(candidate_genes)
        else:
            gene_variances = omics_wide.var(axis=0)
            # Use the full ranking (do not truncate here)
            candidate_ranking = gene_variances.sort_values(ascending=False).index.tolist()
            candidate_genes = set(candidate_ranking)
            omics_data = omics_data[omics_data['entrez_id'].isin(candidate_genes)]
    
    # Update the cd_data object with the candidate-filtered omics data.
    setattr(cd_data, args.omics, omics_data)
    
    # --------------------------
    # Ensure that experiments and the omics data have matching sample IDs
    # --------------------------
    common_ids = set(cd_data.experiments['improve_sample_id'].unique()).intersection(
        set(omics_data['improve_sample_id'].unique())
    )
    cd_data.experiments = cd_data.experiments[cd_data.experiments['improve_sample_id'].isin(common_ids)].reset_index(drop=True)
    omics_data = omics_data[omics_data['improve_sample_id'].isin(common_ids)]
    setattr(cd_data, args.omics, omics_data)
    
    # Final Data Pre-Processing on experiments (unchanged)
    cd_data.experiments.dropna(inplace=True)
    cd_data.experiments = cd_data.experiments[cd_data.experiments.dose_response_metric == "fit_auc"]
    
    # --------------------------
    # Split the data first, then determine overlapping genes and perform final gene-number selection
    # --------------------------
    if args.test_type == "self":
        split = cd_data.train_test_validate(
            split_type=split_method,
            ratio=[8, 1, 1],
            random_state=data_split_seed,
            stratify_by='fit_auc'
        )
        split_sets = [split.train, split.test, split.validate]
    else:
        split = cd_data.split_train_other(
            split_type=split_method,
            ratio=[8, 2],
            random_state=data_split_seed,
            stratify_by='fit_auc'
        )
        split.validate = split.other
        split_sets = [split.train, split.validate]
    
    # Process each split to get the wide-format omics DataFrame.
    # We do not drop columns here so we can compute the common intersection.
    omics_dfs = {}
    for exp_set in split_sets:
        formatted_omics = exp_set.format(data_type=args.omics, shape='wide', inplace=True)
        formatted_omics = np.log1p(formatted_omics)
        scaler = StandardScaler()
        scaled = scaler.fit_transform(formatted_omics)
        scaled_df = pd.DataFrame(scaled, index=formatted_omics.index, columns=formatted_omics.columns)
        omics_dfs[exp_set] = scaled_df

    # Determine the overlapping genes (columns) across all splits.
    common_genes = set.intersection(*(set(df.columns) for df in omics_dfs.values()))
    common_genes = sorted(common_genes)
    print("Overlapping genes across splits:", common_genes)
    
    # Intersect candidate genes with overlapping genes.
    final_candidate = list(candidate_genes.intersection(common_genes))
    # If a candidate ranking exists, sort final_candidate accordingly.
    if candidate_ranking is not None:
        # Preserve the order from candidate_ranking.
        final_candidate = [gene for gene in candidate_ranking if gene in final_candidate]
    # Finally, if more genes remain than desired, truncate.
    if len(final_candidate) > gene_number:
        final_candidate = final_candidate[:gene_number]
    print(f"Final gene set (top {gene_number} genes):", final_candidate)
    
    # Now, for each split, subset the omics DataFrame to the final selected genes.
    for exp_set in split_sets:
        common_df = omics_dfs[exp_set][final_candidate]
        setattr(exp_set, args.omics, common_df)
    
    # --------------------------
    # Process the experiments for each split and create DataLoaders
    # --------------------------
    if args.test_type == "self":
        for exp_set in split_sets:
            print(f"Processing experiments for split: {exp_set}")
            exp_set.experiments = exp_set.format(data_type='experiments', shape='wide', metrics=[dose_response_metric])
            exp_set.experiments = pd.merge(exp_set.experiments, exp_set.drugs, on="improve_drug_id", how='left')
            exp_set.experiments = exp_set.experiments[["improve_sample_id", "improve_drug_id", dose_response_metric, "canSMILES"]]
            
            if encoder in ["gnn", "morganfp", "descriptor"]:
                exp_set.experiments.rename(columns={"canSMILES": "smiles"}, inplace=True)
                
            n_descriptors = None
            feature_names = None
            if encoder in ["descriptor"]:
                exp_set.drug_descriptors = pd.merge(exp_set.drug_descriptors,
                                                    exp_set.drugs,
                                                    on="improve_drug_id",
                                                    how='left')
                exp_set.drug_descriptors = exp_set.drug_descriptors[["improve_sample_id", "canSMILES", "structural_descriptor", "descriptor_value"]]
                exp_set.drug_descriptors.rename(columns={"canSMILES": "smiles"}, inplace=True)
                exp_set.drug_descriptors = exp_set.drug_descriptors[exp_set.drug_descriptors["structural_descriptor"].str.startswith("n")]
                exp_set.drug_descriptors = exp_set.drug_descriptors[~exp_set.drug_descriptors["structural_descriptor"].str.contains("Ring")]
                exp_set.drug_descriptors.drop(columns="improve_sample_id", inplace=True)
                exp_set.drug_descriptors = exp_set.drug_descriptors.pivot_table(
                    index=["smiles"],
                    columns="structural_descriptor",
                    values="descriptor_value",
                    aggfunc='first'
                ).reset_index()
                features = exp_set.drug_descriptors
                n_descriptors = features.shape[1] - 1
                feature_names = features.drop(['smiles'], axis=1).columns.tolist()
                exp_set.experiments = pd.merge(exp_set.experiments, features, on='smiles', how='left')
                n_columns = [col for col in exp_set.experiments.columns if col.startswith("n")]
                exp_set.experiments[n_columns] = exp_set.experiments[n_columns].astype(float)
        
            exp_set.experiments.drop_duplicates(ignore_index=True, inplace=True)
            
        # Create DataLoaders for train, test, and validate splits.
        loaders = {}
        for name, exp_set in zip(['train', 'test', 'validate'], split_sets):
            data_creater = CreateData(
                gexp=exp_set.transcriptomics,
                encoder_type=encoder,
                metric=dose_response_metric,
                data_path="shared_input/",
                feature_names=feature_names
            )
            data = data_creater.create_data(exp_set.experiments)
            loaders[name] = DataLoader(data, batch_size=bs, shuffle=False, drop_last=False)
        
        train_loader = loaders['train']
        test_loader = loaders['test']
        validate_loader = loaders['validate']
        
    else:
        for exp_set in split_sets:
            print(f"Processing experiments for split: {exp_set}")
            exp_set.experiments = exp_set.format(data_type='experiments', shape='wide', metrics=[dose_response_metric])
            exp_set.experiments = pd.merge(exp_set.experiments, exp_set.drugs, on="improve_drug_id", how='left')
            exp_set.experiments = exp_set.experiments[["improve_sample_id", "improve_drug_id", dose_response_metric, "canSMILES"]]
            if encoder in ["gnn", "morganfp", "descriptor"]:
                exp_set.experiments.rename(columns={"canSMILES": "smiles"}, inplace=True)
                
            n_descriptors = None
            feature_names = None
            if encoder in ["descriptor"]:
                exp_set.drug_descriptors = pd.merge(exp_set.drug_descriptors,
                                                    exp_set.drugs,
                                                    on="improve_drug_id",
                                                    how='left')
                exp_set.drug_descriptors = exp_set.drug_descriptors[["improve_sample_id", "canSMILES", "structural_descriptor", "descriptor_value"]]
                exp_set.drug_descriptors.rename(columns={"canSMILES": "smiles"}, inplace=True)
                exp_set.drug_descriptors = exp_set.drug_descriptors[exp_set.drug_descriptors["structural_descriptor"].str.startswith("n")]
                exp_set.drug_descriptors = exp_set.drug_descriptors[~exp_set.drug_descriptors["structural_descriptor"].str.contains("Ring")]
                exp_set.drug_descriptors.drop(columns="improve_sample_id", inplace=True)
                exp_set.drug_descriptors = exp_set.drug_descriptors.pivot_table(
                    index=["smiles"],
                    columns="structural_descriptor",
                    values="descriptor_value",
                    aggfunc='first'
                ).reset_index()
                features = exp_set.drug_descriptors
                n_descriptors = features.shape[1] - 1
                feature_names = features.drop(['smiles'], axis=1).columns.tolist()
                exp_set.experiments = pd.merge(exp_set.experiments, features, on='smiles', how='left')
                n_columns = [col for col in exp_set.experiments.columns if col.startswith("n")]
                exp_set.experiments[n_columns] = exp_set.experiments[n_columns].astype(float)
    
            exp_set.experiments.drop_duplicates(ignore_index=True, inplace=True)
        
        # Create DataLoaders for train and validate splits.
        loaders = {}
        for name, exp_set in zip(['train', 'validate'], split_sets):
            data_creater = CreateData(
                gexp=exp_set.transcriptomics,
                encoder_type=encoder,
                metric=dose_response_metric,
                data_path="shared_input/",
                feature_names=feature_names
            )
            data = data_creater.create_data(exp_set.experiments)
            loaders[name] = DataLoader(data, batch_size=bs, shuffle=False, drop_last=False)
        
        train_loader = loaders['train']
        validate_loader = loaders['validate']
        
    # --------------------------
    # Set up the model, optimizer, loss, and early stopping
    # --------------------------
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    n_genes = len(getattr(split.train, args.omics).columns)
    if encoder == "gnn":
        model = Model(gnn_features=65, n_descriptors=n_descriptors, encoder_type=encoder, n_genes=n_genes).to(device)
    else:
        model = Model(gnn_features=None, n_descriptors=n_descriptors, encoder_type=encoder, n_genes=n_genes).to(device)
    
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    criterion = nn.MSELoss()
    early_stopping = EarlyStopping(patience=10, verbose=True, chkpoint_name=ckpt_path)
    
    # --------------------------
    # Training loop
    # --------------------------
    history = {"train_loss": [], "val_loss": []}
    hist = {"train_rmse": [], "val_rmse": []}
    for epoch in range(0, n_epochs):
        print(f"Training Epoch {epoch} started")
        model.train()
        loss_all = 0
        print("Training model ...")
        for data in train_loader:
            data = data.to(device)
            optimizer.zero_grad()
            output = model(data)
            output = output.reshape(-1,)
            loss = criterion(output, data.y)
            loss.backward()
            optimizer.step()
            loss_all += loss.item()
        print("Calculating training loss")
        train_loss = loss_all / len(train_loader) if len(train_loader) > 0 else 0
        history["train_loss"].append(train_loss)
        
        val_rmse, _, _, _, _, _, _ = test_fn(validate_loader, model, device)
        history["val_loss"].append(val_rmse)
        print("Determining if early stop should occur")
        early_stopping(val_rmse, model)
        if early_stopping.early_stop:
            print("Early stopping")
            break
        
        hist["val_rmse"].append(val_rmse)
        print(f'Epoch: {epoch}, Val_rmse: {val_rmse:.3f}')
        
    # --------------------------
    # Load best model (checkpoint) and test (if in self-test mode)
    # --------------------------
    model.load_state_dict(torch.load(ckpt_path))

    if args.test_type == "self":
        mse, mae, r2, pearson_corr, spearman_corr, _, _ = test_fn(test_loader, model, device)
        test_rmse = np.sqrt(mse)
        print(f"Test results -- RMSE: {test_rmse:.3f}, MAE: {mae:.3f}, R²: {r2:.3f}, Pearson: {pearson_corr:.3f}, Spearman: {spearman_corr:.3f}")
    else:
        test_rmse, mae, r2, pearson_corr, spearman_corr = None, None, None, None, None

    # --------------------------
    # Write (append) results to file (logging all arguments)
    # --------------------------
    header = ("Seed\tDataset\tEpochs\tSplit Type\tEncoder\tTest Type\tGene Selection\tGene Number\tOmics\t"
              "Test RMSE\tTest MAE\tTest R²\tPearson Correlation\tSpearman Correlation\n")
    
    results_file_path = f'results/seed_{data_split_seed}_epoch_{n_epochs}_{args.dataset}_{split_method}_{encoder}_{args.test_type}_{args.gene_selection}_{gene_number}_{args.omics}_train_results_table.txt'
    
    with open(results_file_path, 'a') as file:
        if os.stat(results_file_path).st_size == 0:
            file.write(header)
        if args.test_type == "self":
            file.write(f"{data_split_seed}\t{args.dataset}\t{n_epochs}\t{split_method}\t{encoder}\t{args.test_type}\t"
                       f"{args.gene_selection}\t{gene_number}\t{args.omics}\t{test_rmse:.3f}\t{mae:.3f}\t{r2:.3f}\t"
                       f"{pearson_corr:.3f}\t{spearman_corr:.3f}\n")
        else:
            file.write(f"{data_split_seed}\t{args.dataset}\t{n_epochs}\t{split_method}\t{encoder}\t{args.test_type}\t"
                       f"{args.gene_selection}\t{gene_number}\t{args.omics}\tNA\tNA\tNA\tNA\tNA\n")
    print(f"Results saved to {results_file_path}")

    with open(args.output, 'a') as out_file:
        if os.stat(args.output).st_size == 0:
            out_file.write(header)
        if args.test_type == "self":
            out_file.write(f"{data_split_seed}\t{args.dataset}\t{n_epochs}\t{split_method}\t{encoder}\t{args.test_type}\t"
                           f"{args.gene_selection}\t{gene_number}\t{args.omics}\t{test_rmse:.3f}\t{mae:.3f}\t{r2:.3f}\t"
                           f"{pearson_corr:.3f}\t{spearman_corr:.3f}\n")
        else:
            out_file.write(f"{data_split_seed}\t{args.dataset}\t{n_epochs}\t{split_method}\t{encoder}\t{args.test_type}\t"
                           f"{args.gene_selection}\t{gene_number}\t{args.omics}\tNA\tNA\tNA\tNA\tNA\n")

if __name__ == '__main__':
    main()
