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
  "variance"  : Transcriptomic variance-based selection

A new argument, --gene_number, specifies the top X genes to retain (default: 500)
for gene-selection methods that rank genes.
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
        choices=["landmark", "pca", "mutation", "variance"],
        default="landmark",
        help=("Gene selection method: "
              "'landmark' = graphDRP_landmark_genes_map file (default), "
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
    # The output prefix now includes: dataset, split_type, encoder, test_type, gene_selection, gene_number.
    output_prefix = f"{args.dataset}_{split_method}_{encoder}_{args.test_type}_{args.gene_selection}_{gene_number}"
    ckpt_path = f"./tmp/best_{output_prefix}.pt"
    
    # Create necessary directories
    for folder in ["models", "results", "tmp", "shared_input", "plots"]:
        if not os.path.exists(folder):
            os.makedirs(folder)
            print(f"'{folder}' directory created.")

    # --------------------------
    # Load dataset using coderdata and perform gene selection
    # --------------------------
    print(f"Loading dataset: {args.dataset}")
    cd_data = cd.load(args.dataset, "data")
    
    if args.gene_selection == "landmark":
        # Option "landmark": Use graphDRP_landmark_genes_map file (existing method)
        selected_gene_path = "./shared_input/graphDRP_landmark_genes_map.txt"
        if os.path.exists(selected_gene_path):
            try:
                selected_gene_df = pd.read_csv(selected_gene_path, sep='\t')
                if 'To' not in selected_gene_df.columns:
                    raise ValueError("Column 'To' not found in landmark gene file.")
                selected_genes = set(selected_gene_df['To'])
                cd_data.transcriptomics = cd_data.transcriptomics[cd_data.transcriptomics['entrez_id'].isin(selected_genes)]
            except Exception as e:
                print(f"Error reading landmark gene file: {e}. Proceeding without gene filtering.")
        else:
            print(f"Warning: {selected_gene_path} not found. Proceeding without gene filtering.")
            
    elif args.gene_selection == "pca":
        # Option "pca": PCA-based gene selection
        transcriptomics_df = cd_data.transcriptomics.copy()
        # Identify expression columns (excluding identifier columns)
        expr_cols = [col for col in transcriptomics_df.columns if col not in ['improve_sample_id', 'entrez_id']]
        if not expr_cols:
            print("Warning: No expression column found for PCA gene selection. Skipping gene filtering.")
        else:
            try:
                # Use the data object's format function to obtain wide-format data.
                transcriptomics_wide = cd.format(cd_data, "transcriptomics", "wide")
            except Exception as e:
                print(f"Error formatting transcriptomics data: {e}. Skipping PCA gene selection.")
                transcriptomics_wide = pd.DataFrame()
            if transcriptomics_wide.empty:
                print("Warning: Formatted transcriptomics data is empty. Skipping PCA gene selection.")
            else:
                transcriptomics_wide_filled = transcriptomics_wide.fillna(0)
                scaler = StandardScaler()
                transcriptomics_scaled = scaler.fit_transform(transcriptomics_wide_filled)
                pca = PCA()
                pca.fit(transcriptomics_scaled)
                loadings = pd.DataFrame(
                    pca.components_.T,
                    index=transcriptomics_wide_filled.columns,
                    columns=[f'PC{i+1}' for i in range(pca.n_components_)]
                )
                explained_variance_ratio_cumulative = np.cumsum(pca.explained_variance_ratio_)
                num_pcs = np.where(explained_variance_ratio_cumulative >= 0.95)[0][0] + 1
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
                cd_data.transcriptomics = cd_data.transcriptomics[cd_data.transcriptomics['entrez_id'].isin(pca_top_genes)]
                
    elif args.gene_selection == "mutation":
        # Option "mutation": Mutation counts-based gene selection
        if hasattr(cd_data, 'mutations'):
            mutations_df = cd_data.mutations.copy()
            # Filter out 'Silent' mutations
            non_silent_mutations = mutations_df[mutations_df['variant_classification'] != 'Silent']
            mutation_counts = non_silent_mutations['entrez_id'].value_counts()
            mutation_counts_df = mutation_counts.reset_index()
            mutation_counts_df.columns = ['entrez_id', 'mutation_count']
            mutation_counts_df_sorted = mutation_counts_df.sort_values(by='mutation_count', ascending=False)
            mutation_top_genes = mutation_counts_df_sorted.head(gene_number)['entrez_id'].tolist()
            print(f"Top {gene_number} genes with highest mutation frequency:")
            print(mutation_top_genes)
            cd_data.transcriptomics = cd_data.transcriptomics[cd_data.transcriptomics['entrez_id'].isin(mutation_top_genes)]
        else:
            print("Warning: Mutations data not found. Proceeding without mutation-based gene selection.")
            
    elif args.gene_selection == "variance":
        # Option "variance": Transcriptomic variance-based gene selection
        transcriptomics_df = cd_data.transcriptomics.copy()
        expr_cols = [col for col in transcriptomics_df.columns if col not in ['improve_sample_id', 'entrez_id']]
        if not expr_cols:
            print("Warning: No expression column found for transcriptomic variance selection. Skipping gene filtering.")
        else:
            try:
                transcriptomics_wide = cd.format(cd_data, "transcriptomics", "wide")
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
                cd_data.transcriptomics = cd_data.transcriptomics[cd_data.transcriptomics['entrez_id'].isin(top_genes_by_variance)]
    
    # Ensure that experiments and transcriptomics have matching sample IDs
    common_ids = set(cd_data.experiments['improve_sample_id'].unique()).intersection(
        set(cd_data.transcriptomics['improve_sample_id'].unique())
    )
    cd_data.experiments = cd_data.experiments[cd_data.experiments['improve_sample_id'].isin(common_ids)].reset_index(drop=True)
    cd_data.transcriptomics = cd_data.transcriptomics[cd_data.transcriptomics['improve_sample_id'].isin(common_ids)]
    
    # Final Data Pre-Processing - this should be fixed on the data build side eventually.
    cd_data.experiments.dropna(inplace=True)
    cd_data.experiments = cd_data.experiments[cd_data.experiments.dose_response_metric == "fit_auc"]
    
    # --------------------------
    # Process the data and create DataLoaders
    # --------------------------
    if args.test_type == "self":
        # Create a train/validation/test split
        split = cd_data.train_test_validate(
            split_type=split_method,
            ratio=[8, 1, 1],
            random_state=data_split_seed,
            stratify_by='fit_auc'
        )
        for experiment_set in [split.train, split.test, split.validate]:
            print(f"{experiment_set} in progress")
            # Reformat Transcriptomics Data using the data object's format method.
            experiment_set.transcriptomics = experiment_set.format(data_type='transcriptomics', shape='wide', inplace=True)
            # Log Transform Transcriptomics Data
            experiment_set.transcriptomics = np.log1p(experiment_set.transcriptomics)
            # Scale Transcriptomics Data
            scaler = StandardScaler()
            train_gene_exp_scaled = scaler.fit_transform(experiment_set.transcriptomics)
            train_gene_exp_scaled = pd.DataFrame(train_gene_exp_scaled, 
                                                 index=experiment_set.transcriptomics.index, 
                                                 columns=experiment_set.transcriptomics.columns)
            # Remove columns with NaN values
            experiment_set.transcriptomics = train_gene_exp_scaled.dropna(axis=1)
            
            ### Reformat Experiments Data
            experiment_set.experiments = experiment_set.format(data_type='experiments', shape='wide', metrics=[dose_response_metric])
            # Add Canonical Smiles to Experiments Data
            experiment_set.experiments = pd.merge(experiment_set.experiments, experiment_set.drugs, on="improve_drug_id", how='left')
            experiment_set.experiments = experiment_set.experiments[["improve_sample_id", "improve_drug_id", dose_response_metric, "canSMILES"]]
            
            if encoder in ["gnn", "morganfp", "descriptor"]:
                experiment_set.experiments.rename(columns={"canSMILES": "smiles"}, inplace=True)
                
            n_descriptors = None
            feature_names = None
            if encoder in ["descriptor"]:
                experiment_set.drug_descriptors = pd.merge(experiment_set.drug_descriptors,
                                                           experiment_set.drugs,
                                                           on="improve_drug_id",
                                                           how='left')
                experiment_set.drug_descriptors = experiment_set.drug_descriptors[["improve_drug_id", "canSMILES", "structural_descriptor", "descriptor_value"]]
                experiment_set.drug_descriptors.rename(columns={"canSMILES": "smiles"}, inplace=True)
                experiment_set.drug_descriptors = experiment_set.drug_descriptors[experiment_set.drug_descriptors["structural_descriptor"].str.startswith("n")]
                experiment_set.drug_descriptors = experiment_set.drug_descriptors[~experiment_set.drug_descriptors["structural_descriptor"].str.contains("Ring")]
                experiment_set.drug_descriptors.drop(columns="improve_drug_id", inplace=True)
                experiment_set.drug_descriptors = experiment_set.drug_descriptors.pivot_table(
                    index=["smiles"],
                    columns="structural_descriptor",
                    values="descriptor_value",
                    aggfunc='first'
                ).reset_index()
                features = experiment_set.drug_descriptors
                n_descriptors = features.shape[1] - 1
                feature_names = features.drop(['smiles'], axis=1).columns.tolist()
                experiment_set.experiments = pd.merge(experiment_set.experiments, features, on='smiles', how='left')
                n_columns = [col for col in experiment_set.experiments.columns if col.startswith("n")]
                experiment_set.experiments[n_columns] = experiment_set.experiments[n_columns].astype(float)
        
            experiment_set.experiments.drop_duplicates(ignore_index=True, inplace=True)
            
        # Create DataLoaders for train, test, and validation sets
        data_creater = CreateData(gexp=split.train.transcriptomics,
                                  encoder_type=encoder,
                                  metric=dose_response_metric,
                                  data_path="shared_input/",
                                  feature_names=feature_names)
        train_ds = data_creater.create_data(split.train.experiments)
        train_loader = DataLoader(train_ds, batch_size=bs, shuffle=True, drop_last=True)
        test_ds = data_creater.create_data(split.test.experiments)
        test_loader = DataLoader(test_ds, batch_size=bs, shuffle=False, drop_last=False)
        validate_ds = data_creater.create_data(split.validate.experiments)
        validate_loader = DataLoader(validate_ds, batch_size=bs, shuffle=False, drop_last=False)
    else:
        # external test mode: use only train and validation splits.
        split = cd_data.split_train_other(
            split_type=split_method,
            ratio=[8, 2],
            random_state=data_split_seed,
            stratify_by='fit_auc'
        )
        split.validate = split.other
        for experiment_set in [split.train, split.validate]:
            print(f"{experiment_set} in progress")
            experiment_set.transcriptomics = experiment_set.format(data_type='transcriptomics', shape='wide', inplace=True)
            experiment_set.transcriptomics = np.log1p(experiment_set.transcriptomics)
            scaler = StandardScaler()
            train_gene_exp_scaled = scaler.fit_transform(experiment_set.transcriptomics)
            train_gene_exp_scaled = pd.DataFrame(train_gene_exp_scaled, 
                                                 index=experiment_set.transcriptomics.index, 
                                                 columns=experiment_set.transcriptomics.columns)
            experiment_set.transcriptomics = train_gene_exp_scaled.dropna(axis=1)
            experiment_set.experiments = experiment_set.format(data_type='experiments', shape='wide', metrics=[dose_response_metric])
            experiment_set.experiments = pd.merge(experiment_set.experiments, experiment_set.drugs, on="improve_drug_id", how='left')
            experiment_set.experiments = experiment_set.experiments[["improve_sample_id", "improve_drug_id", dose_response_metric, "canSMILES"]]
            if encoder in ["gnn", "morganfp", "descriptor"]:
                experiment_set.experiments.rename(columns={"canSMILES": "smiles"}, inplace=True)
                
            n_descriptors = None
            feature_names = None
            if encoder in ["descriptor"]:
                experiment_set.drug_descriptors = pd.merge(experiment_set.drug_descriptors,
                                                           experiment_set.drugs,
                                                           on="improve_drug_id",
                                                           how='left')
                experiment_set.drug_descriptors = experiment_set.drug_descriptors[["improve_drug_id", "canSMILES", "structural_descriptor", "descriptor_value"]]
                experiment_set.drug_descriptors.rename(columns={"canSMILES": "smiles"}, inplace=True)
                experiment_set.drug_descriptors = experiment_set.drug_descriptors[experiment_set.drug_descriptors["structural_descriptor"].str.startswith("n")]
                experiment_set.drug_descriptors = experiment_set.drug_descriptors[~experiment_set.drug_descriptors["structural_descriptor"].str.contains("Ring")]
                experiment_set.drug_descriptors.drop(columns="improve_drug_id", inplace=True)
                experiment_set.drug_descriptors = experiment_set.drug_descriptors.pivot_table(
                    index=["smiles"],
                    columns="structural_descriptor",
                    values="descriptor_value",
                    aggfunc='first'
                ).reset_index()
                features = experiment_set.drug_descriptors
                n_descriptors = features.shape[1] - 1
                feature_names = features.drop(['smiles'], axis=1).columns.tolist()
                experiment_set.experiments = pd.merge(experiment_set.experiments, features, on='smiles', how='left')
                n_columns = [col for col in experiment_set.experiments.columns if col.startswith("n")]
                experiment_set.experiments[n_columns] = experiment_set.experiments[n_columns].astype(float)
    
            experiment_set.experiments.drop_duplicates(ignore_index=True, inplace=True)
        
        data_creater = CreateData(gexp=split.train.transcriptomics,
                                  encoder_type=encoder,
                                  metric=dose_response_metric,
                                  data_path="shared_input/",
                                  feature_names=feature_names)
        train_ds = data_creater.create_data(split.train.experiments)
        train_loader = DataLoader(train_ds, batch_size=bs, shuffle=True, drop_last=True)
        validate_ds = data_creater.create_data(split.validate.experiments)
        validate_loader = DataLoader(validate_ds, batch_size=bs, shuffle=False, drop_last=False)

    # --------------------------
    # Set up the model, optimizer, loss, and early stopping
    # --------------------------
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    n_genes = len(split.train.transcriptomics.columns)
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
        
        val_rmse, _, _, _, _ = test_fn(validate_loader, model, device)
        history["val_loss"].append(val_rmse)
        print("Determining if early stop should occur")
        early_stopping(val_rmse, model)
        if early_stopping.early_stop:
            print("Early stopping")
            break
        
        hist["val_rmse"].append(val_rmse)
        print(f'Epoch: {epoch}, Val_rmse: {val_rmse:.3f}')
        if epoch % 33 == 0:
            model_save_path = f'models/{output_prefix}_model_seed_{data_split_seed}_epoch_{epoch}.pt'
            torch.save(model.state_dict(), model_save_path)
            print("Model saved at", model_save_path)
    
    # --------------------------
    # Load best model (checkpoint) and test (if in self-test mode)
    # --------------------------
    model.load_state_dict(torch.load(ckpt_path))
    
    if args.test_type == "self":
        test_rmse, pearson_corr, spearman_corr, _, _ = test_fn(test_loader, model, device)
        print(f"Test results -- RMSE: {test_rmse:.3f}, Pearson: {pearson_corr:.3f}, Spearman: {spearman_corr:.3f}")
    else:
        test_rmse, pearson_corr, spearman_corr = None, None, None
    
    # --------------------------
    # Plot training and validation loss
    # --------------------------
    plot_file = f'plots/seed_{data_split_seed}_epoch_{n_epochs}_{args.dataset}_{split_method}_{encoder}_{args.test_type}_{args.gene_selection}_{gene_number}_auc_train_val_plot.png'
    plt.figure(figsize=(10, 6))
    plt.plot(range(0, len(history["train_loss"])), history["train_loss"], label='Train Loss', linestyle='-', color='b')
    plt.plot(range(0, len(history["val_loss"])), history["val_loss"], label='Val Loss', linestyle='-', color='r')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('Training and Validation Loss Over Epochs')
    plt.legend()
    plt.grid(True)
    plt.savefig(plot_file)
    print(f"Loss plot saved to {plot_file}")
    
    # --------------------------
    # Write (append) results to file
    # --------------------------
    results_file_path = f'results/seed_{data_split_seed}_epoch_{n_epochs}_{args.dataset}_{split_method}_{encoder}_{args.test_type}_{args.gene_selection}_{gene_number}_rna_train_results_table.txt'
    header = ("Seed\tDataset\tSplit Type\tEncoder\tTest Type\tGene Selection\tGene Number\t"
              "Test RMSE\tPearson Correlation\tSpearman Correlation\n")
    with open(results_file_path, 'a') as file:
        if os.stat(results_file_path).st_size == 0:
            file.write(header)
        if args.test_type == "self":
            file.write(f"{data_split_seed}\t{args.dataset}\t{split_method}\t{encoder}\t{args.test_type}\t"
                       f"{args.gene_selection}\t{gene_number}\t{test_rmse:.3f}\t{pearson_corr:.3f}\t{spearman_corr:.3f}\n")
        else:
            file.write(f"{data_split_seed}\t{args.dataset}\t{split_method}\t{encoder}\t{args.test_type}\t"
                       f"{args.gene_selection}\t{gene_number}\tNA\tNA\tNA\n")
    print(f"Results saved to {results_file_path}")
    
    # Also append results to the user-specified output file
    with open(args.output, 'a') as out_file:
        if os.stat(args.output).st_size == 0:
            out_file.write(header)
        if args.test_type == "self":
            out_file.write(f"{data_split_seed}\t{args.dataset}\t{split_method}\t{encoder}\t{args.test_type}\t"
                           f"{args.gene_selection}\t{gene_number}\t{test_rmse:.3f}\t{pearson_corr:.3f}\t{spearman_corr:.3f}\n")
        else:
            out_file.write(f"{data_split_seed}\t{args.dataset}\t{split_method}\t{encoder}\t{args.test_type}\t"
                           f"{args.gene_selection}\t{gene_number}\tNA\tNA\tNA\n")

if __name__ == '__main__':
    main()
