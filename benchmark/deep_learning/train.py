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

A new argument, --gene_number, specifies the top X genes to retain (default: 1000)
for gene-selection methods that rank genes.

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

def set_seed(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    np.random.seed(seed)
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    
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
              "'variance' = Omics variance-based")
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
    
    # Set the seed for the model
    set_seed(data_split_seed)
    
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
    # Compute gene ranking without filtering the global omics data
    # --------------------------
    gene_rank = []
    if args.gene_selection == "landmark":
        selected_gene_path = "./shared_input/graphDRP_landmark_genes_map.txt"
        if os.path.exists(selected_gene_path):
            try:
                selected_gene_df = pd.read_csv(selected_gene_path, sep='\t')
                if 'To' not in selected_gene_df.columns:
                    raise ValueError("Column 'To' not found in landmark gene file.")
                selected_genes = set(selected_gene_df['To'])
                gene_rank = list(selected_genes)  # preserve as list for ordering
            except Exception as e:
                print(f"Error reading landmark gene file: {e}. Proceeding without gene filtering.")
        else:
            print(f"Warning: {selected_gene_path} not found. Proceeding without gene filtering.")
            
    elif args.gene_selection == "random":
        random.seed(args.seed)
        all_genes = omics_data['entrez_id'].unique().tolist()
        if gene_number >= len(all_genes):
            print(f"Requested gene_number {gene_number} is greater than or equal to available genes ({len(all_genes)}). Using all genes.")
            selected_genes = all_genes
        else:
            selected_genes = random.sample(all_genes, gene_number)
        print(f"Random selection established {len(selected_genes)} genes.")
        gene_rank = selected_genes
        
    elif args.gene_selection == "pca":
        temp_df = omics_data.copy()
        expr_cols = [col for col in temp_df.columns if col not in ['improve_sample_id', 'entrez_id']]
        if not expr_cols:
            print("Warning: No expression columns found for PCA gene selection. Skipping gene filtering.")
        else:
            try:
                omics_wide = cd.format(cd_data, args.omics, "wide")
            except Exception as e:
                print(f"Error formatting {args.omics} data: {e}. Skipping PCA gene selection.")
                omics_wide = pd.DataFrame()
            if omics_wide.empty:
                print("Warning: Formatted omics data is empty. Skipping PCA gene selection.")
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
                print(f"PCA-based gene ranking established with {len(pca_top_genes)} genes.")
                gene_rank = pca_top_genes
                
    elif args.gene_selection == "mutation":
        if hasattr(cd_data, 'mutations'):
            mutations_df = cd_data.mutations.copy()
            non_silent_mutations = mutations_df[mutations_df['variant_classification'] != 'Silent']
            mutation_counts = non_silent_mutations['entrez_id'].value_counts()
            mutation_counts_df = mutation_counts.reset_index()
            mutation_counts_df.columns = ['entrez_id', 'mutation_count']
            mutation_counts_df_sorted = mutation_counts_df.sort_values(by='mutation_count', ascending=False)
            mutation_top_genes = mutation_counts_df_sorted.head(gene_number)['entrez_id'].tolist()
            print(f"Mutation-based gene ranking established with top {len(mutation_top_genes)} genes.")
            gene_rank = mutation_top_genes
        else:
            print("Warning: Mutations data not found. Proceeding without mutation-based gene selection.")
            
    elif args.gene_selection == "variance":
        temp_df = omics_data.copy()
        expr_cols = [col for col in temp_df.columns if col not in ['improve_sample_id', 'entrez_id']]
        if not expr_cols:
            print("Warning: No expression columns found for variance selection. Skipping gene filtering.")
        else:
            try:
                omics_wide = cd.format(cd_data, args.omics, "wide")
            except Exception as e:
                print(f"Error formatting {args.omics} data for variance selection: {e}. Skipping variance-based gene selection.")
                omics_wide = pd.DataFrame()
            if omics_wide.empty:
                print("Warning: Formatted omics data is empty. Skipping variance-based gene selection.")
            else:
                gene_variances = omics_wide.var(axis=0)
                top_genes_by_variance = gene_variances.sort_values(ascending=False).head(gene_number).index.tolist()
                print(f"Variance-based gene ranking established with top {len(top_genes_by_variance)} genes.")
                gene_rank = top_genes_by_variance

    # Note: At this point, we have computed gene_rank but we have not filtered the global omics_data.
    # Update the cd_data object with the unfiltered omics data.
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
    
    # Final Data Pre-Processing - this should be fixed on the data build side eventually.
    cd_data.experiments.dropna(inplace=True)
    cd_data.experiments = cd_data.experiments[cd_data.experiments.dose_response_metric == "fit_auc"]
    
    # --------------------------
    # Process the data and create DataLoaders
    # --------------------------
    if args.test_type == "self":
        split = cd_data.train_test_validate(
            split_type=split_method,
            ratio=[8, 1, 1],
            random_state=data_split_seed,
            stratify_by='fit_auc'
        )
        for experiment_set in [split.train, split.test, split.validate]:
            print(f"{experiment_set} in progress")
            formatted_omics = experiment_set.format(data_type=args.omics, shape='wide', inplace=True)
            formatted_omics = np.log1p(formatted_omics)
            scaler = StandardScaler()
            scaled = scaler.fit_transform(formatted_omics)
            scaled_df = pd.DataFrame(scaled, index=formatted_omics.index, columns=formatted_omics.columns)
            formatted_omics = scaled_df.dropna(axis=1)
            setattr(experiment_set, args.omics, formatted_omics)

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
                experiment_set.drug_descriptors = experiment_set.drug_descriptors[["improve_sample_id", "canSMILES", "structural_descriptor", "descriptor_value"]]
                experiment_set.drug_descriptors.rename(columns={"canSMILES": "smiles"}, inplace=True)
                experiment_set.drug_descriptors = experiment_set.drug_descriptors[experiment_set.drug_descriptors["structural_descriptor"].str.startswith("n")]
                experiment_set.drug_descriptors = experiment_set.drug_descriptors[~experiment_set.drug_descriptors["structural_descriptor"].str.contains("Ring")]
                experiment_set.drug_descriptors.drop(columns="improve_sample_id", inplace=True)
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
        
        # --------------------------
        # Perform gene selection after splitting based on overlap across splits and the gene_rank list
        # --------------------------
        split_exps = [split.train, split.test, split.validate]
        common_genes = None
        for exp_set in split_exps:
            current_genes = set(getattr(exp_set, args.omics).columns)
            if common_genes is None:
                common_genes = current_genes
            else:
                common_genes = common_genes.intersection(current_genes)
        print(f"Common genes across splits: {len(common_genes)} found.")
        final_gene_list = [gene for gene in gene_rank if gene in common_genes]
        final_gene_list = final_gene_list[:gene_number]
        print(f"Final selected {len(final_gene_list)} genes based on overlap with gene rank list:")
        print(final_gene_list)
        for exp_set in split_exps:
            data_df = getattr(exp_set, args.omics)
            filtered_df = data_df.loc[:, [col for col in data_df.columns if col in final_gene_list]]
            setattr(exp_set, args.omics, filtered_df)
        
        # Process each split and create loaders
        loaders = {}
        for name, gexp in zip(['train', 'test', 'validate'], [split.train, split.test, split.validate]):
            data_creater = CreateData(
                gexp=gexp.transcriptomics,
                encoder_type=encoder,
                metric=dose_response_metric,
                data_path="shared_input/",
                feature_names=feature_names
            )
            data = data_creater.create_data(gexp.experiments)
            loaders[name] = DataLoader(data, batch_size=bs, shuffle=False, drop_last=False)

        # Unpack loaders
        train_loader = loaders['train']
        test_loader = loaders['test']
        validate_loader = loaders['validate']
        
    else:
        split = cd_data.split_train_other(
            split_type=split_method,
            ratio=[8, 2],
            random_state=data_split_seed,
            stratify_by='fit_auc'
        )
        split.validate = split.other
        for experiment_set in [split.train, split.validate]:
            print(f"{experiment_set} in progress")
            formatted_omics = experiment_set.format(data_type=args.omics, shape='wide', inplace=True)
            formatted_omics = np.log1p(formatted_omics)
            scaler = StandardScaler()
            scaled = scaler.fit_transform(formatted_omics)
            scaled_df = pd.DataFrame(scaled, index=formatted_omics.index, columns=formatted_omics.columns)
            formatted_omics = scaled_df.dropna(axis=1)
            setattr(experiment_set, args.omics, formatted_omics)

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
                experiment_set.drug_descriptors = experiment_set.drug_descriptors[["improve_sample_id", "canSMILES", "structural_descriptor", "descriptor_value"]]
                experiment_set.drug_descriptors.rename(columns={"canSMILES": "smiles"}, inplace=True)
                experiment_set.drug_descriptors = experiment_set.drug_descriptors[experiment_set.drug_descriptors["structural_descriptor"].str.startswith("n")]
                experiment_set.drug_descriptors = experiment_set.drug_descriptors[~experiment_set.drug_descriptors["structural_descriptor"].str.contains("Ring")]
                experiment_set.drug_descriptors.drop(columns="improve_sample_id", inplace=True)
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
            
        # --------------------------
        # Perform gene selection after splitting based on overlap across splits and the gene_rank list
        # --------------------------
        split_exps = [split.train, split.validate]
        common_genes = None
        for exp_set in split_exps:
            current_genes = set(getattr(exp_set, args.omics).columns)
            if common_genes is None:
                common_genes = current_genes
            else:
                common_genes = common_genes.intersection(current_genes)
        print(f"Common genes across splits: {len(common_genes)} found.")
        final_gene_list = [gene for gene in gene_rank if gene in common_genes]
        final_gene_list = final_gene_list[:gene_number]
        print(f"Final selected {len(final_gene_list)} genes based on overlap with gene rank list:")
        print(final_gene_list)
        for exp_set in split_exps:
            data_df = getattr(exp_set, args.omics)
            filtered_df = data_df.loc[:, [col for col in data_df.columns if col in final_gene_list]]
            setattr(exp_set, args.omics, filtered_df)
        
        # Process each split and create loaders
        loaders = {}
        for name, gexp in zip(['train', 'validate'], [split.train, split.validate]):
            data_creater = CreateData(
                gexp=gexp.transcriptomics,
                encoder_type=encoder,
                metric=dose_response_metric,
                data_path="shared_input/",
                feature_names=feature_names
            )
            data = data_creater.create_data(gexp.experiments)
            loaders[name] = DataLoader(data, batch_size=bs, shuffle=False, drop_last=False)

        # Unpack loaders
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
        # if epoch % 33 == 0:
        #     model_save_path = f'models/{output_prefix}_epoch_{epoch}.pt'
        #     torch.save(model.state_dict(), model_save_path)
        #     print("Model saved at", model_save_path)
    # --------------------------
    # Load best model (checkpoint) and test (if in self-test mode)
    # --------------------------
    model.load_state_dict(torch.load(ckpt_path))

    if args.test_type == "self":
        mse, mae, r2, pearson_corr, spearman_corr, target_list, predicted_list = test_fn(test_loader, model, device)
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

    df = pd.DataFrame({
        "improve_sample_id": split.validate.experiments.improve_sample_id.to_list(),
        "improve_drug_id": split.validate.experiments.improve_drug_id.to_list(),
        "target": target_list,
        "predicted": predicted_list,
    })

    # Save the predictions with the same prefix as other results
    predictions_file_path = f"results/{output_prefix}_predictions.csv"
    df.to_csv(predictions_file_path, index_label="index")
    print(f"Predictions saved to {predictions_file_path}")

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



#    header = ("Seed\tDataset\tEpochs\tSplit Type\tEncoder\tTest Type\tGene Selection\tGene Number\tOmics\t"
#               "Test RMSE\tTest MAE\tTest R²\tPearson Correlation\tSpearman Correlation\n")
    
#     results_file_path = f'results/seed_{data_split_seed}_epoch_{n_epochs}_{args.dataset}_{split_method}_{encoder}_{args.test_type}_{args.gene_selection}_{gene_number}_{args.omics}_train_results_table.txt'
   
#     df = pd.DataFrame({
#     "improve_sample_id": split.validate.experiments.improve_sample_id.to_list(),
#     "improve_drug_id": split.validate.experiments.improve_drug_id.to_list(),
#     "target": target_list,
#     "predicted": predicted_list,
#     })



# df.to_csv("test_first_predictions.csv", index_label="index")
#     with open(results_file_path, 'a') as file:
#         if os.stat(results_file_path).st_size == 0:
#             file.write(header)
#         if args.test_type == "self":
#             file.write(f"{data_split_seed}\t{args.dataset}\t{n_epochs}\t{split_method}\t{encoder}\t{args.test_type}\t"
#                        f"{args.gene_selection}\t{gene_number}\t{args.omics}\t{test_rmse:.3f}\t{mae:.3f}\t{r2:.3f}\t"
#                        f"{pearson_corr:.3f}\t{spearman_corr:.3f}\n")
#         else:
#             file.write(f"{data_split_seed}\t{args.dataset}\t{n_epochs}\t{split_method}\t{encoder}\t{args.test_type}\t"
#                        f"{args.gene_selection}\t{gene_number}\t{args.omics}\tNA\tNA\tNA\tNA\tNA\n")
#     print(f"Results saved to {results_file_path}")

#     with open(args.output, 'a') as out_file:
#         if os.stat(args.output).st_size == 0:
#             out_file.write(header)
#         if args.test_type == "self":
#             out_file.write(f"{data_split_seed}\t{args.dataset}\t{n_epochs}\t{split_method}\t{encoder}\t{args.test_type}\t"
#                            f"{args.gene_selection}\t{gene_number}\t{args.omics}\t{test_rmse:.3f}\t{mae:.3f}\t{r2:.3f}\t"
#                            f"{pearson_corr:.3f}\t{spearman_corr:.3f}\n")
#         else:
#             out_file.write(f"{data_split_seed}\t{args.dataset}\t{n_epochs}\t{split_method}\t{encoder}\t{args.test_type}\t"
#                            f"{args.gene_selection}\t{gene_number}\t{args.omics}\tNA\tNA\tNA\tNA\tNA\n")

if __name__ == '__main__':
    main()
