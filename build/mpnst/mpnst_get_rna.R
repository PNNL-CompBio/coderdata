# Load required libraries
library(data.table)
library(biomaRt)
library(synapser)

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if a token was provided
if (length(args) == 0) {
  stop("No token provided. Usage: Rscript my_script.R <PAT>", call. = FALSE)
}

# Set your personal access token
PAT <- args[1]

# TESTING
# PAT <- "eyJ0eXAiOiJKV1QiLCJraWQiOiJXN05OOldMSlQ6SjVSSzpMN1RMOlQ3TDc6M1ZYNjpKRU9VOjY0NFI6VTNJWDo1S1oyOjdaQ0s6RlBUSCIsImFsZyI6IlJTMjU2In0.eyJhY2Nlc3MiOnsic2NvcGUiOlsidmlldyIsImRvd25sb2FkIiwibW9kaWZ5Il0sIm9pZGNfY2xhaW1zIjp7fX0sInRva2VuX3R5cGUiOiJQRVJTT05BTF9BQ0NFU1NfVE9LRU4iLCJpc3MiOiJodHRwczovL3JlcG8tcHJvZC5wcm9kLnNhZ2ViYXNlLm9yZy9hdXRoL3YxIiwiYXVkIjoiMCIsIm5iZiI6MTcwMzYwNzc5MiwiaWF0IjoxNzAzNjA3NzkyLCJqdGkiOiI0NjY1Iiwic3ViIjoiMzM3OTU1NiJ9.AclZ5HYjf3BZv5Kmph2CVMa3cKBNE8tFCAzWx2OiodlRb9iIZ2AQqB8vD8AK7QUMVnMb0kgvmugERhaYfZr_CKmG67z94iJqNDRZeOi5QLw52a95TS9w3KdqSyImSiLyKh3eewCxSldVKbA1S8fqqs0MBcTvU_DcDvLlIKo8jKz2ifFAPnFZmlCgMbPjCdyd6nidbTApCvANq4G_KdUIJTgIuZ_k5r7sGbv_cMnD-cDqm3W76OQ4fzRK8e7kaslZuDMMd1tmZQ_5yOndwlXc3nQrm5wumxQopJZpH_cHaRIbkUKjsiC9eLjoCZTzaDrKg2vePoWaCLq9W6ehif4h9A" # TEST
# setwd("/Users/moon515/Library/CloudStorage/OneDrive-PNNL/Documents/coderdata/build/mpnst")

# Log in to Synapse
synLogin(authToken = PAT)

# Define the Ensembl mart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Path to the directory to save .sf files
path <- "./tmp"
dir.create(path, showWarnings = FALSE)

# Read the sample mapping CSV and genes.csv
samples_df <- fread("synapse_NF-MPNST_samples.csv")

# Remove rows with missing or empty 'RNASeq' values
samples_df <- samples_df[samples_df$RNASeq != ""]

# gene mapping table
genes_df <- fread("../genes.csv")

# Extract the RNASeq Synapse IDs
rna_seq_ids <- samples_df$RNASeq
#print(samples_df)

# Download .sf files from Synapse
for (syn_id in rna_seq_ids) {
  # Define the file path
  file_path <- file.path(path, paste0(syn_id, "_quant.sf"))

  # Check if the file already exists to avoid re-downloading
  if (!file.exists(file_path)) {
                                        # Download the file from Synapse
    print(paste('Getting',syn_id))
    syn_file <- synGet(syn_id)
    file.copy(syn_file$path, file_path)
  }
}

# Extract unique ENSG IDs from the genes_df
unique_ensg_ids <- unique(genes_df$other_id)

# Get the mapping from ENST to ENSG
mapping <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id'),
                 filters = 'ensembl_gene_id',
                 values = unique_ensg_ids,
                 mart = ensembl)

# Initialize a list to store data tables
list_of_data_tables <- list()

# Iterate over each .sf file in the directory
sf_files <- list.files(path, pattern = "_quant.sf$", full.names = TRUE)
for (file_name in sf_files) {
  # Read the .sf file
  sf_df <- fread(file_name)
  # Pre-process sf_df: Remove version numbers from ENST IDs
  sf_df[, Name := gsub("\\..*", "", Name)]

  # Merge sf_df with mapping
  sf_df <- merge(sf_df, mapping, by.x = "Name", by.y = "ensembl_transcript_id")

  # Get the sample ID from the file name
  sample_id <- gsub("_quant\\.sf", "", basename(file_name))

  # Aggregate TPMs by gene ID (ENSG)
  sf_df_aggregated <- sf_df[, .(TPM = sum(TPM)), by = .(ensembl_gene_id)]
  sf_df_aggregated$sample_id <- sample_id
  # Append to the list
  list_of_data_tables[[length(list_of_data_tables) + 1]] <- sf_df_aggregated
}

# Combine all data tables into one
combined_df <- rbindlist(list_of_data_tables)

# Create a mapping table from sample_id to improve_sample_id
sample_id_mapping <- samples_df[, .(sample_id = RNASeq, improve_sample_id = ID)]

# Merge combined_df with sample_id_mapping to get improve_sample_id
combined_df <- merge(combined_df, sample_id_mapping, by = "sample_id", all.x = TRUE, no.dups = FALSE)

# Add fixed values for 'source' and 'study'
combined_df[, source := "synapse"]
combined_df[, study := "MPNST"]

# Rename columns to match the target structure
setnames(combined_df, old = c("ensembl_gene_id", "TPM"), new = c("entrez_id", "transcriptomics"))

# Reorder columns as per the target file structure
combined_df <- combined_df[, .(improve_sample_id, transcriptomics, entrez_id, source, study)]

# Create a unique mapping for entrez_id and ENSG
unique_mapping <- genes_df[, .SD[1], by = .(other_id)]

# Merge with combined_df to update entrez_id
combined_df <- merge(combined_df, unique_mapping, by.x = "entrez_id", by.y = "other_id", all.x = TRUE)

# Rename the columns and keep only the necessary ones
setnames(combined_df, old = "entrez_id", new = "ENSG")
setnames(combined_df, old = "entrez_id.y", new = "entrez_id")

# Reorder columns as per your requirement
combined_df <- combined_df[, .(improve_sample_id, transcriptomics, entrez_id, source, study)]

# Display the head of the updated dataframe
head(combined_df)

# Write the combined data to a CSV file
fwrite(combined_df, "MPNST_RNA_seq.csv")

# Check if the final CSV file exists and has data
if(file.exists("MPNST_RNA_seq.csv") && (file.info("MPNST_RNA_seq.csv")$size > 0)) {
  # Remove the 'tmp' directory and its contents
  unlink(path, recursive = TRUE, force = TRUE)
  cat("The 'tmp' directory has been successfully removed.\n")
} else {
  cat("There was an error in the script execution. 'tmp' directory not removed.\n")
}
