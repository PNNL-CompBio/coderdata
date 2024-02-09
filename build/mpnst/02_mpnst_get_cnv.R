# Load required libraries
library(data.table)
library(dplyr)
library(synapser)

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if a token was provided
if (length(args) == 0) {
  stop("No token provided. Usage: Rscript my_script.R <PAT>", call. = FALSE)
}

# Set your personal access token
PAT <- args[1]

# Log in to Synapse
synLogin(authToken = PAT)

# Path to the directory to save .cnn files
path <- "./tmp_CNN"
dir.create(path, showWarnings = FALSE)

# Read the sample mapping CSV and genes.csv
samples_df <- fread("mpnst/synapse_NF-MPNST_samples.csv")

# Remove rows with missing or empty 'Copy_number' values
samples_df <- samples_df[samples_df$Copy_number != ""]
# Create a sample mapping vector
sample_mapping <- setNames(samples_df$ID, samples_df$Copy_number)
coding_genes_df <- fread("knowledge_based_isoform_selection_v1.3.txt")
coding_genes_df %>% filter(coding == "coding") -> coding_genes_df
unique_genes <- unique(coding_genes_df$gene_name)

# Extract the CNN Synapse IDs
cnn_seq_ids <- samples_df$Copy_number

# Download .cnn files from Synapse
for (syn_id in cnn_seq_ids) {
  # Define the file path
  file_path <- file.path(path, paste0(syn_id, ".cnn"))
  
  # Check if the file already exists to avoid re-downloading
  if (!file.exists(file_path)) {
    # Download the file from Synapse
    syn_file <- synGet(syn_id)
    file.copy(syn_file$path, file_path)
  }
}

# Initialize a list to store data tables
list_of_cnn_tables <- list()

# Iterate over each .cnn file in the directory
cnn_files <- list.files(path, pattern = "\\.cnn$", full.names = TRUE)
for (file_name in cnn_files) {
  # Read the .cnn file
  cnn_df <- fread(file_name)

  # Get the sample ID from the file name and improve it
  sample_id <- gsub("\\.cnn", "", basename(file_name))
  improve_sample_id <- sample_mapping[sample_id]

  # Check if 'gene' and 'log2' columns exist
  if (all(c("gene", "log2") %in% colnames(cnn_df))) {
    # Split the 'gene' column into individual genes and create repeated rows
    long_df <- cnn_df[, strsplit(as.character(gene), ","), by = .(chromosome, start, end, depth, log2)]
    filtered_df <- long_df %>% filter(V1 %in% unique_genes) # get only protein coding genes and remove empty gene symbols
    filtered_df <- filtered_df[, .(entrez_id = V1,
                           improve_sample_id = improve_sample_id,
                           CNV = log2,
                           copy_call = NA,
                           source = "synapse",
                           study = "mpnst"),
                       by = .(chromosome, start, end, depth)]
  } else {
    # Handle case where the necessary columns do not exist
    next
  }

  # Append to the list
  list_of_cnn_tables[[length(list_of_cnn_tables) + 1]] <- filtered_df
}

# Combine all data tables into one
combined_cnn_df <- rbindlist(list_of_cnn_tables)

# Reorder the columns
combined_cnn_df <- combined_cnn_df[, .(entrez_id, improve_sample_id, CNV, copy_call, source, study)]

# Convert the 'entrez_id' gene IDs to entrez IDs from gene.csv file.
# gene mapping table
genes_df <- fread("genes.csv")
unique_mapping <- genes_df[, .SD[1], by = .(gene_symbol)]
tmp <- merge(combined_cnn_df, unique_mapping, by.x = "entrez_id", by.y = "gene_symbol", all.x = TRUE)

# Rename the columns and keep only the necessary ones
setnames(tmp, old = "entrez_id", new = "gene_symbol")
setnames(tmp, old = "entrez_id.y", new = "entrez_id")

# Reorder columns as per your requirement
tmp <- tmp[, .(entrez_id, improve_sample_id, CNV, copy_call, source, study)]

#deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
tmp[, copy_call := fcase(
  CNV < 0.5210507, "deep del",
  CNV >= 0.5210507 & CNV < 0.7311832, "het loss",
  CNV >= 0.7311832 & CNV < 1.214125, "diploid",
  CNV >= 1.214125 & CNV < 1.422233, "gain",
  CNV >= 1.422233, "amp",
  default = NA_character_  # Assigns NA for values outside the specified ranges
)]

# Write the combined data to a CSV file
fwrite(tmp, "mpnst/MPNST_cnn_mutation_seq.csv")

# Check if the final CSV file exists and has data
if(file.exists("mpnst/MPNST_cnn_mutation_seq.csv") && (file.info("mpnst/MPNST_cnn_mutation_seq.csv")$size > 0)) {
  # Remove the 'tmp' directory and its contents
  unlink(path, recursive = TRUE, force = TRUE)
  cat("The 'tmp_CNN' directory has been successfully removed.\n")
} else {
  cat("There was an error in the script execution. 'tmp_CNN' directory not removed.\n")
}

