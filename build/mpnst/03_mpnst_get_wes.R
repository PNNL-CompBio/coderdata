# Load required libraries
library(data.table)
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

# Path to the directory to save .sf files
path <- "./tmp_WES"
dir.create(path, showWarnings = FALSE)

# Read the sample mapping CSV
samples_df <- fread("mpnst/synapse_NF-MPNST_samples.csv")

# Remove rows with missing or empty 'Mutations' values
samples_df <- samples_df[samples_df$Mutations != ""]
# Create a sample mapping vector
sample_mapping <- setNames(samples_df$ID, samples_df$Mutations)

# Extract the WES Synapse IDs
wes_seq_ids <- samples_df$Mutations

# Download .maf files from Synapse
for (syn_id in wes_seq_ids) {
  # Define the file path
  file_path <- file.path(path, paste0(syn_id, ".maf"))
  
  # Check if the file already exists to avoid re-downloading
  if (!file.exists(file_path)) {
    # Download the file from Synapse
    syn_file <- synGet(syn_id)
    file.copy(syn_file$path, file_path)
  }
}

# Initialize a list to store data tables
list_of_maf_tables <- list()

# Iterate over each .wes file in the directory
# maf_files <- list.files(path, pattern = ".maf$", full.names = TRUE)
maf_files <- list.files(path, pattern = "\\.maf$", full.names = TRUE)
for (file_name in maf_files) {
  # Read the .maf file
  maf_df <- fread(file_name)

  # Get the sample ID from the file name
  sample_id <- gsub("\\.maf", "", basename(file_name))
  improve_sample_id <- sample_mapping[sample_id]
  
  # Assuming 'Variant_Classification' and 'Hugo_Symbol' are columns of interest
  if (all(c("Variant_Classification", "Hugo_Symbol") %in% colnames(maf_df))) {
    # Convert to long format for the columns of interest
    long_df <- maf_df[, .(improve_sample_id = improve_sample_id,
                         mutations = paste(Feature, HGVSc, sep = ":"),
                         entrez_id = Entrez_Gene_Id,
                         variant_classification = Variant_Classification, 
                         source = "synapse", 
                         study = "MPNST"), 
                     by = .(Variant_Classification, Hugo_Symbol)]
  } else {
    # Handle case where the necessary columns do not exist
    next
  }

  # Append to the list
  list_of_maf_tables[[length(list_of_maf_tables) + 1]] <- long_df
}

# Combine all data tables into one
combined_maf_df <- rbindlist(list_of_maf_tables)

# Reorder the columns
combined_maf_df <- combined_maf_df[, .(improve_sample_id, mutations, entrez_id, variant_classification, source, study)]

# Write the combined data to a CSV file
fwrite(combined_maf_df, "mpnst/MPNST_WES_mutation_seq.csv")

# Check if the final CSV file exists and has data
if(file.exists("mpnst/MPNST_WES_mutation_seq.csv") && (file.info("mpnst/MPNST_WES_mutation_seq.csv")$size > 0)) {
  # Remove the 'tmp' directory and its contents
  unlink(path, recursive = TRUE, force = TRUE)
  cat("The 'tmp_WES' directory has been successfully removed.\n")
} else {
  cat("There was an error in the script execution. 'tmp_WES' directory not removed.\n")
}
