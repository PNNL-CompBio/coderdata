# This script generate a new sample table based on pervious beatAML improved sample ID
# It will take the maximum value of beatAML improved sample ID and continue from ID count from there
# Load required libraries
library(data.table)

main <- fread("mpnst/NF_MPNST_samples.csv")
previous_aml <- fread("beatAML/beataml_samples.csv")
max_id <- max(previous_aml$improve_sample_id)
main$improve_sample_id <- seq(from = max_id + 1, length.out = nrow(main))

synapse_main <- fread("mpnst/synapse_NF-MPNST_samples.csv")
# Step 1: Create a dictionary from 'main'
id_dict <- setNames(main$improve_sample_id, main$other_id)

# Step 2: Update 'ID' in 'synapse_main'
synapse_main$ID <- id_dict[synapse_main$Sample]

# Handling NA values if any mismatch occurs (Optional based on your data integrity)
# If there are NAs generated, you might need to check for unmatched keys
# synapse_main$ID[is.na(synapse_main$ID)] <- -1  # Assign a placeholder like -1 for unmatched rows

# Step 3: Save the updated 'synapse_main'
fwrite(synapse_main, "mpnst/synapse_NF-MPNST_samples.csv")
fwrite(main, "mpnst/NF_MPNST_samples.csv") # updated sample file


