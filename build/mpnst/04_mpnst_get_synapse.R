# Load required libraries
library(data.table)
library(dplyr)
library(stringr) # for string manipulation
library(synapser)

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if a token was provided
if (length(args) == 0) {
  stop("No token provided. Usage: Rscript my_script.R <PAT>", call. = FALSE)
}

# Set your personal access token
PAT <- args[1]
synLogin(authToken = PAT)

dt <- fread("mpnst/no_combo_manifest_updated_2_5_24.csv")
dt %>% select(id, individualID, experimentId, experimentalCondition, name) -> selected_dt

# Modify the extract_date_hour function to return a named vector
extract_date_hour <- function(experiment_id) {
  pattern <- "(\\d{6})_?(\\d{2,3})?"
  matches <- str_match(experiment_id, pattern)
  date <- matches[, 2]
  hour <- matches[, 3]
  date[is.na(date)] <- NA  # Replace with NA instead of blank
  hour[is.na(hour)] <- 48  # Replace with 48 instead of blank (default)
  return(c(date = date, hour = hour))
}

# Apply the function and create a data frame directly
time_hour_df <- as.data.frame(t(sapply(selected_dt$experimentId, extract_date_hour)), stringsAsFactors = FALSE)

# Check the structure of time_hour_df
print(str(time_hour_df))

# If the structure is correct, add to selected_dt
selected_dt <- cbind(selected_dt, time_hour_df)

main <- fread("mpnst/synapse_NF-MPNST_samples.csv")

exp_input_anno <- data.table(
  Sample = selected_dt$individualID,
  source = "nf_data_portal",
  study = selected_dt$date,
  time = as.numeric(selected_dt$hour),
  time_unit = "hours",
  synapseID = selected_dt$id,
  drug = selected_dt$experimentalCondition
)


# download list of synapseID
# Extract the Synapse IDs
seq_ids <- exp_input_anno$synapseID
path <- "./tmp_drug"
dir.create(path, showWarnings = FALSE)

# Download .csv files from Synapse
for (syn_id in seq_ids) {
  # Define the file path
  file_path <- file.path(path, paste0(syn_id, ".csv"))
  # Check if the file already exists
  if (!file.exists(file_path)) {
    # File does not exist, proceed to download
    syn_file <- synGet(syn_id)
    file.copy(syn_file$path, file_path)
  } else {
    # File already exists, skip download
    cat("File for", syn_id, "already exists. Skipping download.\n")
  }
}

# all downloaded data is available in tmp_drug
# Set the path to the directory
dir_path <- "tmp_drug"

# List all CSV files in the directory
file_list <- list.files(path = dir_path, pattern = "\\.csv$", full.names = TRUE)

# Initialize an empty data frame for storing the combined data
combined_data <- data.frame()
# mapping table for synapseID to ID
map_improved_id <- data.table(improved_id = main$ID,
                              individualID = main$Sample)  
# Ensure that map_improved_id has appropriate column names
setnames(map_improved_id, c("improved_id", "individualID"))

# Process each file
for (file in file_list) {
  # Read the file
  # data <- fread(file)
  cat("processing ", file, "\n")
  data <- fread(file)
  # Extract Synapse ID from the filename
  synapse_id <- str_extract(basename(file), "syn\\d+")
  # Find the experimentId, drug, and individualID in selected_dt
  selected_info <- selected_dt[id == synapse_id, .(individualID, experimentId, experimentalCondition)]
  # Join with map_improved_id to find improved_id
  joined_info <- selected_info[map_improved_id, on = .(individualID), nomatch = 0L]
  
  # # Find the experimentId in selected_dt
  # experiment_id <- selected_dt[id == synapse_id, experimentId]
  # # Find the drug in selected_dt
  # drug <- selected_dt[id == synapse_id, experimentalCondition]
  # # Find the sample in selected_dt
  # id <- selected_dt[id == synapse_id, individualID]
  # # Find improved_id based on id
  # improved_sample <- ??? # use the map_improved_id
  
  # Apply transformations
  output <- data %>% 
    filter(response_type == "percent viability") %>%
    mutate(
        improve_sample_id = joined_info$improved_id,
        DOSE = exp(dosage),
           GROWTH = response / 100,
           source = "nf_data_portal",
           CELL = joined_info$individualID,
           Drug = joined_info$experimentalCondition,
           study = joined_info$experimentId) %>% 
    select(improve_sample_id,DOSE,GROWTH,source,CELL,Drug,study)
  
  # Append the processed data to the combined data frame
  combined_data <- rbind(combined_data, output)
}

# # Process each file
# for (file in file_list) {
#   # Read the file
#   data <- fread(file)
#   # Extract Synapse ID from the filename
#   synapse_id <- str_extract(basename(file), "syn\\d+")
#   # Find the improve_sample_id in selected_dt
#   improve_sample <- selected_dt[id == synapse_id, experimentId]
#   # Find the experimentId in selected_dt
#   experiment_id <- selected_dt[id == synapse_id, experimentId]
#   # Find the drug in selected_dt
#   drug <- selected_dt[id == synapse_id, experimentalCondition]
#   # Find the drug in selected_dt
#   id <- selected_dt[id == synapse_id, individualID]
  
#   # Apply transformations
#   output <- data %>% 
#     filter(response_type == "percent viability") %>%
#     mutate(improve_sample_id = ,
#            DOSE = exp(dosage),
#            GROWTH = response / 100,
#            source = "nf_data_portal",
#            CELL = id,
#            DRUG = drug,
#            STUDY = experiment_id) %>% 
#     select(DOSE,GROWTH,source,CELL,DRUG,STUDY)
  
#   # Append the processed data to the combined data frame
#   combined_data <- rbind(combined_data, output)
# }

# Write the combined data to a TSV file
combined_tsv_file_name <- "combined_data.tsv"
fwrite(combined_data, file.path(dir_path, combined_tsv_file_name), sep = "\t")

# NEXT STEP
# python utils/fit_curve.py --input ./tmp_drug/combined_data.tsv