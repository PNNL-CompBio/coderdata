# run remap utils python file to finish this pipeline
# run with dockerfile and see whether the program works as intended

# Load required libraries
library(data.table)
library(dplyr)
library(stringr) # for string manipulation

# auc etc... info
auc_table <- fread("tmp_drug/zzz.out.0")
# annotation table
anno_table <- fread("mpnst/no_combo_manifest_updated_2_5_24.csv")
# read previous sample table
sample_table <- fread("mpnst/NF_MPNST_samples.csv")
sample_table %>% select(other_id, improve_sample_id) %>%  head()

# Join the tables and remove rows with NA in improve_sample_id
output_table <- auc_table %>%
  left_join(sample_table, by = c("improve_sample_id" = "improve_sample_id")) %>%
  select(improve_sample_id, everything()) %>%
  filter(!is.na(improve_sample_id)) %>%   # Remove rows with NA in improve_sample_id
  mutate(Drug = Drug) %>% 
  select(Drug, improve_sample_id, source, study, auc, ic50, ec50, ec50se, einf, hs, aac1, auc1, dss1) # reorder column names

# parse and add time and time_units
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
time_hour_df <- as.data.frame(t(sapply(output_table$STUDY, extract_date_hour)), stringsAsFactors = FALSE)
output_table %>% mutate(time = paste0(time_hour_df$date,"_",time_hour_df$hour),
                        time_unit = "YYMMDD_hours") -> final

# update to improve drug ID
source("./utils/mapDrugsToPubchem.R")
drug_table <- buildDrugTable(final$Drug)
# search and assign
# 1. Aggregate drug_mapping_table to find the most common improve_drug_id for each chem_name
aggregated_drug_table <- as.data.table(drug_table)[, .N, by = .(chem_name, improve_drug_id)]
setkey(aggregated_drug_table, chem_name, N)
aggregated_drug_table <- aggregated_drug_table[, .SD[which.max(N)], by = chem_name]

# 2. Join tables
setkey(final, Drug)
setkey(aggregated_drug_table, chem_name)
joined_dt <- final[aggregated_drug_table, nomatch = 0]

# 3. Assign improved drug ID
final[joined_dt, Drug := improve_drug_id]

# save it as a compressed file
# Assuming your final data frame is named 'final'
fwrite(final, file = "experiments.tsv.gz", sep = "\t", quote = FALSE, compress = "gzip")

# move "experiments.tsv.gz" to "/tmp/experiments.tsv.gz"
# file.rename("experiments.tsv.gz", "/tmp/experiments.tsv.gz")
# move "drugs.tsv.gz" to "drugs.tsv.gz"
# file.rename("drugs.tsv.gz", "/tmp/drugs.tsv.gz")

# improve drug ID to SMI
source("./utils/remapDrugsToSmiles_mpnst.R")

# Rename 'drugs.tsv.gz' to 'MPNST_drugs.tsv.gz'
result1 <- file.rename("drugs.tsv.gz", "mpnst/MPNST_drugs.tsv.gz")

# Rename 'experiments.tsv.gz' to 'MPNST_experiments.tsv.gz'
result2 <- file.rename("experiments.tsv.gz", "mpnst/MPNST_experiments.csv.gz")

# Check if the final output file exists and has data
if(file.exists("mpnst/MPNST_drugs.tsv.gz") && (file.exists("mpnst/MPNST_experiments.csv.gz") > 0)) {
  # Remove the 'tmp' directory and its contents
  unlink("tmp_drug", recursive = TRUE, force = TRUE)
  cat("The 'tmp_drug' directory has been successfully removed.\n")
} else {
  cat("There was an error in the script execution. 'tmp_drug' directory not removed.\n")
}

