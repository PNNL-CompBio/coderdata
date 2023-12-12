# Get sample info for LINCS data
# Author: Belinda B. Garana
# Email: belinda.garana@pnnl.gov
# Created 2023-12-11
# Last edit: 2023-12-11

library(readr)
library(dplyr)

#### 1. Read in existing samples.csv file ####
samples = readr::read_csv('https://figshare.com/ndownloader/files/40576103',
                          quote='"')

#### 2. Read in LINCS samples information ####
LINCS.link <- 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE101nnn/GSE101406/suppl/GSE101406%5FBroad%5FLINCS%5Fcell%5Finfo.txt.gz'
LINCS.info <- readr::read_delim(LINCS.link, "\t")

#### 3. Create new IMPROVE IDs for new samples ####
common_name <- LINCS.info[!(LINCS.info$cell_id %in% samples$other_id) &
                          !(LINCS.info$cell_id %in% samples$other_names) &
                          !(LINCS.info$cell_id %in% samples$common_name), ]$cell_id
species <- "Homo sapiens (Human)"
improve_sample_id <- seq(max(samples$improve_sample_id)+1, length.out=length(common_name))
id_source <- "LINCS"
other_names <- NA
model_type <- "cell line"

cancer_type <- c()
other_id <- c()
for (i in 1:length(common_name)) {
  cancer_type[i] <- LINCS.info[LINCS.info$cell_id == common_name[i], ]$primary_site
  other_id[i] <- LINCS.info[LINCS.info$cell_id == common_name[i], ]$provider_catalog_id
}

new.samples <- data.frame(common_name, cancer_type, other_names, species,
                          improve_sample_id, id_source, other_id, model_type)

# replace -666 with NA
new.samples[new.samples$other_id == "-666", ]$other_id <- NA

#### 4. Generate new samples.csv file ####
old.samples <- dplyr::distinct(samples[samples$other_id %in% LINCS.info$cell_id |
                            samples$other_names %in% LINCS.info$cell_id |
                            samples$common_name %in% LINCS.info$cell_id, ])
LINCS.samples <- rbind(old.samples, new.samples)
write.csv(LINCS.samples, "lincs_samples.csv", row.names = FALSE)
