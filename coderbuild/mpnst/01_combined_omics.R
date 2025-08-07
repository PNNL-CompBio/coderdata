#!/usr/bin/env Rscript

# Combined MPNST & MPNST-PDX Data Extraction Script
# This script unifies data extraction for PDX, Tumor, and Xenograft-Derived Organoid samples.

# Load required libraries
library(data.table)
library(synapser)
library(dplyr)
library(tidyr)

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 01_combined_omics.R <PAT> <samples.csv> <genes.csv>", call. = FALSE)
}
PAT      <- args[1]
samples  <- args[2]
genes    <- args[3]

# Log in to Synapse
token <- PAT
synLogin(authToken = token)

# Read sample mapping and gene mapping
samples_df <- fread(samples) %>%
  select(improve_sample_id, common_name, model_type) %>%
  distinct()
genes_df <- fread(genes)

# Subset by model type
pdx_samps   <- filter(samples_df, model_type == "patient derived xenograft")
tumor_samps<- filter(samples_df, model_type == "tumor")
mt_samps    <- filter(samples_df, model_type == "xenograft derived organoid")  # These end up being the same as pdx_samps in the manifest.

# Retrieve manifest table from Synapse
manifest <- synTableQuery("select * from syn53503360")$asDataFrame() %>%
  rename(common_name = Sample)

# Build sample tables
pdx_data <- manifest %>%
  select(common_name, starts_with("PDX")) %>%
  left_join(pdx_samps, by = "common_name") %>%
  select(improve_sample_id, common_name, model_type,
         RNASeq = PDX_RNASeq,
         Mutations = PDX_Somatic_Mutations,
         CopyNumber = PDX_CNV,
         Proteomics = PDX_Proteomics) %>%
  filter(!is.na(improve_sample_id))

tumor_data <- manifest %>%
  select(common_name, starts_with("Tumor")) %>%
  left_join(tumor_samps, by = "common_name") %>%
  select(improve_sample_id, common_name, model_type,
         RNASeq = Tumor_RNASeq,
         Mutations = Tumor_Somatic_Mutations,
         CopyNumber = Tumor_CNV) %>%
  mutate(Proteomics = "") %>%
  filter(!is.na(improve_sample_id))

mt_data <- manifest %>%                     #Note, this is the same as pdx_data but I think we default to "xenograft derived organoid" if present (based on original files)
  select(common_name, starts_with("PDX")) %>%
  left_join(mt_samps, by = "common_name") %>%
  select(improve_sample_id, common_name, model_type,
         RNASeq = PDX_RNASeq,
         Mutations = PDX_Somatic_Mutations,
         CopyNumber = PDX_CNV,
         Proteomics = PDX_Proteomics) %>%
  filter(!is.na(improve_sample_id))

# Combine all sample tables
dcombined <- bind_rows(pdx_data, tumor_data, mt_data) %>% distinct()
print("dcombined:")
print(dcombined)

# Helper to assign study label based on model_type
study_label <- function(type) {
  case_when(
    type == "patient derived xenograft"     ~ "MPNST PDX",
    type == "tumor"                          ~ "MPNST Tumor",
    type == "xenograft derived organoid"     ~ "MPNST PDX MT",
    TRUE                                       ~ "MPNST"
  )
}

# Helper to pick metadata based on sample ID and column
pick_meta <- function(id, column) {
  # columns are  {"Proteomics","RNASeq","Mutations","CopyNumber"}
  if (any(tumor_data[[column]] == id, na.rm = TRUE)) {
    sdf <- tumor_data %>% filter(.data[[column]] == id) %>% slice(1)
  } else if (any(mt_data[[column]] == id, na.rm = TRUE)) {
    sdf <- mt_data    %>% filter(.data[[column]] == id) %>% slice(1)
  } else if (any(pdx_data[[column]] == id, na.rm = TRUE)) {
    sdf <- pdx_data   %>% filter(.data[[column]] == id) %>% slice(1)
  } else {
    return(NULL)
  }
  list(
    sample_id  = sdf$improve_sample_id,
    model_type = sdf$model_type
  )
}

# Safe extraction: only return non-empty data frames
i_safe_extract <- function(df, sample_id, source_val, study_val) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df$improve_sample_id <- sample_id
  df$source            <- source_val
  df$study             <- study_val
  df
}

# 1) Proteomics
proteomics_list <- lapply(
  setdiff(dcombined$Proteomics, c("", NA, "NA")),
  function(id) {
    meta <- pick_meta(id, "Proteomics")
    if (is.null(meta)) return(NULL)

    df <- tryCatch(
      fread(synGet(id)$path) %>%
        rename(gene_symbol = Gene) %>%
        left_join(genes_df, by = "gene_symbol") %>%
        select(entrez_id, proteomics = logRatio) %>%
        filter(!is.na(entrez_id), proteomics != 0) %>%
        distinct(),
      error = function(e) NULL
    )
    i_safe_extract(
      df,
      meta$sample_id,
      "NF Data Portal",
      study_label(meta$model_type)
    )
  }
)
proteomics <- bind_rows(proteomics_list)
fwrite(proteomics, file.path("/tmp", "mpnst_proteomics.csv"))
message("Wrote combined proteomics")


# 2) Transcriptomics (PDX, Tumor, and Organoid / MT which comes from PDX..)
transcriptomics_list <- lapply(
  setdiff(dcombined$RNASeq, c("", NA, "NA")),
  function(id) {
    meta <- pick_meta(id, "RNASeq")
    if (is.null(meta)) return(NULL)

    df <- tryCatch({
      fread(synGet(id)$path) %>%
        separate(Name, into = c("other_id","vers"), sep = "\\.") %>%
        select(-vers) %>%
        left_join(genes_df) %>%
        select(entrez_id, transcriptomics = TPM) %>%
        filter(!is.na(entrez_id), transcriptomics != 0) %>%
        distinct()
    }, error = function(e) NULL)

    i_safe_extract(
      df,
      meta$sample_id,
      "NF Data Portal",
      study_label(meta$model_type)
    )
  }
)
transcriptomics <- bind_rows(transcriptomics_list)
fwrite(transcriptomics, file.path("/tmp", "mpnst_transcriptomics.csv"))
message("Wrote combined transcriptomics")


# 3) Mutations (WES)
wes_list <- lapply(
  setdiff(dcombined$Mutations, c("", NA, "NA")),
  function(id) {
    meta <- pick_meta(id, "Mutations")
    if (is.null(meta)) return(NULL)

    clean_id <- gsub('[\"\\[\\]]', '', id)
    df <- tryCatch(
      fread(synGet(clean_id)$path) %>%
        select(entrez_id = Entrez_Gene_Id,
               mutation               = HGVSc,
               variant_classification = Variant_Classification) %>%
        filter(entrez_id %in% genes_df$entrez_id) %>%
        distinct(),
      error = function(e) NULL
    )

    i_safe_extract(
      df,
      meta$sample_id,
      "NF Data Portal",
      study_label(meta$model_type)
    )
  }
)
wes <- bind_rows(wes_list)
fwrite(wes, file.path("/tmp", "mpnst_mutations.csv"))
message("Wrote combined mutations")


# 4) Copy Number Variation (CNV)
cnv_list <- lapply(
  setdiff(dcombined$CopyNumber, c("", NA, "NA")),
  function(id) {
    meta <- pick_meta(id, "CopyNumber")
    if (is.null(meta)) return(NULL)

    clean_id <- gsub('[\"\\[\\]]', '', id)
    raw <- tryCatch(fread(synGet(clean_id)$path), error = function(e) NULL)
    if (is.null(raw)) return(NULL)

    df_long <- raw %>%
      separate_rows(gene, sep = ",") %>%
      rename(gene_symbol = gene) %>%
      left_join(genes_df, by = "gene_symbol") %>%
      filter(!is.na(entrez_id)) %>%
      select(entrez_id, log2) %>%
      distinct() %>%
      mutate(copy_number = 2^log2) %>%
      select(-log2)

    df <- df_long %>%
      mutate(copy_call = case_when(
        copy_number < 0.5210507 ~ "deep del",
        copy_number < 0.7311832 ~ "het loss",
        copy_number < 1.214125  ~ "diploid",
        copy_number < 1.422233  ~ "gain",
        TRUE                    ~ "amp"
      ))

    i_safe_extract(
      df,
      meta$sample_id,
      "NF Data Portal",
      study_label(meta$model_type)
    )
  }
)
cnv <- bind_rows(cnv_list)
fwrite(cnv, file.path("/tmp", "mpnst_copy_number.csv"))
message("Wrote combined copy number")


message("All combined data files created.")
