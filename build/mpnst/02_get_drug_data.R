#!/usr/bin/env Rscript

# Combined Drug List Extraction for MPNST & MPNST‑PDX

library(data.table)
library(dplyr)
library(stringr)
library(synapser)
library(reticulate)

# 0) Args & login
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript combined_drug_list.R <newdrugfile.tsv> [old_drugs.tsv,...]", call.=FALSE)
}
newdrugfile  <- args[1]
newdrugfile <- file.path(newdrugfile)
olddrugfiles <- if (length(args)>=2 && nzchar(args[2])) args[2] else NA

token <- Sys.getenv("SYNAPSE_AUTH_TOKEN")
if (token == "") stop("Please set SYNAPSE_AUTH_TOKEN in your environment", call.=FALSE)
synLogin(authToken = token)

# 1) Fetch manifest
manifest <- synTableQuery("select * from syn53503360")$asDataFrame() %>%
  rename(common_name = Sample)

# 2) PDX‑sourced drugs via annotations
pdx_df <- manifest %>%
  select(common_name, PDX_Drug_Data) %>%
  distinct() %>%
  filter(!is.na(PDX_Drug_Data))

pdx_ids <- unique(unlist(strsplit(pdx_df$PDX_Drug_Data, ",")))
pdx_ids <- pdx_ids[ pdx_ids != "" & !is.na(pdx_ids) & pdx_ids != "NA" ]

get_pdx_drugs <- function(synid) {
  # Query the metadata table for this file's experimentalCondition
  q <- sprintf(
    "select experimentalCondition from syn21993642 where id='%s'",
    synid
  )
  df <- synTableQuery(q)$asDataFrame()
  if (nrow(df)==0) return(character(0))
  # Split on semicolon, lowercase and drop empties
  conds <- unlist(strsplit(df$experimentalCondition, ";"))
  tolower(conds[conds!=""])
}

pdx_drugs <- unique(unlist(lapply(pdx_ids, get_pdx_drugs)))
pdx_drugs <- setdiff(pdx_drugs, "control")


# 3) MicroTissue‑sourced drugs via table "children"
mts_df <- manifest %>%
  select(common_name, MicroTissueDrugFolder) %>%
  filter(!is.na(MicroTissueDrugFolder))

mts_ids <- unique(unlist(strsplit(mts_df$MicroTissueDrugFolder, ",")))
mts_ids <- mts_ids[mts_ids != "" & !is.na(mts_ids) & mts_ids != "NA"]

get_mts_drugs <- function(parentId) {
  q <- sprintf("select experimentalCondition from syn21993642 where parentId='%s'", parentId)
  synTableQuery(q)$asDataFrame() %>%
    pull(experimentalCondition) %>%
    unique() %>%
    tolower()
}

mts_drugs <- unique(unlist(lapply(mts_ids, get_mts_drugs)))

# 4) Combine and fix bad names
all_drugs <- unique(c(pdx_drugs, mts_drugs))
all_drugs[all_drugs == "pd901"] <- "pd-0325901"
message("Combined drug list: ", paste(all_drugs, collapse=", "))

# 5) Read old‑drug files or initialize empty
if (!is.na(olddrugfiles)) {
  paths <- strsplit(olddrugfiles, ",")[[1]] %>% trimws()
  old_list <- lapply(paths, function(f) {
    if (file.exists(f)) fread(f, sep="\t", header=TRUE) else {
      warning("Missing old‑drug file: ", f)
      NULL
    }
  })
  old_list <- Filter(Negate(is.null), old_list)
  if (length(old_list) > 0) {
    olddrugs <- unique(rbindlist(old_list, use.names=TRUE, fill=TRUE))
    message("Read ", nrow(olddrugs), " old drug records")
  } else {
    olddrugs <- data.table(
      improve_drug_id=integer(), chem_name=character(),
      pubchem_id=character(), canSMILES=character(),
      InChIKey=character(), formula=character(), weight=numeric()
    )
    message("No valid old data; using empty template")
  }
} else {
  olddrugs <- data.table(
    improve_drug_id=integer(), chem_name=character(),
    pubchem_id=character(), canSMILES=character(),
    InChIKey=character(), formula=character(), weight=numeric()
  )
  message("No old‑drug files provided; starting fresh")
}

# 6) Write placeholder
fwrite(olddrugs, newdrugfile, sep="\t", quote=FALSE)
message("Wrote placeholder to ", newdrugfile)

# 7) Augment via Python
ignore_file <- "/tmp/combined_drugs_ignore_chems.txt"
use_python("/opt/venv/bin/python3", required=TRUE)
# use_python("/Users/jaco059/miniconda3/bin/python3", required=TRUE)

# source_python("build/utils/pubchem_retrieval.py")
source_python("pubchem_retrieval.py")

# prepare prev_drug_filepaths argument for python helper (NULL if none)
prev_paths <- if (!is.na(olddrugfiles)) olddrugfiles else NULL
update_dataframe_and_write_tsv(
  unique_names           = all_drugs,
  output_filename        = newdrugfile,
  ignore_chems           = ignore_file,
  batch_size             = as.integer(50),
  isname                 = TRUE,
  prev_drug_filepaths    = prev_paths,
  restrict_to_raw_names  = all_drugs
)


# 8) Final filter & save
tab       <- fread(newdrugfile, sep="\t", header=TRUE)
final_tab <- unique(tab)
fwrite(final_tab, newdrugfile, sep="\t", quote=FALSE)
message("Wrote full synonyms list to ", newdrugfile)