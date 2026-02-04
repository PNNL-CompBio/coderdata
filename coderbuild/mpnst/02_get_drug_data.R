# 02_get_drug_data.R
#!/usr/bin/env Rscript

# Combined Drug List Extraction for MPNST & MPNST-PDX

library(data.table)
library(dplyr)
library(stringr)
library(synapser)
library(reticulate)

# ----------------------------
# DEBUG helper
# ----------------------------
debug_file_summary <- function(path, label) {
  message("\n========== DEBUG: ", label, " ==========")
  message("Path: ", path)
  if (is.null(path) || is.na(path) || !file.exists(path)) {
    message("File missing.")
    return(invisible(NULL))
  }

  df <- tryCatch(fread(path, sep="\t", header=TRUE), error=function(e) {
    message("fread failed: ", conditionMessage(e))
    return(NULL)
  })
  if (is.null(df)) return(invisible(NULL))

  message("Rows: ", nrow(df), "  Cols: ", ncol(df))
  message("Columns: ", paste(names(df), collapse=", "))

  if ("chem_name" %in% names(df)) {
    cn <- unique(tolower(trimws(as.character(df$chem_name))))
    message("Unique chem_name (first 25): ", paste(head(cn, 25), collapse=", "))
    message("Has chem_name == 'verteporfin'? ", any(cn == "verteporfin", na.rm = TRUE))
    message("Has chem_name == 'wu713d62n9'? ", any(cn == "wu713d62n9", na.rm = TRUE))
  }
  if ("pubchem_id" %in% names(df)) {
    pc <- unique(trimws(as.character(df$pubchem_id)))
    message("Unique pubchem_id (first 10): ", paste(head(pc, 10), collapse=", "))
  }

  invisible(df)
}

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

# 2) PDX-sourced drugs via annotations
pdx_df <- manifest %>%
  select(common_name, PDX_Drug_Data) %>%
  distinct() %>%
  filter(!is.na(PDX_Drug_Data))

pdx_ids <- unique(unlist(strsplit(pdx_df$PDX_Drug_Data, ",")))
pdx_ids <- pdx_ids[ pdx_ids != "" & !is.na(pdx_ids) & pdx_ids != "NA" ]

get_pdx_drugs <- function(synid) {
  q <- sprintf(
    "select experimentalCondition from syn21993642 where id='%s'",
    synid
  )
  df <- synTableQuery(q)$asDataFrame()
  if (nrow(df)==0) return(character(0))
  conds <- unlist(strsplit(df$experimentalCondition, ";"))
  tolower(conds[conds!=""])
}

pdx_drugs <- unique(unlist(lapply(pdx_ids, get_pdx_drugs)))
pdx_drugs <- setdiff(pdx_drugs, "control")

# 3) MicroTissue-sourced drugs via table "children"
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

# 3.5) IMPROVE-sourced drugs via improve_sample_id (syn69801348)
improve_sample_drugs <- character(0)
improve_path <- tryCatch(synGet("syn69801348")$path, error=function(e) NA)
if (!is.na(improve_path) && file.exists(improve_path)) {
  improve_tab <- fread(improve_path, sep="\t", header=TRUE, fill=TRUE)
  if ("improve_sample_id" %in% names(improve_tab)) {
    ids <- improve_tab$improve_sample_id %>% as.character() %>% unique()
    ids <- ids[!is.na(ids) & ids != ""]
    drug_part <- sub("^.*_", "", ids)
    drug_tokens <- unlist(strsplit(drug_part, "\\+"))
    improve_sample_drugs <- unique(tolower(trimws(drug_tokens)))
    improve_sample_drugs <- improve_sample_drugs[improve_sample_drugs != ""]
  }
}

# 3.6) Add drugs from treatment lookup list (syn64608501) Full_name column
lookup_drugs <- character(0)
lookup_path <- tryCatch(synGet("syn64608501")$path, error=function(e) NA)
if (!is.na(lookup_path) && file.exists(lookup_path)) {
  lookup_dt <- tryCatch(
    fread(lookup_path, sep = "\t", header = TRUE, fill = TRUE),
    error = function(e) NULL
  )

  if (is.null(lookup_dt)) {
    lookup_dt <- tryCatch(
      as.data.table(readxl::read_excel(lookup_path)),
      error = function(e) NULL
    )
  }

  if (!is.null(lookup_dt)) {
    full_col <- NULL
    if ("Full_name" %in% names(lookup_dt)) full_col <- "Full_name"
    if (is.null(full_col) && "Full Name" %in% names(lookup_dt)) full_col <- "Full Name"

    if (!is.null(full_col)) {
      lookup_drugs <- unique(tolower(trimws(as.character(lookup_dt[[full_col]]))))
      lookup_drugs <- lookup_drugs[!is.na(lookup_drugs) & nzchar(lookup_drugs)]
      lookup_drugs <- setdiff(lookup_drugs, c("control", "dimethyl_sulfoxide", "dimethyl sulfoxide", "dmso"))
    } else {
      message("\n========== DEBUG: DRUG LOOKUP LIST (syn64608501) ==========")
      message("WARNING: Could not find Full_name column in lookup file. Columns: ",
              paste(names(lookup_dt), collapse = ", "))
    }
  } else {
    message("\n========== DEBUG: DRUG LOOKUP LIST (syn64608501) ==========")
    message("WARNING: Could not read lookup file at: ", lookup_path)
  }
} else {
  message("\n========== DEBUG: DRUG LOOKUP LIST (syn64608501) ==========")
  message("WARNING: lookup file missing/unavailable for syn64608501")
}

message("\n========== DEBUG: DRUG LOOKUP LIST (syn64608501) ==========")
message("lookup_path: ", ifelse(is.na(lookup_path), "<NA>", lookup_path))
message("lookup_drugs count: ", length(lookup_drugs))
if (length(lookup_drugs) > 0) {
  message("First 40 lookup_drugs: ", paste(head(lookup_drugs, 40), collapse = ", "))
}

# 4) Combine and fix bad names
all_drugs <- unique(c(pdx_drugs, mts_drugs, improve_sample_drugs, lookup_drugs))
all_drugs[all_drugs == "pd901"] <- "pd-0325901"
all_drugs[all_drugs == "rmc4630"] <- "rmc-4630"

# IMPORTANT: drop NA/empty
all_drugs <- all_drugs[!is.na(all_drugs) & nzchar(all_drugs)]

# ---- NEW: replace verteporfin with WU713D62N9 so we can do ONE name-based call
# (keep exact casing you provided)
if (any(all_drugs == "verteporfin", na.rm = TRUE)) {
  all_drugs[all_drugs == "verteporfin"] <- "WU713D62N9"
}

message("\n========== DEBUG: DRUG LISTS ==========")
message("all_drugs count: ", length(all_drugs))
message("Contains 'verteporfin' in all_drugs? ", any(all_drugs == "verteporfin", na.rm=TRUE))
message("Contains 'WU713D62N9' in all_drugs? ", any(all_drugs == "WU713D62N9", na.rm=TRUE))
message("First 40 all_drugs: ", paste(head(all_drugs, 40), collapse=", "))

# 5) Read old-drug files or initialize empty
if (!is.na(olddrugfiles)) {
  paths <- strsplit(olddrugfiles, ",")[[1]] %>% trimws()
  old_list <- lapply(paths, function(f) {
    if (file.exists(f)) fread(f, sep="\t", header=TRUE) else {
      warning("Missing old-drug file: ", f)
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
  message("No old-drug files provided; starting fresh")
}

# 6) Write placeholder
fwrite(olddrugs, newdrugfile, sep="\t", quote=FALSE)
message("Wrote placeholder to ", newdrugfile)

# 7) Augment via Python (single call)
ignore_file <- "/tmp/combined_drugs_ignore_chems.txt"

message("\n========== DEBUG: IGNORE FILE ==========")
message("ignore_file: ", ignore_file)
message("ignore_file exists? ", file.exists(ignore_file))
if (file.exists(ignore_file)) {
  ign <- tolower(trimws(readLines(ignore_file, warn=FALSE)))
  message("ignore_file lines (first 50): ", paste(head(ign, 50), collapse=" | "))
  message("ignore contains 'verteporfin'? ", any(ign == "verteporfin"))
  message("ignore contains 'wu713d62n9'? ", any(ign == "wu713d62n9"))
}

use_python("/opt/venv/bin/python3", required=TRUE)
source_python("pubchem_retrieval.py")

prev_paths <- if (!is.na(olddrugfiles)) olddrugfiles else NULL

# reticulate-safe list conversion
name_drugs_py <- as.list(as.character(all_drugs))

if (length(all_drugs) > 0) {
  message("\n========== DEBUG: RUN SINGLE NAME-BASED PYTHON CALL ==========")
  message("Calling update_dataframe_and_write_tsv with isname=TRUE for ", length(all_drugs), " names.")
  message("First 30 all_drugs: ", paste(head(all_drugs, 30), collapse=", "))

  update_dataframe_and_write_tsv(
    unique_names           = name_drugs_py,
    output_filename        = newdrugfile,
    ignore_chems           = ignore_file,
    batch_size             = as.integer(50),
    isname                 = TRUE,
    prev_drug_filepaths    = prev_paths,
    restrict_to_raw_names  = name_drugs_py
  )

  debug_file_summary(newdrugfile, "AFTER SINGLE NAME-BASED PYTHON CALL (newdrugfile)")
} else {
  message("\n========== DEBUG: SKIP PYTHON CALL (all_drugs empty) ==========")
}

# 8) Final unique & save (no CID merge)
message("\n========== DEBUG: MAIN OUTPUT BEFORE FINAL UNIQUE ==========")
debug_file_summary(newdrugfile, "MAIN OUTPUT BEFORE FINAL UNIQUE (newdrugfile)")

tab <- fread(newdrugfile, sep="\t", header=TRUE)
final_tab <- unique(tab)

message("\n========== DEBUG: FINAL CHECK ==========")
message("Final rows after unique(): ", nrow(final_tab))
if ("chem_name" %in% names(final_tab)) {
  cnf <- unique(tolower(trimws(as.character(final_tab$chem_name))))
  message("FINAL has verteporfin? ", any(cnf == "verteporfin", na.rm=TRUE))
  message("FINAL has wu713d62n9? ", any(cnf == "wu713d62n9", na.rm=TRUE))
}

fwrite(final_tab, newdrugfile, sep="\t", quote=FALSE)
message("Wrote full synonyms list to ", newdrugfile)
