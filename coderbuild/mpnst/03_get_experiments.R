library(data.table)
library(synapser)
library(dplyr)
library(stringr)
library(readr)
library(readxl)
library(tidyr)

# Check that correct number of arguments are present
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript 03_get_experiments.R <PAT> <samples.csv> <drugfile.tsv> <out_prefix>", call. = FALSE)
}
PAT        <- args[1]
samples    <- args[2]
drugfile   <- args[3]
out_prefix <- args[4]

synLogin(authToken = PAT)

# Read in sampes file
samples_df <- fread(samples) %>%
  select(improve_sample_id, common_name, model_type) %>%
  distinct()

pdx_samps <- filter(samples_df, model_type == "patient derived xenograft")
mt_samps  <- filter(samples_df, model_type == "xenograft derived organoid")

# Get manifest table from Synapse
manifest <- synTableQuery("select * from syn53503360")$asDataFrame() %>%
  rename(common_name = Sample) %>%
  as.data.table()

# Helper Function to extract date and hour from experiment ID
extract_date_hour <- function(experiment_id) {
  pattern <- "(\\d{6})_?(\\d{2,3})?"
  m <- str_match(experiment_id, pattern)
  date <- m[,2]; hour <- m[,3]
  date[is.na(date)] <- NA
  hour[is.na(hour)] <- 48
  list(date = date, hour = hour)
}

# ────────────────────────────────────────────────
# MicroTissue Experiments 
# ────────────────────────────────────────────────

getDrugDataByParent <- function(parid, sampleId) {
  q <- sprintf(
    "select id,name,experimentalCondition,parentId from syn21993642 where parentId='%s'",
    parid
  )
  qtab <- synTableQuery(q)$asDataFrame() %>%
    filter(!is.na(experimentalCondition), name != "synapse_storage_manifest.csv") %>%
    select(id, name, experimentalCondition)
  do.call(rbind, lapply(qtab$id, function(x) {
    info <- filter(qtab, id == x)
    d    <- extract_date_hour(info$name)
    fread(synGet(x)$path) %>%
      filter(response_type == "percent viability") %>%
      transmute(
        improve_sample_id = sampleId,
        DOSE              = (10^dosage) * 1e6,
        GROWTH            = response,
        source            = "NF Data Portal",
        chem_name         = compound_name,
        study             = paste0("MT ", d$date, " exp"),
        time              = d$hour
      )
  }))
}

# Create map of MicroTissue Drug Folders
mts_map <- manifest %>%
  select(common_name, MicroTissueDrugFolder) %>%
  inner_join(mt_samps, by = "common_name") %>%
  separate_rows(MicroTissueDrugFolder, sep = ",") %>%
  # keep exactly what old script did: drop only "NA" and actual NA
  filter(
    !is.na(MicroTissueDrugFolder),
    MicroTissueDrugFolder != "NA"
  ) %>%
  select(
    improve_sample_id,
    folder = MicroTissueDrugFolder
  )

# Fetch all MicroTissue drug response data
mt_data <- do.call(rbind, lapply(seq_len(nrow(mts_map)), function(i) {
  sample_id <- mts_map$improve_sample_id[i]
  folder    <- mts_map$folder[i]
  getDrugDataByParent(folder, sample_id)
}))

drug_map <- fread(drugfile) %>%
  select(improve_drug_id, chem_name) %>%
  distinct()

# Clean up drug names and join with drug_map
mt_curve <- mt_data %>%
  mutate(
    chem_name = tolower(chem_name),
    chem_name = ifelse(chem_name == "pd901", "pd-0325901", chem_name)
  ) %>%
  left_join(drug_map, by = "chem_name") %>%
  filter(!is.na(improve_drug_id)) %>%
  transmute(
    source             = source,
    improve_sample_id  = improve_sample_id,
    Drug               = improve_drug_id,
    study              = study,
    time               = time,
    time_unit          = "hours",
    DOSE               = DOSE,
    GROWTH             = GROWTH
  )

# Run curve fitting, Write MicroTissue curve data
fwrite(mt_curve, file.path("/tmp", paste0(out_prefix, "_mt_curve_data.tsv")), sep = "\t")

message("Wrote MT curve data")

# Write MT experiments file
system(sprintf(
  "/opt/venv/bin/python fit_curve.py --input %s --output %s",
  paste0("/tmp/", out_prefix, "_mt_curve_data.tsv"),
  paste0("/tmp/", out_prefix, "_mt_experiments")
))
file.rename(
  paste0("/tmp/", out_prefix, "_mt_experiments.0"),
  paste0("/tmp/", out_prefix, "_mt_experiments.tsv")
)
message("Wrote MT experiments")

# ────────────────────────────────────────────────
# PDX Experiments
# ────────────────────────────────────────────────

# Create a map of PDX Drug Data
# This will be used to fetch the drug data for each PDX sample
pdx_map <- do.call(rbind, lapply(seq_len(nrow(manifest)), function(i) {
  row <- manifest[i, ]
  samp <- pdx_samps[pdx_samps$common_name == row$common_name, ]
  if (nrow(samp)==0 || is.na(row$PDX_Drug_Data) || row$PDX_Drug_Data %in% c("", "NA"))
    return(NULL)
  ids <- strsplit(row$PDX_Drug_Data, ",")[[1]]
  ids <- trimws(ids[ids!=""])
  data.frame(
    improve_sample_id = samp$improve_sample_id,
    child_id          = ids,
    stringsAsFactors  = FALSE
  )
}))

# Create a dataframe of PDX metadata
pdx_meta <- do.call(rbind, lapply(seq_len(nrow(pdx_map)), function(i) {
  sid <- pdx_map$improve_sample_id[i]
  cid <- pdx_map$child_id[i]
  pid <- synGet(cid)$parentId
  if (is.null(pid) || pid=="") stop("no parentId for ", cid)
  data.frame(
    improve_sample_id = sid,
    child_id          = cid,
    parentId          = pid,
    stringsAsFactors  = FALSE
  )
}))

all_pdx <- do.call(rbind, lapply(seq_len(nrow(pdx_meta)), function(i) {
  m   <- pdx_meta[i, ]
  pth <- synGet(m$child_id)$path
  raw <- if (grepl("\\.xlsx?$", pth)) read_xlsx(pth) else read_csv(pth)

  # detect second‐drug column
  sec_opts  <- c("compound 2_name", "compound_2_name")
  drug2_col <- intersect(sec_opts, names(raw))[1]
  compound2 <- if (!is.na(drug2_col)) raw[[drug2_col]] else NA_character_

  df <- data.frame(
    child_id                  = m$child_id,
    specimen_id               = raw$specimen_id,
    compound_name             = raw$compound_name,
    compound_2_name           = compound2,
    experimental_time_point   = raw$experimental_time_point,
    experimental_time_point_unit = raw$experimental_time_point_unit,
    assay_value               = raw$assay_value,
    stringsAsFactors = FALSE
  )

  df <- within(df, {
    drug1     <- tolower(trimws(compound_name))
    drug2     <- tolower(trimws(compound_2_name))
    treatment <- ifelse(
      is.na(drug1) | drug1 %in% c("", "na", "n/a", "nan"),
      "control",
      ifelse(!is.na(drug2) & drug2 != "",
             paste(drug1, drug2, sep = "+"),
             drug1
      )
    )
    time      <- experimental_time_point
    time_unit <- experimental_time_point_unit
    volume    <- assay_value
  })

  df[ , c("child_id", "specimen_id", "treatment", "time", "time_unit", "volume")]
}))

# join on parentId and sample
pdx_data <- merge(all_pdx, pdx_meta, by="child_id") 

pdx_data <- subset(pdx_data, duplicated(child_id) | TRUE)  
pdx_data <- within(pdx_data, {
  experiment <- parentId
  model_id   <- improve_sample_id
})

# Filter out experiments missing a control
has_ctl     <- tapply(pdx_data$treatment == "control", pdx_data$experiment, any)
no_ctl_exps <- names(has_ctl)[!has_ctl]
pdx_data <- pdx_data[pdx_data$experiment %in% names(has_ctl)[has_ctl], ]

# Reorder final columns
pdx_data <- pdx_data[ , c("experiment","specimen_id","treatment",
                          "time","time_unit","volume","model_id")]

# Correct doxorubinsin typo across all data
pdx_data$treatment <- gsub("doxorubinsin",
                           "doxorubicin",
                           pdx_data$treatment,
                           ignore.case = TRUE)

# Drop any remaining NA rows
pdx_data <- na.omit(pdx_data)

# write & fit
fwrite(pdx_data, file.path("/tmp", paste0(out_prefix, "_pdx_curve_data.tsv")), sep = "\t")


message("Wrote PDX curve data")

system(sprintf(
  "/opt/venv/bin/python calc_pdx_metrics.py %s --drugfile %s --outprefix %s --source 'NF Data Portal' --study 'MPNST PDX'",
  paste0("/tmp/", out_prefix, "_pdx_curve_data.tsv"),
  drugfile,
  paste0("/tmp/", out_prefix, "_pdx")
))



message("Wrote PDX experiments to ", "/tmp/",  out_prefix, "_pdx_experiments.tsv and combinations")


# ────────────────────────────────────────────────
# Combine all Experiments
# ────────────────────────────────────────────────

# Read MicroTissue experiments
mt_exp <- fread(paste0("/tmp/", out_prefix, "_mt_experiments.tsv")) %>%
  mutate(
    dose_response_value = as.character(dose_response_value)
  )

# Read PDX experiments
pdx_exp <- fread(paste0("/tmp/", out_prefix, "_pdx_experiments.tsv")) %>%   
  mutate(
    dose_response_value = as.character(dose_response_value)
  )

# Join experiments into one.
all_exp <- bind_rows(mt_exp, pdx_exp)

# Write out Experiments
fwrite(all_exp, paste0("/tmp/", out_prefix, "_experiments.tsv"), sep = "\t")
message("Wrote combined experiments: /tmp/", out_prefix, "_experiments.tsv")


# Rename the Drug Combination data file to fit schema naming
file.rename(
  paste0("/tmp/", out_prefix, "_pdx_combinations.tsv"),
  paste0("/tmp/", out_prefix, "_combinations.tsv")
)