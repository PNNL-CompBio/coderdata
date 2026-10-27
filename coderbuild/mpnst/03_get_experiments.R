# This is the original version of the script before adding treated-combo support.
# Not actually called in this branch, but kept for reference.

# #  03_get_experiments.R
# #!/usr/bin/env Rscript

# library(data.table)
# library(synapser)
# library(dplyr)
# library(stringr)
# library(readr)
# library(readxl)
# library(tidyr)

# # ============================================================
# # 03_get_experiments.R
# # - MT (original): unchanged logic, but forces time -> integer
# # - MT (treated): robust parsing + sample duplication mapping
# # - PDX: unchanged
# # ============================================================

# # -----------------------
# # args
# # -----------------------
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) != 4) {
#   stop("Usage: Rscript 03_get_experiments.R <PAT> <samples.csv> <drugfile.tsv> <out_prefix>", call. = FALSE)
# }
# PAT        <- args[1]
# samples    <- args[2]
# drugfile   <- args[3]
# out_prefix <- args[4]

# synLogin(authToken = PAT)

# # -----------------------
# # helpers
# # -----------------------
# norm_str <- function(x) tolower(trimws(as.character(x)))

# # normalize drug strings for joining
# norm_drug <- function(x) {
#   x <- norm_str(x)
#   x <- ifelse(x == "pd901", "pd-0325901", x)
#   x
# }

# # make column names consistent across files
# std_names <- function(nms) {
#   nms <- tolower(nms)
#   nms <- gsub("[^a-z0-9]+", "_", nms)
#   nms <- gsub("_+", "_", nms)
#   nms <- gsub("^_|_$", "", nms)
#   nms
# }

# # choose first existing column from candidates
# pick_col <- function(df, candidates) {
#   candidates <- std_names(candidates)
#   nm <- names(df)
#   hit <- intersect(candidates, nm)
#   if (length(hit) == 0) return(NA_character_)
#   hit[1]
# }

# # safe annotation get
# get_ann <- function(synid) {
#   a <- tryCatch(synGetAnnotations(synid)$annotations, error = function(e) NULL)
#   if (is.null(a)) return(list())
#   a
# }

# # -----------------------
# # load inputs
# # -----------------------
# samples_all <- fread(samples) %>% as.data.frame()

# samples_df <- samples_all %>%
#   dplyr::select(improve_sample_id, common_name, model_type) %>%
#   distinct()

# pdx_samps <- dplyr::filter(samples_df, model_type == "patient derived xenograft")
# mt_samps  <- dplyr::filter(samples_df, model_type == "xenograft derived organoid")

# drug_map <- fread(drugfile) %>%
#   dplyr::select(improve_drug_id, chem_name) %>%
#   distinct() %>%
#   mutate(chem_name = norm_drug(chem_name))

# # -----------------------
# # manifest (unchanged)
# # -----------------------
# manifest <- synTableQuery("select * from syn53503360")$asDataFrame() %>%
#   rename(common_name = Sample) %>%
#   as.data.table()

# extract_date_hour <- function(experiment_id) {
#   pattern <- "(\\d{6})_?(\\d{2,3})?"
#   m <- str_match(experiment_id, pattern)
#   date <- m[,2]; hour <- m[,3]
#   date[is.na(date)] <- NA
#   hour[is.na(hour)] <- 48
#   # force integer hour
#   hour <- suppressWarnings(as.integer(hour))
#   hour[is.na(hour)] <- 48L
#   list(date = date, hour = hour)
# }

# # ============================================================
# # 1) MT ORIGINAL (same logic; ensure time integer)
# # ============================================================

# getDrugDataByParent <- function(parid, sampleId) {
#   q <- sprintf(
#     "select id,name,experimentalCondition,parentId from syn21993642 where parentId='%s'",
#     parid
#   )
#   qtab <- synTableQuery(q)$asDataFrame() %>%
#     filter(!is.na(experimentalCondition), name != "synapse_storage_manifest.csv") %>%
#     select(id, name, experimentalCondition)

#   do.call(rbind, lapply(qtab$id, function(x) {
#     info <- filter(qtab, id == x)
#     d    <- extract_date_hour(info$name)

#     fread(synGet(x)$path) %>%
#       filter(response_type == "percent viability") %>%
#       transmute(
#         improve_sample_id = sampleId,
#         DOSE              = (10^dosage) * 1e6,
#         GROWTH            = response,
#         source            = "NF Data Portal",
#         chem_name         = compound_name,
#         study             = paste0("MT ", d$date, " exp"),
#         time              = as.integer(d$hour)
#       )
#   }))
# }

# mts_map <- manifest %>%
#   select(common_name, MicroTissueDrugFolder) %>%
#   inner_join(mt_samps, by = "common_name") %>%
#   separate_rows(MicroTissueDrugFolder, sep = ",") %>%
#   filter(!is.na(MicroTissueDrugFolder), MicroTissueDrugFolder != "NA") %>%
#   select(improve_sample_id, folder = MicroTissueDrugFolder)

# mt_data <- do.call(rbind, lapply(seq_len(nrow(mts_map)), function(i) {
#   sample_id <- mts_map$improve_sample_id[i]
#   folder    <- mts_map$folder[i]
#   getDrugDataByParent(folder, sample_id)
# }))

# mt_curve_orig <- mt_data %>%
#   mutate(
#     chem_name = norm_drug(chem_name),
#     time      = suppressWarnings(as.integer(time))
#   ) %>%
#   left_join(drug_map, by = "chem_name") %>%
#   filter(!is.na(improve_drug_id)) %>%
#   transmute(
#     source             = source,
#     improve_sample_id  = improve_sample_id,
#     Drug               = improve_drug_id,
#     study              = study,
#     time               = as.integer(time),
#     time_unit          = "hours",
#     DOSE               = DOSE,
#     GROWTH             = GROWTH
#   )

# message("Built MT original curve rows: ", nrow(mt_curve_orig))

# # ============================================================
# # 2) MT TREATED (robust parsing + duplication mapping)
# # ============================================================

# # treated mapping from samples.csv:
# # extract drug_key from other_id using regex (no vectorized grepl patterns)
# treated_lookup <- samples_all %>%
#   filter(!is.na(other_id), grepl("_treated_microtissue$", other_id)) %>%
#   mutate(
#     base_common_name = common_name,
#     other_id_clean   = as.character(other_id),
#     # drug_key = part between first "_" and "_<N>hr_treated_microtissue"
#     drug_key = sub("^.*?_", "", other_id_clean),
#     drug_key = sub("_\\d+hr_treated_microtissue$", "", drug_key),
#     drug_key = norm_drug(drug_key)
#   ) %>%
#   select(base_common_name, drug_key, improve_sample_id) %>%
#   distinct()

# # generic fallback: xenograft derived organoid row for base_common_name (exclude treated rows)
# generic_xdo_lookup <- samples_all %>%
#   filter(
#     model_type == "xenograft derived organoid",
#     !is.na(common_name),
#     (is.na(other_id) | !grepl("_treated_microtissue$", other_id))
#   ) %>%
#   select(base_common_name = common_name, generic_improve_sample_id = improve_sample_id) %>%
#   distinct()

# # enumerate candidate treated MT files from the same index table
# mt_index <- synTableQuery("select id,name,experimentalCondition,parentId from syn21993642")$asDataFrame() %>%
#   filter(name != "synapse_storage_manifest.csv") %>%
#   select(id, name, experimentalCondition)

# # build treated curve rows for one file
# get_treated_curve_rows_one <- function(file_id, fallback_condition = NA_character_) {

#   ann <- get_ann(file_id)

#   # always map sample base name from individualID
#   base_common_name <- NA_character_
#   if (!is.null(ann$individualID) && length(ann$individualID) > 0) {
#     base_common_name <- as.character(ann$individualID[[1]])
#   }
#   if (is.na(base_common_name) || base_common_name == "") return(NULL)

#   # drug for treated sample selection: experimentalCondition (annotation preferred)
#   drug_raw <- NA_character_
#   if (!is.null(ann$experimentalCondition) && length(ann$experimentalCondition) > 0) {
#     drug_raw <- as.character(ann$experimentalCondition[[1]])
#   } else if (!is.na(fallback_condition)) {
#     drug_raw <- as.character(fallback_condition)
#   }
#   if (is.na(drug_raw) || drug_raw == "") return(NULL)

#   drug_key <- norm_drug(drug_raw)

#   # load file + standardize columns
#   pth <- synGet(file_id)$path
#   dt <- suppressWarnings(fread(pth))
#   if (nrow(dt) == 0) return(NULL)

#   setnames(dt, names(dt), std_names(names(dt)))

#   # pick key columns robustly
#   col_resp_type <- pick_col(dt, c("response_type", "response type"))
#   col_resp      <- pick_col(dt, c("response", "value"))
#   col_dose      <- pick_col(dt, c("dosage", "dose", "log_m", "logm"))
#   col_cmpd      <- pick_col(dt, c("compound_name", "compound", "drug", "compoundname"))

#   if (is.na(col_resp_type) || is.na(col_resp) || is.na(col_dose)) return(NULL)

#   # percent viability only
#   dt <- dt[get(col_resp_type) == "percent viability", ]
#   if (nrow(dt) == 0) return(NULL)

#   # map improve_drug_id:
#   # 1) annotation experimentalCondition
#   # 2) fallback to file compound_name if present
#   improve_drug_id <- NA_character_

#   dm1 <- drug_map %>% filter(chem_name == drug_key)
#   if (nrow(dm1) > 0) improve_drug_id <- dm1$improve_drug_id[1]

#   if (is.na(improve_drug_id) && !is.na(col_cmpd)) {
#     chem2 <- norm_drug(dt[[col_cmpd]][1])
#     dm2 <- drug_map %>% filter(chem_name == chem2)
#     if (nrow(dm2) > 0) improve_drug_id <- dm2$improve_drug_id[1]
#   }

#   # if still unmapped drug, skip (can be relaxed later when you provide more mapping rules)
#   if (is.na(improve_drug_id) || improve_drug_id == "") return(NULL)

#   # sample mapping:
#   # - if treated match exists: DUPLICATE across all improve_sample_id matches
#   # - else fallback to generic XDO sample for that individualID
#   treated_matches <- treated_lookup %>%
#     filter(base_common_name == !!base_common_name, drug_key == !!drug_key)

#   if (nrow(treated_matches) == 0) {
#     gen <- generic_xdo_lookup %>% filter(base_common_name == !!base_common_name)
#     if (nrow(gen) == 0) return(NULL)
#     target_ids <- gen$generic_improve_sample_id
#   } else {
#     target_ids <- treated_matches$improve_sample_id
#   }

#   # build output duplicated per target id
#   out <- do.call(rbind, lapply(target_ids, function(sid) {
#     data.frame(
#       source            = "NF Data Portal",
#       improve_sample_id = sid,
#       Drug              = improve_drug_id,
#       study             = "MT treated",
#       time              = 48L,
#       time_unit         = "hours",
#       DOSE              = (10^(as.numeric(dt[[col_dose]]))) * 1e6,
#       GROWTH            = as.numeric(dt[[col_resp]]),
#       stringsAsFactors  = FALSE
#     )
#   }))

#   out
# }

# mt_curve_treated_list <- lapply(seq_len(nrow(mt_index)), function(i) {
#   get_treated_curve_rows_one(
#     file_id = mt_index$id[i],
#     fallback_condition = mt_index$experimentalCondition[i]
#   )
# })

# mt_curve_treated <- do.call(rbind, mt_curve_treated_list)

# if (is.null(mt_curve_treated) || nrow(mt_curve_treated) == 0) {
#   mt_curve_treated <- data.frame(
#     source = character(),
#     improve_sample_id = integer(),
#     Drug = character(),
#     study = character(),
#     time = integer(),
#     time_unit = character(),
#     DOSE = numeric(),
#     GROWTH = numeric(),
#     stringsAsFactors = FALSE
#   )
# }

# mt_curve_treated <- mt_curve_treated %>%
#   mutate(time = as.integer(time))

# message("Built MT treated curve rows: ", nrow(mt_curve_treated))

# # ============================================================
# # 3) Combine MT and fit once
# # ============================================================

# # ensure time is integer on both sides (prevents your bind_rows crash)
# mt_curve_orig    <- mt_curve_orig %>% mutate(time = as.integer(time))
# mt_curve_treated <- mt_curve_treated %>% mutate(time = as.integer(time))

# mt_curve_all <- bind_rows(mt_curve_orig, mt_curve_treated)

# fwrite(mt_curve_all, file.path("/tmp", paste0(out_prefix, "_mt_curve_data.tsv")), sep = "\t")
# message("Wrote MT curve data: /tmp/", out_prefix, "_mt_curve_data.tsv")

# system(sprintf(
#   "/opt/venv/bin/python fit_curve.py --input %s --output %s",
#   paste0("/tmp/", out_prefix, "_mt_curve_data.tsv"),
#   paste0("/tmp/", out_prefix, "_mt_experiments")
# ))

# file.rename(
#   paste0("/tmp/", out_prefix, "_mt_experiments.0"),
#   paste0("/tmp/", out_prefix, "_mt_experiments.tsv")
# )
# message("Wrote MT experiments")

# # ============================================================
# # 4) PDX (UNCHANGED)
# # ============================================================

# pdx_map <- do.call(rbind, lapply(seq_len(nrow(manifest)), function(i) {
#   row <- manifest[i, ]
#   samp <- pdx_samps[pdx_samps$common_name == row$common_name, ]
#   if (nrow(samp)==0 || is.na(row$PDX_Drug_Data) || row$PDX_Drug_Data %in% c("", "NA"))
#     return(NULL)
#   ids <- strsplit(row$PDX_Drug_Data, ",")[[1]]
#   ids <- trimws(ids[ids!=""])
#   data.frame(
#     improve_sample_id = samp$improve_sample_id,
#     child_id          = ids,
#     stringsAsFactors  = FALSE
#   )
# }))

# pdx_meta <- do.call(rbind, lapply(seq_len(nrow(pdx_map)), function(i) {
#   sid <- pdx_map$improve_sample_id[i]
#   cid <- pdx_map$child_id[i]
#   pid <- synGet(cid)$parentId
#   if (is.null(pid) || pid=="") stop("no parentId for ", cid)
#   data.frame(
#     improve_sample_id = sid,
#     child_id          = cid,
#     parentId          = pid,
#     stringsAsFactors  = FALSE
#   )
# }))

# all_pdx <- do.call(rbind, lapply(seq_len(nrow(pdx_meta)), function(i) {
#   m   <- pdx_meta[i, ]
#   pth <- synGet(m$child_id)$path
#   raw <- if (grepl("\\.xlsx?$", pth)) read_xlsx(pth) else read_csv(pth)

#   sec_opts  <- c("compound 2_name", "compound_2_name")
#   drug2_col <- intersect(sec_opts, names(raw))[1]
#   compound2 <- if (!is.na(drug2_col)) raw[[drug2_col]] else NA_character_

#   df <- data.frame(
#     child_id                  = m$child_id,
#     specimen_id               = raw$specimen_id,
#     compound_name             = raw$compound_name,
#     compound_2_name           = compound2,
#     experimental_time_point   = raw$experimental_time_point,
#     experimental_time_point_unit = raw$experimental_time_point_unit,
#     assay_value               = raw$assay_value,
#     stringsAsFactors = FALSE
#   )

#   df <- within(df, {
#     drug1     <- tolower(trimws(compound_name))
#     drug2     <- tolower(trimws(compound_2_name))
#     treatment <- ifelse(
#       is.na(drug1) | drug1 %in% c("", "na", "n/a", "nan"),
#       "control",
#       ifelse(!is.na(drug2) & drug2 != "",
#              paste(drug1, drug2, sep = "+"),
#              drug1
#       )
#     )
#     time      <- experimental_time_point
#     time_unit <- experimental_time_point_unit
#     volume    <- assay_value
#   })

#   df[ , c("child_id", "specimen_id", "treatment", "time", "time_unit", "volume")]
# }))

# pdx_data <- merge(all_pdx, pdx_meta, by="child_id")

# pdx_data <- subset(pdx_data, duplicated(child_id) | TRUE)
# pdx_data <- within(pdx_data, {
#   experiment <- parentId
#   model_id   <- improve_sample_id
# })

# has_ctl     <- tapply(pdx_data$treatment == "control", pdx_data$experiment, any)
# no_ctl_exps <- names(has_ctl)[!has_ctl]
# pdx_data <- pdx_data[pdx_data$experiment %in% names(has_ctl)[has_ctl], ]

# pdx_data <- pdx_data[ , c("experiment","specimen_id","treatment",
#                           "time","time_unit","volume","model_id")]

# pdx_data$treatment <- gsub("doxorubinsin",
#                            "doxorubicin",
#                            pdx_data$treatment,
#                            ignore.case = TRUE)

# pdx_data <- na.omit(pdx_data)

# fwrite(pdx_data, file.path("/tmp", paste0(out_prefix, "_pdx_curve_data.tsv")), sep = "\t")
# message("Wrote PDX curve data")

# system(sprintf(
#   "/opt/venv/bin/python calc_pdx_metrics.py %s --drugfile %s --outprefix %s --source 'NF Data Portal' --study 'MPNST PDX'",
#   paste0("/tmp/", out_prefix, "_pdx_curve_data.tsv"),
#   drugfile,
#   paste0("/tmp/", out_prefix, "_pdx")
# ))

# message("Wrote PDX experiments to /tmp/", out_prefix, "_pdx_experiments.tsv and combinations")

# # ============================================================
# # 5) Combine all Experiments (unchanged)
# # ============================================================

# mt_exp <- fread(paste0("/tmp/", out_prefix, "_mt_experiments.tsv")) %>%
#   mutate(dose_response_value = as.character(dose_response_value))

# pdx_exp <- fread(paste0("/tmp/", out_prefix, "_pdx_experiments.tsv")) %>%
#   mutate(dose_response_value = as.character(dose_response_value))

# all_exp <- bind_rows(mt_exp, pdx_exp)

# fwrite(all_exp, paste0("/tmp/", out_prefix, "_experiments.tsv"), sep = "\t")
# message("Wrote combined experiments: /tmp/", out_prefix, "_experiments.tsv")

# file.rename(
#   paste0("/tmp/", out_prefix, "_pdx_combinations.tsv"),
#   paste0("/tmp/", out_prefix, "_combinations.tsv")
# )
