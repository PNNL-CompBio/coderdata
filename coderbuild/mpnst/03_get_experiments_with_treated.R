#!/usr/bin/env Rscript
# 03_get_experiments_with_treated_combos.R
#
# DEBUG VERSION
# - Adds verbose print statements throughout MT + treated-microtissue code paths + helpers
# - Adds safe readers/guards so empty MT outputs don’t crash downstream (instead they print diagnostics)
#

suppressPackageStartupMessages({
  library(data.table)
  library(synapser)
  library(dplyr)
  library(stringr)
  library(readr)
  library(readxl)
  library(tidyr)
})

# ============================================================
# FAST DEBUG SWITCH - skip ahead to later treated section (save time debugging)
# ============================================================
SKIP_FAST <- FALSE  # set to TRUE to skip MT/PDX/combine sections

# ============================================================
# DEBUG CONTROLS
# ============================================================
DEBUG <- TRUE
DEBUG_MAX_ITEMS <- 10   # limit loop prints
DEBUG_HEAD_N    <- 6

# -----------------------
# args
# -----------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript 03_get_experiments_with_treated_combos.R <PAT> <samples.csv> <drugfile.tsv> <out_prefix>", call. = FALSE)
}
PAT        <- args[1]
samples    <- args[2]
drugfile   <- args[3]
out_prefix <- args[4]

# -----------------------
# helpers
# -----------------------
cat0 <- function(...) cat(..., "\n", sep = "")
hr  <- function() cat0(strrep("─", 80))

dbg_time <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
dbg <- function(...) { if (isTRUE(DEBUG)) cat0(dbg_time(), " | ", paste0(...)) }
dbg_hr <- function(title = NULL) {
  if (!isTRUE(DEBUG)) return(invisible(NULL))
  hr()
  if (!is.null(title)) cat0(dbg_time(), " | ", title)
  hr()
}

dbg_df <- function(df, label = "df") {
  if (!isTRUE(DEBUG)) return(invisible(NULL))
  if (is.null(df)) {
    dbg("[", label, "] is NULL")
    return(invisible(NULL))
  }
  dbg("[", label, "] rows=", nrow(df), " cols=", ncol(df))
  dbg("[", label, "] cols: ", paste(names(df), collapse = ", "))
  if (nrow(df) > 0) print(utils::head(df, DEBUG_HEAD_N))
  invisible(NULL)
}

dbg_file <- function(path, label = "file") {
  if (!isTRUE(DEBUG)) return(invisible(NULL))
  dbg("[", label, "] path=", path)
  if (is.null(path) || is.na(path) || path == "" || !file.exists(path)) {
    dbg("[", label, "] MISSING")
    return(invisible(NULL))
  }
  sz <- file.info(path)$size
  dbg("[", label, "] exists size=", sz)
  if (!is.na(sz) && sz > 0) {
    dt <- tryCatch(data.table::fread(path, nrows = DEBUG_HEAD_N), error = function(e) NULL)
    if (!is.null(dt)) {
      dbg("[", label, "] peek rows=", nrow(dt), " cols=", ncol(dt))
      dbg("[", label, "] cols: ", paste(names(dt), collapse = ", "))
      print(dt)
    } else {
      dbg("[", label, "] peek fread failed")
    }
  }
  invisible(NULL)
}

dbg_stop <- function(...) {
  dbg("FATAL: ", paste0(...))
  stop(paste0(...), call. = FALSE)
}

# safer fread: returns data.frame() if missing/empty/unreadable
safe_fread_df <- function(path, label = "safe_fread") {
  if (is.null(path) || is.na(path) || path == "" || !file.exists(path)) {
    dbg("[", label, "] missing: ", path)
    return(data.frame())
  }
  sz <- file.info(path)$size
  if (is.na(sz) || sz == 0) {
    dbg("[", label, "] EMPTY (0 bytes): ", path)
    return(data.frame())
  }
  dt <- tryCatch(data.table::fread(path), error = function(e) {
    dbg("[", label, "] fread ERROR: ", conditionMessage(e))
    return(NULL)
  })
  if (is.null(dt)) return(data.frame())
  as.data.frame(dt)
}

norm_str <- function(x) tolower(trimws(as.character(x)))

norm_drug <- function(x) {
  x <- norm_str(x)
  x <- ifelse(x == "pd901", "pd-0325901", x)
  x
}

safe_token <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("[^A-Za-z0-9_\\-\\.]+", "_", x)
  x
}

warn_first_by <- function(df, group_cols, label) {
  if (is.null(df) || nrow(df) == 0) {
    dbg("[", label, "] warn_first_by: input empty")
    return(df)
  }

  counts <- df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(n = n(), .groups = "drop")

  dup_n <- sum(counts$n > 1)
  if (dup_n > 0) {
    warning(sprintf("[%s] Found %d duplicated key group(s); using first() for each group.", label, dup_n))
    if (isTRUE(DEBUG)) {
      dbg("[", label, "] duplicated groups (top):")
      print(utils::head(counts[counts$n > 1, ], DEBUG_HEAD_N))
    }
  }

  df %>%
    group_by(across(all_of(group_cols))) %>%
    slice(1) %>%
    ungroup()
}

extract_date_hour <- function(experiment_id) {
  pattern <- "(\\d{6})_?(\\d{2,3})?"
  m <- str_match(experiment_id, pattern)
  date <- m[,2]; hour <- m[,3]
  date[is.na(date)] <- NA
  hour[is.na(hour)] <- 48
  list(date = date, hour = hour)
}

canon_name <- function(x) {
  x <- norm_drug(x)
  x <- gsub("\\s+", "", x)
  x <- gsub("_", "", x)
  x <- gsub("-", "", x)
  x
}

canon_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("\\s+", "", x)
  x <- gsub("_", "", x)
  x <- gsub("-", "", x)
  x
}

split_plus_once <- function(x) {
  x <- as.character(x)
  if (!grepl("\\+", x)) return(c(NA_character_, NA_character_))
  parts <- strsplit(x, "\\+", fixed = FALSE)[[1]]
  if (length(parts) < 2) return(c(NA_character_, NA_character_))
  if (length(parts) > 2) {
    a <- parts[1]
    b <- paste(parts[-1], collapse = "+")
    return(c(a, b))
  }
  c(parts[1], parts[2])
}

parse_combo_key <- function(combo_key) {
  x <- as.character(combo_key)
  if (!grepl("_", x) || !grepl("\\+", x)) {
    return(list(base = NA_character_, drug1 = NA_character_, drug2 = NA_character_))
  }
  base <- sub("^([^_]+)_.*$", "\\1", x)
  combo_str <- sub("^[^_]+_", "", x)
  dp <- split_plus_once(combo_str)
  list(base = base, drug1 = dp[1], drug2 = dp[2])
}

escape_regex <- function(x) {
  x <- as.character(x)
  x <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
  x
}

# ============================================================
# fallback mapping (regex search in chem_name)
# ============================================================
fallback_id_from_chem_name <- function(pattern_raw, chem_tbl) {
  if (is.null(pattern_raw) || is.na(pattern_raw)) return(NA_character_)
  p <- trimws(as.character(pattern_raw))
  if (p == "" || tolower(p) %in% c("na","n/a","nan")) return(NA_character_)

  patt <- escape_regex(p)
  hits <- chem_tbl %>%
    filter(grepl(patt, chem_name, ignore.case = TRUE)) %>%
    pull(improve_drug_id)

  if (length(hits) > 0) return(as.character(hits[1]))
  NA_character_
}

fill_ids_with_fallback <- function(df, raw_col, id_col, chem_tbl, label) {
  stopifnot(raw_col %in% names(df), id_col %in% names(df))

  missing_idx <- which(is.na(df[[id_col]]) | df[[id_col]] == "")
  if (length(missing_idx) == 0) {
    dbg(label, ": no missing ", id_col, " to fill.")
    return(df)
  }

  dbg(label, ": attempting fallback fill for ", length(missing_idx),
      " missing ", id_col, " values using regex search in chem_name.")

  filled <- vapply(df[[raw_col]][missing_idx], fallback_id_from_chem_name, character(1), chem_tbl = chem_tbl)

  replace_idx <- missing_idx[!is.na(filled)]
  if (length(replace_idx) > 0) {
    df[[id_col]][replace_idx] <- filled[!is.na(filled)]
  }

  dbg(label, ": fallback filled ", length(replace_idx),
      " of ", length(missing_idx), " missing ", id_col, ".")
  df
}

# ============================================================
# treated token extraction (supports underscores)
# ============================================================
base_to_otherid_prefix <- function(common_name) {
  gsub("-", "_", trimws(as.character(common_name)))
}

extract_treated_token_from_other_id <- function(other_id, common_name) {
  oid <- trimws(as.character(other_id))
  cn  <- trimws(as.character(common_name))

  if (is.na(oid) || oid == "" || is.na(cn) || cn == "") {
    if (isTRUE(DEBUG)) dbg("[token_extract] invalid inputs oid/cn: oid=", oid, " cn=", cn)
    return(NA_character_)
  }

  treated_suffix_re <- "_(8hr|24hr)_treated_microtissue$"
  no_suffix <- sub(treated_suffix_re, "", oid)
  if (identical(no_suffix, oid)) {
    if (isTRUE(DEBUG)) dbg("[token_extract] suffix not found: ", oid)
    return(NA_character_)
  }

  base_prefixes <- unique(c(cn, gsub("-", "_", cn)))
  drug <- no_suffix
  stripped <- FALSE
  for (bp in base_prefixes) {
    drug2 <- sub(paste0("^", escape_regex(bp), "_"), "", drug)
    if (!identical(drug2, drug)) {
      drug <- drug2
      stripped <- TRUE
      break
    }
  }
  if (!stripped && isTRUE(DEBUG)) {
    dbg("[token_extract] base prefix NOT stripped (unexpected?) oid=", oid, " cn=", cn, " no_suffix=", no_suffix)
  }

  drug <- trimws(drug)
  if (drug == "" || tolower(drug) %in% c("na","n/a","nan")) {
    if (isTRUE(DEBUG)) dbg("[token_extract] empty/NA drug after parse: oid=", oid, " cn=", cn, " parsed=", drug)
    return(NA_character_)
  }
  drug
}

# ============================================================
# START
# ============================================================
dbg_hr("SCRIPT START")
dbg("SKIP_FAST=", SKIP_FAST)
dbg("samples=", samples)
dbg("drugfile=", drugfile)
dbg("out_prefix=", out_prefix)

synLogin(authToken = PAT)
dbg("Synapse login OK")

# -----------------------
# load inputs
# -----------------------
dbg_hr("LOADING INPUTS")
samples_all <- fread(samples) %>% as.data.frame()
dbg_df(samples_all, "samples_all")
if (!all(c("improve_sample_id","common_name","model_type") %in% names(samples_all))) {
  dbg_stop("samples.csv missing required columns: improve_sample_id, common_name, model_type")
}

samples_df <- samples_all %>%
  select(improve_sample_id, common_name, model_type) %>%
  distinct()
dbg_df(samples_df, "samples_df")

pdx_samps <- filter(samples_df, model_type == "patient derived xenograft")
mt_samps  <- filter(samples_df, model_type == "xenograft derived organoid")
dbg_df(pdx_samps, "pdx_samps")
dbg_df(mt_samps,  "mt_samps")

drug_map <- fread(drugfile) %>%
  select(improve_drug_id, chem_name) %>%
  distinct() %>%
  mutate(
    improve_drug_id = as.character(improve_drug_id),
    chem_name       = norm_drug(chem_name)
  )
dbg_df(drug_map, "drug_map")
dbg("drug_map unique chem_name=", length(unique(drug_map$chem_name)))

manifest <- synTableQuery("select * from syn53503360")$asDataFrame() %>%
  rename(common_name = Sample) %>%
  as.data.table()
dbg_df(as.data.frame(manifest), "manifest")
if (!all(c("common_name","MicroTissueDrugFolder") %in% names(manifest))) {
  dbg("manifest columns: ", paste(names(manifest), collapse=", "))
  dbg_stop("manifest missing MicroTissueDrugFolder and/or common_name")
}

# ============================================================
# OPTIONAL: MT / PDX / Combine
# ============================================================
if (!SKIP_FAST) {

  # ────────────────────────────────────────────────
  # MicroTissue Experiments
  # ────────────────────────────────────────────────
  dbg_hr("MT SECTION")

  getDrugDataByParent <- function(parid, sampleId) {
    dbg_hr(paste0("getDrugDataByParent(parid=", parid, ", sampleId=", sampleId, ")"))
    q <- sprintf(
      "select id,name,experimentalCondition,parentId from syn21993642 where parentId='%s'",
      parid
    )
    dbg("MT query: ", q)

    qtab <- tryCatch({
      synTableQuery(q)$asDataFrame()
    }, error = function(e) {
      dbg("[MT] synTableQuery ERROR for parentId=", parid, " : ", conditionMessage(e))
      return(NULL)
    })

    if (is.null(qtab) || nrow(qtab) == 0) {
      dbg("[MT] qtab empty for parentId=", parid)
      return(NULL)
    }

    dbg_df(qtab, "MT qtab raw")
    qtab <- qtab %>%
      filter(!is.na(experimentalCondition), name != "synapse_storage_manifest.csv") %>%
      select(id, name, experimentalCondition)

    dbg_df(qtab, "MT qtab filtered")
    if (nrow(qtab) == 0) {
      dbg("[MT] qtab filtered to 0 rows for parentId=", parid, " (all NA experimentalCondition?)")
      return(NULL)
    }

    # Loop over child files
    out_list <- vector("list", length = nrow(qtab))
    for (ii in seq_len(nrow(qtab))) {
      x <- qtab$id[ii]
      info <- qtab[ii, , drop = FALSE]
      d <- extract_date_hour(info$name)

      dbg_hr(paste0("[MT child ", ii, "/", nrow(qtab), "] synID=", x, " name=", info$name))
      sget <- tryCatch(synGet(x), error = function(e) {
        dbg("[MT] synGet ERROR for ", x, ": ", conditionMessage(e))
        return(NULL)
      })
      if (is.null(sget)) {
        out_list[[ii]] <- NULL
        next
      }

      pth <- sget$path
      dbg_file(pth, label = paste0("MT child file ", x))

      raw <- tryCatch({
        data.table::fread(pth)
      }, error = function(e) {
        dbg("[MT] fread ERROR for ", x, " path=", pth, " : ", conditionMessage(e))
        return(NULL)
      })
      if (is.null(raw) || nrow(raw) == 0) {
        dbg("[MT] raw empty for ", x)
        out_list[[ii]] <- NULL
        next
      }

      dbg("[MT] raw rows=", nrow(raw), " cols=", ncol(raw))
      dbg("[MT] raw cols: ", paste(names(raw), collapse=", "))
      if ("response_type" %in% names(raw)) {
        dbg("[MT] response_type top:")
        print(utils::head(sort(table(raw$response_type), decreasing = TRUE), DEBUG_HEAD_N))
      } else {
        dbg("[MT] WARNING: raw missing response_type column")
      }

      # Must have these columns for downstream
      needed <- c("response_type","dosage","response","compound_name")
      miss <- setdiff(needed, names(raw))
      if (length(miss) > 0) {
        dbg("[MT] MISSING required cols in ", x, ": ", paste(miss, collapse=", "))
        out_list[[ii]] <- NULL
        next
      }

      raw2 <- raw %>%
        filter(response_type == "percent viability")

      dbg("[MT] rows after filter(response_type=='percent viability') = ", nrow(raw2))
      if (nrow(raw2) == 0) {
        dbg("[MT] NOTE: no 'percent viability' rows in ", x, " (check response_type values)")
        out_list[[ii]] <- NULL
        next
      }

      # Transmute final
      tmp <- raw2 %>%
        transmute(
          improve_sample_id = sampleId,
          DOSE              = (10^dosage) * 1e6,
          GROWTH            = response,
          source            = "NF Data Portal",
          chem_name         = compound_name,
          study             = paste0("MT ", d$date, " exp"),
          time              = d$hour
        )

      dbg_df(tmp, label = paste0("MT tmp ", x))
      out_list[[ii]] <- tmp
    }

    out <- dplyr::bind_rows(out_list)
    dbg_df(out, label = paste0("MT OUT parentId=", parid))
    if (nrow(out) == 0) {
      dbg("[MT] parentId=", parid, " produced 0 rows total")
      return(NULL)
    }
    out
  }

  # Build MT folder map
  mts_map <- manifest %>%
    select(common_name, MicroTissueDrugFolder) %>%
    inner_join(mt_samps, by = "common_name") %>%
    separate_rows(MicroTissueDrugFolder, sep = ",") %>%
    filter(!is.na(MicroTissueDrugFolder), MicroTissueDrugFolder != "NA") %>%
    select(improve_sample_id, folder = MicroTissueDrugFolder)

  dbg_df(as.data.frame(mts_map), "mts_map")
  dbg("mts_map unique folders=", length(unique(mts_map$folder)))

  if (nrow(mts_map) == 0) {
    dbg_stop("MT: mts_map is empty. Likely manifest join failed or MicroTissueDrugFolder missing/NA.")
  }

  # Pull MT data for each folder
  mt_list <- vector("list", length = nrow(mts_map))
  for (i in seq_len(nrow(mts_map))) {
    sample_id <- mts_map$improve_sample_id[i]
    folder    <- mts_map$folder[i]

    if (isTRUE(DEBUG) && i <= DEBUG_MAX_ITEMS) {
      dbg_hr(paste0("MT MAP ROW ", i, "/", nrow(mts_map), " sample_id=", sample_id, " folder=", folder))
    }

    mt_list[[i]] <- getDrugDataByParent(folder, sample_id)
  }

  mt_data <- dplyr::bind_rows(mt_list)
  dbg_df(mt_data, "mt_data")

  # If mt_data empty, dump context + stop (this is the earliest real failure)
  if (nrow(mt_data) == 0) {
    dbg("MT FAILURE: mt_data is empty. Writing debug dumps.")
    dbg_out1 <- file.path("/tmp", paste0(out_prefix, "_DEBUG_mts_map.tsv"))
    dbg_out2 <- file.path("/tmp", paste0(out_prefix, "_DEBUG_manifest_mt_subset.tsv"))
    fwrite(as.data.table(mts_map), dbg_out1, sep = "\t")
    fwrite(as.data.table(manifest %>% as.data.frame() %>% filter(common_name %in% mt_samps$common_name)),
           dbg_out2, sep = "\t")
    dbg("Wrote: ", dbg_out1)
    dbg("Wrote: ", dbg_out2)
    dbg_stop("MT produced 0 rows. Investigate printed logs above + debug TSVs in /tmp.")
  }

  # Join to drug map
  mt_curve <- mt_data %>%
    mutate(chem_name = norm_drug(chem_name))

  # show join coverage BEFORE filtering NA improve_drug_id
  mt_joined <- mt_curve %>%
    left_join(drug_map, by = "chem_name")

  dbg_df(mt_joined, "mt_joined")
  n_unmapped <- sum(is.na(mt_joined$improve_drug_id))
  dbg("MT join unmapped chem_name rows=", n_unmapped, " / ", nrow(mt_joined))
  if (n_unmapped > 0) {
    top_unmapped <- mt_joined %>%
      filter(is.na(improve_drug_id)) %>%
      count(chem_name, sort = TRUE) %>%
      head(20)
    dbg("MT top unmapped chem_name (up to 20):")
    print(top_unmapped)
    dbg_unmapped_path <- file.path("/tmp", paste0(out_prefix, "_DEBUG_mt_unmapped_chem_name.tsv"))
    fwrite(as.data.table(mt_joined %>% filter(is.na(improve_drug_id))), dbg_unmapped_path, sep = "\t")
    dbg("Wrote: ", dbg_unmapped_path)
  }

  mt_curve <- mt_joined %>%
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

  dbg_df(mt_curve, "mt_curve (final)")

  out_mt_curve_path <- file.path("/tmp", paste0(out_prefix, "_mt_curve_data.tsv"))
  fwrite(mt_curve, out_mt_curve_path, sep = "\t")
  dbg_file(out_mt_curve_path, "MT curve output")
  message("Wrote MT curve data")

  if (nrow(mt_curve) == 0) {
    dbg_stop("MT curve is EMPTY after join/filter. Root cause is likely chem_name mismatch vs drugfile.tsv. See DEBUG_mt_unmapped_chem_name.tsv.")
  }

  # Fit curves
  out_mt_prefix <- paste0("/tmp/", out_prefix, "_mt_experiments")
  cmd <- sprintf("/opt/venv/bin/python fit_curve.py --input %s --output %s",
                 out_mt_curve_path, out_mt_prefix)
  dbg_hr("RUN fit_curve.py (MT)")
  dbg("CMD: ", cmd)
  rc <- system(cmd)
  dbg("fit_curve.py return code: ", rc)

  # fit_curve.py writes "<prefix>.0"
  mt_exp0 <- paste0(out_mt_prefix, ".0")
  dbg_file(mt_exp0, "MT experiments .0 (pre-rename)")
  if (!file.exists(mt_exp0) || file.info(mt_exp0)$size == 0) {
    dbg_stop("fit_curve.py did not produce a non-empty MT experiments file: ", mt_exp0,
             " (check fit_curve.py logs + MT curve input)")
  }

  mt_exp_path <- paste0("/tmp/", out_prefix, "_mt_experiments.tsv")
  ok_rename <- file.rename(mt_exp0, mt_exp_path)
  dbg("rename ", mt_exp0, " -> ", mt_exp_path, " ok=", ok_rename)
  dbg_file(mt_exp_path, "MT experiments TSV")
  message("Wrote MT experiments")

  # ────────────────────────────────────────────────
  # PDX Experiments
  # ────────────────────────────────────────────────
  dbg_hr("PDX SECTION")

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
  if (is.null(pdx_map) || nrow(pdx_map) == 0) {
    dbg("PDX: pdx_map empty (no PDX drug data)")
  } else {
    dbg_df(pdx_map, "pdx_map")
  }

  pdx_meta <- if (!is.null(pdx_map) && nrow(pdx_map) > 0) {
    do.call(rbind, lapply(seq_len(nrow(pdx_map)), function(i) {
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
  } else {
    data.frame()
  }

  all_pdx <- if (nrow(pdx_meta) > 0) {
    do.call(rbind, lapply(seq_len(nrow(pdx_meta)), function(i) {
      m   <- pdx_meta[i, ]
      pth <- synGet(m$child_id)$path
      raw <- if (grepl("\\.xlsx?$", pth)) read_xlsx(pth) else read_csv(pth)

      sec_opts  <- c("compound 2_name", "compound_2_name")
      drug2_col <- intersect(sec_opts, names(raw))[1]
      compound2 <- if (!is.na(drug2_col)) raw[[drug2_col]] else NA_character_

      df <- data.frame(
        child_id                     = m$child_id,
        specimen_id                  = raw$specimen_id,
        compound_name                = raw$compound_name,
        compound_2_name              = compound2,
        experimental_time_point      = raw$experimental_time_point,
        experimental_time_point_unit = raw$experimental_time_point_unit,
        assay_value                  = raw$assay_value,
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
  } else {
    data.frame()
  }

  pdx_data <- if (nrow(all_pdx) > 0) merge(all_pdx, pdx_meta, by="child_id") else data.frame()
  if (nrow(pdx_data) > 0) {
    pdx_data <- within(pdx_data, {
      experiment <- parentId
      model_id   <- improve_sample_id
    })

    has_ctl <- tapply(pdx_data$treatment == "control", pdx_data$experiment, any)
    pdx_data <- pdx_data[pdx_data$experiment %in% names(has_ctl)[has_ctl], ]

    pdx_data <- pdx_data[ , c("experiment","specimen_id","treatment",
                              "time","time_unit","volume","model_id")]

    pdx_data$treatment <- gsub("doxorubinsin", "doxorubicin", pdx_data$treatment, ignore.case = TRUE)
    pdx_data <- na.omit(pdx_data)

    out_pdx_curve <- file.path("/tmp", paste0(out_prefix, "_pdx_curve_data.tsv"))
    fwrite(pdx_data, out_pdx_curve, sep = "\t")
    dbg_file(out_pdx_curve, "PDX curve output")
    message("Wrote PDX curve data")

    cmd2 <- sprintf(
      "/opt/venv/bin/python calc_pdx_metrics.py %s --drugfile %s --outprefix %s --source 'NF Data Portal' --study 'MPNST PDX'",
      out_pdx_curve,
      drugfile,
      paste0("/tmp/", out_prefix, "_pdx")
    )
    dbg("CMD: ", cmd2)
    rc2 <- system(cmd2)
    dbg("calc_pdx_metrics.py return code: ", rc2)
    message("Wrote PDX experiments to /tmp/", out_prefix, "_pdx_experiments.tsv and combinations")
  } else {
    dbg("PDX: pdx_data empty; skipping PDX metrics")
  }

  # ────────────────────────────────────────────────
  # Combine all Experiments
  # ────────────────────────────────────────────────
  dbg_hr("COMBINE MT + PDX EXPERIMENTS")

  mt_exp_path2  <- paste0("/tmp/", out_prefix, "_mt_experiments.tsv")
  pdx_exp_path2 <- paste0("/tmp/", out_prefix, "_pdx_experiments.tsv")

  mt_exp <- safe_fread_df(mt_exp_path2, label = "read MT experiments")
  pdx_exp <- safe_fread_df(pdx_exp_path2, label = "read PDX experiments")

  dbg_df(mt_exp, "mt_exp (read)")
  dbg_df(pdx_exp, "pdx_exp (read)")

  # Make sure dose_response_value exists so mutate won't crash
  if (!("dose_response_value" %in% names(mt_exp)))  mt_exp$dose_response_value  <- NA
  if (!("dose_response_value" %in% names(pdx_exp))) pdx_exp$dose_response_value <- NA

  mt_exp  <- mt_exp  %>% mutate(dose_response_value = as.character(dose_response_value))
  pdx_exp <- pdx_exp %>% mutate(dose_response_value = as.character(dose_response_value))

  all_exp <- bind_rows(mt_exp, pdx_exp)
  dbg_df(all_exp, "all_exp (combined)")

  out_all <- paste0("/tmp/", out_prefix, "_experiments.tsv")
  fwrite(all_exp, out_all, sep = "\t")
  dbg_file(out_all, "combined experiments TSV")
  message("Wrote combined experiments: /tmp/", out_prefix, "_experiments.tsv")

  # rename combinations
  pdx_combo_src <- paste0("/tmp/", out_prefix, "_pdx_combinations.tsv")
  combo_dst     <- paste0("/tmp/", out_prefix, "_combinations.tsv")
  if (file.exists(pdx_combo_src)) {
    okc <- file.rename(pdx_combo_src, combo_dst)
    dbg("rename combos ", pdx_combo_src, " -> ", combo_dst, " ok=", okc)
  } else {
    dbg("PDX combinations file missing (ok if no PDX): ", pdx_combo_src)
  }

} else {
  dbg_hr("SKIP_FAST=TRUE")
  cat0("SKIP_FAST=TRUE: skipping MT + PDX + combine sections (debugging treated combos only).")
  cat0("Will still read: samples.csv, drugfile.tsv, and syn69801348.")
}

# ============================================================
# Treated-combination experiments from syn69801348 (auc ONLY - note we can change to fit_auc if needed)
# ============================================================
TREATED_COMBO_SYNID <- "syn69801348"
dbg_hr(paste0("TREATED COMBO SECTION from ", TREATED_COMBO_SYNID))

combo_path <- synGet(TREATED_COMBO_SYNID)$path
dbg_file(combo_path, "syn69801348 downloaded file")
combo_raw  <- fread(combo_path)

needed_cols <- c("source","improve_sample_id","improve_drug_id","study","time","time_unit","dose_response_metric","dose_response_value")
missing_cols <- setdiff(needed_cols, names(combo_raw))
if (length(missing_cols) > 0) {
  stop("syn69801348 is missing required columns: ", paste(missing_cols, collapse=", "), call. = FALSE)
}

dbg_df(as.data.frame(combo_raw), "combo_raw (pre-filter)")

combo_raw <- combo_raw %>%
  mutate(
    combo_key            = as.character(improve_sample_id),
    improve_sample_id    = as.character(improve_sample_id),
    improve_drug_id      = as.character(improve_drug_id),
    study                = as.character(study),
    time                 = suppressWarnings(as.integer(time)),
    time_unit            = as.character(time_unit),
    dose_response_metric = as.character(dose_response_metric),
    dose_response_value  = suppressWarnings(as.numeric(dose_response_value))
  ) %>%
  filter(tolower(dose_response_metric) == "auc") %>%    # This is auc only - can change to fit_auc if needed.
  filter(!is.na(time), !is.na(dose_response_value)) %>%
  distinct()

dbg("syn69801348 rows (auc): ", nrow(combo_raw))
dbg("combo_raw cols: ", paste(names(combo_raw), collapse=", "))
dbg("combo_raw head():")
print(head(combo_raw, 6))

# ---- parse combo_key into base_sample + drug1/drug2 from combo_key
if (nrow(combo_raw) > 0) {
  parsed <- lapply(combo_raw$combo_key, parse_combo_key)
  combo_raw$base_sample <- vapply(parsed, `[[`, character(1), "base")
  combo_raw$drug1_raw   <- vapply(parsed, `[[`, character(1), "drug1")
  combo_raw$drug2_raw   <- vapply(parsed, `[[`, character(1), "drug2")
} else {
  combo_raw$base_sample <- character()
  combo_raw$drug1_raw <- character()
  combo_raw$drug2_raw <- character()
}

parse_fail <- combo_raw %>% filter(is.na(base_sample) | is.na(drug1_raw) | is.na(drug2_raw))
if (nrow(parse_fail) > 0) {
  warning("Some combo_key values failed parsing (showing up to 10):")
  print(head(parse_fail$combo_key, 10))
  dbg <- file.path("/tmp", paste0(out_prefix, "_treated_parse_fail.tsv"))
  fwrite(as.data.table(parse_fail), dbg, sep = "\t")
  cat0("Wrote: ", dbg)
}
combo_raw <- combo_raw %>% filter(!is.na(base_sample), !is.na(drug1_raw), !is.na(drug2_raw))
dbg("rows after combo_key parse filter: ", nrow(combo_raw))

combo_raw <- combo_raw %>%
  mutate(
    drug1_canon  = canon_name(drug1_raw),
    drug2_canon  = canon_name(drug2_raw),
    drug_canon   = canon_name(improve_drug_id),
    is_combo_row = grepl("\\+", improve_drug_id)
  )

# ============================================================
# Treated samples mapping
# ============================================================
dbg_hr("TREATED SAMPLES MAPPING")

treated_samples_long <- samples_all %>%
  mutate(
    other_id          = as.character(other_id),
    model_type        = as.character(model_type),
    improve_sample_id = as.integer(improve_sample_id),
    common_name       = as.character(common_name)
  ) %>%
  filter(
    model_type == "xenograft derived organoid",
    !is.na(other_id),
    grepl("hr_treated_microtissue$", other_id)
  ) %>%
  mutate(
    base_common_name = common_name,
    base_prefix      = base_to_otherid_prefix(common_name),
    treated_token    = mapply(extract_treated_token_from_other_id, other_id, common_name),
    treated_token_u  = toupper(safe_token(treated_token))
  ) %>%
  select(
    base_common_name,
    treated_token,
    treated_token_u,
    treated_other_id = other_id,
    treated_improve_sample_id = improve_sample_id
  ) %>%
  distinct()

dbg("treated_samples_long rows: ", nrow(treated_samples_long))
dbg("treated_samples_long head():")
print(head(treated_samples_long, 10))

failed_token <- treated_samples_long %>% filter(is.na(treated_token) | treated_token == "")
if (nrow(failed_token) > 0) {
  warning("Some treated other_id rows failed treated_token extraction (showing up to 20):")
  print(head(failed_token, 20))
  dbg_fail_tok <- file.path("/tmp", paste0(out_prefix, "_treated_token_extract_fail.tsv"))
  fwrite(as.data.table(failed_token), dbg_fail_tok, sep = "\t")
  cat0("Wrote token extraction failures: ", dbg_fail_tok)
}

# ============================================================
# Build chem_tbl + mapping token -> improve_drug_id
# ============================================================
dbg_hr("TREATED TOKEN -> DRUG ID MAPPING")

chem_tbl <- drug_map %>%
  transmute(
    improve_drug_id = as.character(improve_drug_id),
    chem_name       = as.character(chem_name),
    chem_canon      = canon_name(chem_name),
    id_canon        = canon_id(improve_drug_id),
    abbrev3         = substr(toupper(gsub("[^A-Za-z]+","",chem_name)), 1, 3),
    chem_token_u    = toupper(safe_token(chem_name))
  ) %>%
  distinct()

dbg_df(chem_tbl, "chem_tbl")

treated_tokens <- treated_samples_long %>%
  distinct(treated_token, treated_token_u) %>%
  mutate(
    treated_token_canon = canon_name(treated_token),
    treated_token3 = substr(toupper(gsub("[^A-Za-z]+","",treated_token)), 1, 3)
  )

dbg_df(treated_tokens, "treated_tokens")

# These were updated on synapse so this isn't needed now.
manual_override <- c(
  # "DOX" = "doxorubicin",
  # "TNO" = "TNO155"
)

override_tbl <- data.frame(
  treated_token_u = character(0),
  override_val    = character(0),
  override_id     = character(0),
  stringsAsFactors = FALSE
)

if (length(manual_override) > 0) {
  override_tbl <- data.frame(
    treated_token_u = names(manual_override),
    override_val    = as.character(manual_override),
    stringsAsFactors = FALSE
  )
  resolve_override_to_id <- function(v) {
    if (is.na(v) || trimws(v) == "") return(NA_character_)
    v0 <- as.character(v)
    hit_id <- chem_tbl$improve_drug_id[chem_tbl$id_canon == canon_id(v0)]
    if (length(hit_id) > 0) return(hit_id[1])
    hit_nm <- chem_tbl$improve_drug_id[chem_tbl$chem_canon == canon_name(v0)]
    if (length(hit_nm) > 0) return(hit_nm[1])
    NA_character_
  }
  override_tbl$override_id <- vapply(override_tbl$override_val, resolve_override_to_id, character(1))
}

map_by_canon <- treated_tokens %>%
  left_join(chem_tbl %>% select(improve_drug_id, chem_canon),
            by = c("treated_token_canon" = "chem_canon")) %>%
  transmute(treated_token_u, canon_id = improve_drug_id)

map_by_token <- treated_tokens %>%
  left_join(chem_tbl %>% select(improve_drug_id, chem_token_u),
            by = c("treated_token_u" = "chem_token_u")) %>%
  transmute(treated_token_u, token_id = improve_drug_id)

map_by_abbrev3 <- treated_tokens %>%
  left_join(chem_tbl %>% select(improve_drug_id, abbrev3),
            by = c("treated_token3" = "abbrev3"), relationship = "many-to-many") %>%
  transmute(treated_token_u, abbr_id = improve_drug_id)

map_by_canon   <- warn_first_by(map_by_canon,   c("treated_token_u"), "treated_drug_map_canon")
map_by_token   <- warn_first_by(map_by_token,   c("treated_token_u"), "treated_drug_map_token")
map_by_abbrev3 <- warn_first_by(map_by_abbrev3, c("treated_token_u"), "treated_drug_map_abbrev3")

drug_token_map <- treated_tokens %>%
  left_join(override_tbl %>% select(treated_token_u, override_id), by = "treated_token_u") %>%
  left_join(map_by_canon,  by = "treated_token_u") %>%
  left_join(map_by_token,  by = "treated_token_u") %>%
  left_join(map_by_abbrev3, by = "treated_token_u") %>%
  mutate(mapped_improve_drug_id = dplyr::coalesce(override_id, canon_id, token_id, abbr_id)) %>%
  select(treated_token, treated_token_u, mapped_improve_drug_id)

drug_token_map <- warn_first_by(drug_token_map, c("treated_token_u"), "treated_drug_map_final")

mapping_out <- file.path("/tmp", paste0(out_prefix, "_treated_drug_mapping.tsv"))
fwrite(as.data.table(drug_token_map), mapping_out, sep = "\t")
dbg("Wrote treated drug mapping debug: ", mapping_out)
dbg_df(drug_token_map, "drug_token_map")

treated_samples_mapped <- treated_samples_long %>%
  left_join(drug_token_map %>% select(treated_token_u, mapped_improve_drug_id), by = "treated_token_u") %>%
  filter(!is.na(mapped_improve_drug_id)) %>%
  mutate(mapped_improve_drug_id = as.character(mapped_improve_drug_id))

dbg("treated_samples_mapped rows: ", nrow(treated_samples_mapped))
dbg_df(treated_samples_mapped, "treated_samples_mapped")

# ============================================================
# Build pairing tables
# ============================================================
dbg_hr("PAIRING auc SINGLE vs COMBO")

ctx_cols <- c("combo_key","base_sample","time","source","study","time_unit")

combo_ctx <- combo_raw %>%
  filter(is_combo_row) %>%
  transmute(
    combo_key, base_sample,
    drug1_raw, drug2_raw, drug1_canon, drug2_canon,
    time, source, study, time_unit,
    auc_combo = dose_response_value
  )

single_ctx <- combo_raw %>%
  filter(!is_combo_row) %>%
  transmute(
    combo_key, base_sample,
    time, source, study, time_unit,
    drugA_raw   = improve_drug_id,
    drugA_canon = drug_canon,
    auc_A       = dose_response_value
  )

dbg_df(combo_ctx,  "combo_ctx")
dbg_df(single_ctx, "single_ctx")

combo_ctx  <- warn_first_by(combo_ctx,  ctx_cols, "syn69801348_combo_auc")
single_ctx <- warn_first_by(single_ctx, c(ctx_cols,"drugA_canon"), "syn69801348_single_auc")

paired <- single_ctx %>%
  inner_join(combo_ctx, by = ctx_cols) %>%
  mutate(
    drugB_raw = dplyr::case_when(
      drugA_canon == drug1_canon ~ drug2_raw,
      drugA_canon == drug2_canon ~ drug1_raw,
      TRUE ~ NA_character_
    ),
    drugB_canon = canon_name(drugB_raw),
    delta_auc = auc_A - auc_combo
  )

dbg("paired rows: ", nrow(paired))
dbg_df(paired, "paired")

dbg_pairs <- file.path("/tmp", paste0(out_prefix, "_treated_pairs_debug.tsv"))
fwrite(as.data.table(paired), dbg_pairs, sep = "\t")
dbg("Wrote paired debug: ", dbg_pairs)

if (nrow(paired) == 0) {
  dbg_stop("paired is empty. Investigate /tmp/*_treated_pairs_debug.tsv and parse failures.")
}

# ============================================================
# Map canonical -> improve_drug_id (primary), then fallback regex
# ============================================================
dbg_hr("MAP drugA/drugB to improve_drug_id")

chem_to_id <- chem_tbl %>% distinct(chem_canon, improve_drug_id)
chem_to_id <- warn_first_by(chem_to_id, c("chem_canon"), "chem_to_id")

paired <- paired %>%
  left_join(chem_to_id %>% rename(drugA_id = improve_drug_id), by = c("drugA_canon" = "chem_canon")) %>%
  left_join(chem_to_id %>% rename(drugB_id = improve_drug_id), by = c("drugB_canon" = "chem_canon"))

missingA_before <- sum(is.na(paired$drugA_id))
missingB_before <- sum(is.na(paired$drugB_id))
dbg("Missing drugA_id before fallback: ", missingA_before, " / ", nrow(paired))
dbg("Missing drugB_id before fallback: ", missingB_before, " / ", nrow(paired))

paired <- fill_ids_with_fallback(paired, raw_col = "drugA_raw", id_col = "drugA_id", chem_tbl = chem_tbl, label = "drugA_id")
paired <- fill_ids_with_fallback(paired, raw_col = "drugB_raw", id_col = "drugB_id", chem_tbl = chem_tbl, label = "drugB_id")

missingA_after <- sum(is.na(paired$drugA_id))
missingB_after <- sum(is.na(paired$drugB_id))
dbg("Missing drugA_id after fallback: ", missingA_after, " / ", nrow(paired))
dbg("Missing drugB_id after fallback: ", missingB_after, " / ", nrow(paired))

fallback_dbg <- paired %>%
  transmute(
    combo_key, base_sample, time, source, study, time_unit,
    drugA_raw, drugA_canon, drugA_id,
    drugB_raw, drugB_canon, drugB_id,
    delta_auc,
    drugA_id_missing = is.na(drugA_id),
    drugB_id_missing = is.na(drugB_id)
  )

fallback_dbg_path <- file.path("/tmp", paste0(out_prefix, "_treated_drug_id_fallback_debug.tsv"))
fwrite(as.data.table(fallback_dbg), fallback_dbg_path, sep = "\t")
dbg("Wrote drug-id fallback debug: ", fallback_dbg_path)

dbg("paired with IDs mapped head():")
print(head(paired %>% select(combo_key, base_sample, time, drugA_raw, drugA_id, drugB_raw, drugB_id, delta_auc), 10))

# ============================================================
# Join to treated samples
# ============================================================
dbg_hr("JOIN paired -> treated samples")

treated_joined <- paired %>%
  filter(!is.na(drugA_id), !is.na(drugB_id)) %>%
  inner_join(
    treated_samples_mapped %>%
      transmute(
        base_sample = base_common_name,
        drugA_id = mapped_improve_drug_id,
        treated_improve_sample_id,
        treated_other_id,
        treated_token
      ),
    by = c("base_sample","drugA_id")
  )

dbg("treated_joined rows: ", nrow(treated_joined))
dbg_df(treated_joined, "treated_joined")

dbg_joined <- file.path("/tmp", paste0(out_prefix, "_treated_joined_debug.tsv"))
fwrite(as.data.table(treated_joined), dbg_joined, sep = "\t")
dbg("Wrote treated_joined debug: ", dbg_joined)

if (nrow(treated_joined) == 0) {
  dbg("EARLY EXIT: treated_joined is empty. Likely token->drugA_id mapping mismatch.")
  dbg("Example base_sample from paired: ", paste(head(unique(paired$base_sample), 10), collapse=", "))
  dbg("Example base_common_name from treated_samples_mapped: ", paste(head(unique(treated_samples_mapped$base_common_name), 10), collapse=", "))
  dbg("Example drugA_id from paired: ", paste(head(unique(paired$drugA_id), 10), collapse=", "))
  dbg("Example mapped_improve_drug_id: ", paste(head(unique(treated_samples_mapped$mapped_improve_drug_id), 10), collapse=", "))
  quit(save="no", status=0)
}

# ============================================================
# Outputs
# ============================================================
dbg_hr("WRITE TREATED OUTPUTS")

treated_out_ids <- treated_joined %>%
  transmute(
    source               = source,
    improve_sample_id    = as.integer(treated_improve_sample_id),
    improve_drug_id      = as.character(drugB_id),
    study                = study,
    time                 = as.integer(time),
    time_unit            = time_unit,
    dose_response_metric = "delta_auc",
    dose_response_value  = as.numeric(delta_auc)
  ) %>%
  distinct()

treated_ids_path <- file.path("/tmp", paste0(out_prefix, "_experiment_combinations_treated.tsv"))
fwrite(as.data.table(treated_out_ids), treated_ids_path, sep = "\t")
dbg("Wrote treated combination experiments (IDs): ", treated_ids_path, " [rows=", nrow(treated_out_ids), "]")
dbg_file(treated_ids_path, "treated IDs TSV")

treated_out_names <- treated_joined %>%
  transmute(
    source               = source,
    improve_sample_id    = as.character(combo_key),
    improve_drug_id      = as.character(drugB_raw),
    study                = study,
    time                 = as.integer(time),
    time_unit            = time_unit,
    dose_response_metric = "delta_auc",
    dose_response_value  = as.numeric(delta_auc)
  ) %>%
  distinct()

treated_names_path <- file.path("/tmp", paste0(out_prefix, "_experiment_combinations_treated_names.tsv"))
fwrite(as.data.table(treated_out_names), treated_names_path, sep = "\t")
dbg("Wrote treated combination experiments (NAMES): ", treated_names_path, " [rows=", nrow(treated_out_names), "]")
dbg_file(treated_names_path, "treated NAMES TSV")

# ============================================================
# FINAL: append treated-combo experiments into main experiments
# ============================================================
dbg_hr("FINAL APPEND treated -> experiments")

exp_path <- file.path("/tmp", paste0(out_prefix, "_experiments.tsv"))
trt_path <- treated_ids_path

if (!file.exists(exp_path)) stop("Missing experiments file: ", exp_path, call. = FALSE)
if (!file.exists(trt_path)) stop("Missing treated-combos file: ", trt_path, call. = FALSE)

exp_dt <- data.table::fread(exp_path)
trt_dt <- data.table::fread(trt_path)

dbg("exp_dt rows=", nrow(exp_dt), " cols=", ncol(exp_dt))
dbg("trt_dt rows=", nrow(trt_dt), " cols=", ncol(trt_dt))

missing_cols <- setdiff(names(exp_dt), names(trt_dt))
for (cc in missing_cols) trt_dt[[cc]] <- NA

extra_cols <- setdiff(names(trt_dt), names(exp_dt))
if (length(extra_cols) > 0) trt_dt <- trt_dt[, setdiff(names(trt_dt), extra_cols), with = FALSE]

data.table::setcolorder(trt_dt, names(exp_dt))

for (cc in intersect(c("improve_sample_id", "time"), names(exp_dt))) {
  exp_dt[[cc]] <- suppressWarnings(as.integer(exp_dt[[cc]]))
  trt_dt[[cc]] <- suppressWarnings(as.integer(trt_dt[[cc]]))
}
for (cc in intersect(c("dose_response_value"), names(exp_dt))) {
  exp_dt[[cc]] <- suppressWarnings(as.numeric(exp_dt[[cc]]))
  trt_dt[[cc]] <- suppressWarnings(as.numeric(trt_dt[[cc]]))
}

before_n <- nrow(exp_dt)
add_n    <- nrow(trt_dt)

final_dt <- data.table::rbindlist(list(exp_dt, trt_dt), use.names = TRUE, fill = TRUE)
final_dt <- unique(final_dt)

bak_path <- paste0(exp_path, ".bak")
file.copy(exp_path, bak_path, overwrite = TRUE)
data.table::fwrite(final_dt, exp_path, sep = "\t")

message(sprintf(
  "Appended treated combos into %s: +%d rows (was %d, now %d). Backup: %s",
  exp_path, add_n, before_n, nrow(final_dt), bak_path
))

dbg_hr("SCRIPT END")
