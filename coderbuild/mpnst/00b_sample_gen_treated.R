# 00_sample_gen_treated.R
# This script generates treated microtissue samples and appends them to an existing samples file.

suppressPackageStartupMessages({
  library(data.table)
  library(synapser)
  library(readxl)
})

# -----------------------
# logging helpers
# -----------------------
.ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
.p  <- function(...) cat(.ts(), paste0(...), "\n")

print_df_head <- function(df, n = 8) {
  if (is.null(df) || nrow(df) == 0) {
    .p("[PRINT] <0 rows>")
    return(invisible(NULL))
  }
  n <- min(n, nrow(df))
  print(utils::head(df, n))
  invisible(NULL)
}

# -----------------------
# constants
# -----------------------
CONST_CANCER  <- "Malignant peripheral nerve sheath tumor"
CONST_SPECIES <- "Homo sapiens (Human)"
CONST_MODEL   <- "xenograft derived organoid"
CONST_SOURCE  <- "NF Data Portal"

# NOTE: updated to correctly preserve hyphens in IDs (MN-2, JH-2-002, WU-225)
safe_token <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("[^[:alnum:]_.-]+", "_", x)
  x
}

# -----------------------
# Canonicalize individual IDs to hyphen style (MN-2, JH-2-002, WU-225)
# -----------------------
canonicalize_individual_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[_\\s]+", "-", x)   # underscores/spaces -> hyphen
  x <- gsub("-+", "-", x)        # collapse repeated hyphens

  # Normalize common patterns (only changes if it matches)
  x <- sub("^MN-?([0-9A-Za-z]+)$", "MN-\\1", x)
  x <- sub("^WU-?([0-9A-Za-z]+)$", "WU-\\1", x)
  x <- sub("^JH-?([0-9]+)-?([0-9A-Za-z]+)$", "JH-\\1-\\2", x)

  x
}

# -----------------------
# Migrate existing treated other_id prefixes from underscore->hyphen
# This is because some samples were mis-labeled.
# -----------------------
migrate_treated_other_id_prefix <- function(other_id) {
  x <- as.character(other_id)
  is_treated <- grepl("_treated_microtissue$", x)

  y <- x
  y[is_treated] <- sub("^MN_([0-9A-Za-z]+)_", "MN-\\1_", y[is_treated])
  y[is_treated] <- sub("^WU_([0-9A-Za-z]+)_", "WU-\\1_", y[is_treated])
  y[is_treated] <- sub("^JH_([0-9]+)_([0-9A-Za-z]+)_", "JH-\\1-\\2_", y[is_treated])

  y
}

# -----------------------
#  Timepoint normalizer (handles "24 HR", "8HR ", etc)
# -----------------------
normalize_timepoint <- function(tp) {
  tp <- trimws(as.character(tp))
  tpu <- toupper(tp)
  tpu <- gsub("\\s+", "", tpu)

  out <- tp
  out[tpu %in% c("8HR", "8H")]    <- "8h"
  out[tpu %in% c("24HR", "24H")]  <- "24h"
  out[tpu %in% c("0HR", "0H")]    <- "0h"

  is_h <- grepl("^[0-9]+h$", tolower(out))
  out[is_h] <- tolower(out[is_h])
  out
}

hours_from_timepoint <- function(tp_norm) {
  tp_norm <- tolower(trimws(as.character(tp_norm)))
  hrs <- gsub("h$", "", tp_norm)
  suppressWarnings(as.integer(hrs))
}

# -----------------------
# Drop controls: DMSO + untreated (and close variants)
# -----------------------
is_dmso <- function(drug) {
  d <- tolower(trimws(as.character(drug)))
  grepl("\\bdmso\\b", d) || grepl("dimethyl", d) || grepl("sulfoxide", d)
}

is_untreated <- function(drug) {
  d <- tolower(trimws(as.character(drug)))
  grepl("\\buntreated\\b", d) || grepl("no[_\\s-]*treat", d) || grepl("\\bvehicle\\b", d)
}

is_control <- function(drug) {
  isTRUE(is_dmso(drug)) || isTRUE(is_untreated(drug))
}

# ============================================================
# Builder function:
#   - Drops DMSO + untreated rows
#   - Aggregates replicates: ONE row per treated sample
#   - other_id  = {individual_id}_{drug}_{hours}hr_treated_microtissue
#   - other_names = comma-separated replicate specimen IDs
#
# KEY FIX:
#   - After first aggregation, compute other_id_val then
#     RE-aggregate by other_id_val to merge replicate IDs when
#     sanitization causes collisions (e.g. RMC-4630 vs RMC 4630).
# ============================================================
build_from_standard_map <- function(std_map, orig_samples, batch_key, next_id, require_base = TRUE) {
  std_map <- as.data.table(std_map)

  # Canonicalize IDs so other_id uses MN-2 / JH-2-002 / WU-225 formatting
  std_map[, individual_id := canonicalize_individual_id(individual_id)]

  # base_exists check (optional)
  std_map[, base_exists := individual_id %in% orig_samples$common_name]

  # drop controls early (DMSO + untreated)
  std_map[, drop_control := vapply(drug, is_control, logical(1))]
  n_ctrl <- sum(std_map$drop_control %in% TRUE)
  if (n_ctrl > 0) .p("[", batch_key, "] dropping control rows (DMSO/untreated): ", n_ctrl)
  std_map <- std_map[drop_control == FALSE | is.na(drop_control)]

  debug_tbl <- std_map[, .(
    batch        = batch_key,
    specimen_id,
    individual_id,
    drug,
    timepoint,
    base_exists
  )]

  keep_tbl <- if (require_base) std_map[base_exists == TRUE] else copy(std_map)

  .p("[", batch_key, "] rows after drop controls: ", nrow(std_map))
  .p("[", batch_key, "] require_base=", require_base)
  .p("[", batch_key, "] rows kept pre-aggregate: ", nrow(keep_tbl))
  .p("[", batch_key, "] preview:")
  print_df_head(keep_tbl, 8)

  if (nrow(keep_tbl) == 0) {
    return(list(
      new_samples_tbl = data.table(),
      map_tbl         = data.table(),
      debug_tbl       = debug_tbl,
      next_id         = next_id
    ))
  }

  # ---- First aggregation: collapse exact (individual_id, drug, timepoint) rows
  agg1 <- keep_tbl[
    ,
    .(specimen_ids = sort(unique(trimws(as.character(specimen_id))))),
    by = .(individual_id, drug, timepoint)
  ]

  # Remove NA/blank specimen IDs inside groups
  agg1[, specimen_ids := lapply(specimen_ids, function(v) v[!is.na(v) & nzchar(v)])]
  agg1 <- agg1[lengths(specimen_ids) > 0]

  agg1[, other_names := vapply(specimen_ids, function(v) paste(v, collapse = ", "), character(1))]
  agg1[, hours := hours_from_timepoint(timepoint)]
  agg1[is.na(hours), hours := -1L]

  # ---- Compute other_id_val per group (this is what ends up in samples.csv)
  agg1[, other_id_val := paste0(
    safe_token(individual_id), "_",
    safe_token(drug), "_",
    as.integer(hours), "hr_treated_microtissue"
  )]

  # ---- SECOND aggregation: merge replicate IDs when other_id_val collides
  agg <- agg1[
    ,
    .(
      individual_id = individual_id[1],
      drug          = drug[1],
      timepoint     = timepoint[1],
      hours         = hours[1],
      other_names   = paste(
        sort(unique(trimws(unlist(strsplit(paste(other_names, collapse = ", "), ",\\s*"))))),
        collapse = ", "
      )
    ),
    by = .(other_id_val)
  ]

  .p("[", batch_key, "] aggregate groups after merge-by-other_id: ", nrow(agg))
  .p("[", batch_key, "] aggregate head:")
  print_df_head(agg[, .(other_id_val, individual_id, drug, timepoint, other_names, hours)], 8)

  # fast lookup of existing treated keys in orig_samples
  orig_treated_keys <- paste0(as.character(orig_samples$other_id), "||", as.character(orig_samples$model_type))

  map_tbl <- data.table(
    batch                     = character(),
    individual_id             = character(),
    drug                      = character(),
    timepoint                 = character(),
    hours                     = integer(),
    other_id                  = character(),
    other_names               = character(),
    treated_improve_sample_id = integer()
  )

  new_samples_tbl <- data.table(
    other_id           = character(),
    common_name        = character(),
    other_id_source    = character(),
    other_names        = character(),
    cancer_type        = character(),
    species            = character(),
    model_type         = character(),
    improve_sample_id  = integer()
  )

  for (i in seq_len(nrow(agg))) {
    indiv <- as.character(agg$individual_id[i])
    drug  <- as.character(agg$drug[i])
    tp    <- as.character(agg$timepoint[i])
    hrs   <- as.integer(agg$hours[i])
    repnames <- as.character(agg$other_names[i])

    other_id_val <- as.character(agg$other_id_val[i])
    key <- paste0(other_id_val, "||", CONST_MODEL)

    # Skip if already exists in samples.csv
    if (key %in% orig_treated_keys) {
      .p("[", batch_key, "][SKIP] already exists other_id=", other_id_val)
      next
    }

    treated_id <- next_id
    next_id <- next_id + 1L

    new_samples_tbl <- rbind(
      new_samples_tbl,
      data.table(
        other_id          = other_id_val,
        common_name       = indiv,
        other_id_source   = CONST_SOURCE,
        other_names       = repnames,
        cancer_type       = CONST_CANCER,
        species           = CONST_SPECIES,
        model_type        = CONST_MODEL,
        improve_sample_id = treated_id
      ),
      fill = TRUE
    )

    map_tbl <- rbind(
      map_tbl,
      data.table(
        batch                     = batch_key,
        individual_id             = indiv,
        drug                      = drug,
        timepoint                 = tp,
        hours                     = hrs,
        other_id                  = other_id_val,
        other_names               = repnames,
        treated_improve_sample_id = treated_id
      ),
      fill = TRUE
    )

    .p("[", batch_key, "][ADD] ", other_id_val, " | reps={", repnames, "} -> ", treated_id)
  }

  list(
    new_samples_tbl = new_samples_tbl,
    map_tbl         = map_tbl,
    debug_tbl       = debug_tbl,
    next_id         = next_id
  )
}

write_batch_outputs <- function(batch_key, debug_tbl, map_tbl) {
  out_debug_path <- file.path("/tmp", sprintf("mpnst_samples_treated_%s_debug.csv", batch_key))
  out_map_path   <- file.path("/tmp", sprintf("mpnst_samples_treated_%s_map.csv", batch_key))

  .p("[WRITE] ", batch_key, " debug -> ", out_debug_path)
  fwrite(debug_tbl, out_debug_path)
  .p("[WRITE] OK")

  .p("[WRITE] ", batch_key, " map -> ", out_map_path)
  fwrite(map_tbl, out_map_path)
  .p("[WRITE] OK")
}

# ============================================================
# Batch 1
# ============================================================
BATCH1_KEY <- "batch1"
BATCH1_SAMPLEMAP_SYNID <- "syn64608368"

fetch_batch1_map <- function() {
  .p("[BATCH1] synGet(", BATCH1_SAMPLEMAP_SYNID, ") ...")
  pth <- synGet(BATCH1_SAMPLEMAP_SYNID)$path
  dt <- as.data.table(fread(pth))
  .p("[BATCH1] map dim: ", nrow(dt), " x ", ncol(dt))
  print_df_head(dt, 6)
  dt
}

normalize_batch1_map <- function(raw_map) {
  dt <- copy(as.data.table(raw_map))
  dt[, specimenID   := trimws(as.character(specimenID))]
  dt[, individualID := trimws(as.character(individualID))]
  dt[, Drug         := trimws(as.character(Drug))]
  dt[, Timepoint    := trimws(as.character(Timepoint))]

  dt <- dt[!is.na(specimenID) & nzchar(specimenID) & !is.na(individualID) & nzchar(individualID)]

  dt[, .(
    specimen_id   = specimenID,
    individual_id = individualID,
    drug          = Drug,
    timepoint     = normalize_timepoint(Timepoint)
  )]
}

run_batch1 <- function(orig_samples, next_id) {
  build_from_standard_map(normalize_batch1_map(fetch_batch1_map()), orig_samples, BATCH1_KEY, next_id, require_base = TRUE)
}

# ============================================================
# Batch 2
# ============================================================
BATCH2_KEY <- "batch2"
BATCH2_SAMPLEMAP_SYNID <- "syn66302373"

fetch_batch2_map <- function() {
  .p("[BATCH2] synGet(", BATCH2_SAMPLEMAP_SYNID, ") ...")
  pth <- synGet(BATCH2_SAMPLEMAP_SYNID)$path
  dt <- as.data.table(fread(pth))
  .p("[BATCH2] map dim: ", nrow(dt), " x ", ncol(dt))
  print_df_head(dt, 6)
  dt
}

normalize_batch2_map <- function(raw_map) {
  dt <- copy(as.data.table(raw_map))
  dt[, specimenID   := trimws(as.character(specimenID))]
  dt[, individualID := trimws(as.character(individualID))]
  dt[, Drug         := trimws(as.character(Drug))]
  dt[, Timepoint    := trimws(as.character(Timepoint))]

  # MN2 -> MN-2
  dt[individualID == "MN2", individualID := "MN-2"]

  dt <- dt[!is.na(specimenID) & nzchar(specimenID) & !is.na(individualID) & nzchar(individualID)]

  dt[, .(
    specimen_id   = specimenID,
    individual_id = individualID,
    drug          = Drug,
    timepoint     = normalize_timepoint(Timepoint)
  )]
}

run_batch2 <- function(orig_samples, next_id) {
  build_from_standard_map(normalize_batch2_map(fetch_batch2_map()), orig_samples, BATCH2_KEY, next_id, require_base = TRUE)
}

# ============================================================
# Batch 3
# ============================================================
BATCH3_KEY <- "batch3"
BATCH3_SAMPLEMAP_SYNID <- "syn66050299"

fetch_batch3_map <- function() {
  .p("[BATCH3] synGet(", BATCH3_SAMPLEMAP_SYNID, ") ...")
  pth <- synGet(BATCH3_SAMPLEMAP_SYNID)$path
  dt <- as.data.table(fread(pth))
  .p("[BATCH3] map dim: ", nrow(dt), " x ", ncol(dt))
  print_df_head(dt, 6)
  dt
}

normalize_batch3_map <- function(raw_map) {
  dt <- copy(as.data.table(raw_map))
  dt[, specimenID   := trimws(as.character(specimenID))]
  dt[, individualID := trimws(as.character(individualID))]
  dt[, Drug         := trimws(as.character(Drug))]
  dt[, Timepoint    := trimws(as.character(Timepoint))]

  dt <- dt[!is.na(specimenID) & nzchar(specimenID) & !is.na(individualID) & nzchar(individualID)]

  dt[, .(
    specimen_id   = specimenID,
    individual_id = individualID,
    drug          = Drug,
    timepoint     = normalize_timepoint(Timepoint)
  )]
}

run_batch3 <- function(orig_samples, next_id) {
  build_from_standard_map(normalize_batch3_map(fetch_batch3_map()), orig_samples, BATCH3_KEY, next_id, require_base = TRUE)
}

# ============================================================
# Batch 4 (Excel)
# ============================================================
BATCH4_KEY <- "batch4"
BATCH4_ANNOT_SYNID <- "syn68884278"
BATCH4_SHEET_NAME <- "Samplelist"

find_col <- function(nms, candidates) {
  nms_clean  <- tolower(gsub("\\s+", " ", trimws(nms)))
  cand_clean <- tolower(gsub("\\s+", " ", trimws(candidates)))
  idx <- match(cand_clean, nms_clean, nomatch = 0)
  idx <- idx[idx > 0]
  if (length(idx) > 0) return(nms[idx[1]])
  NULL
}

fetch_batch4_map <- function() {
  .p("[BATCH4] synGet(", BATCH4_ANNOT_SYNID, ") ...")
  pth <- synGet(BATCH4_ANNOT_SYNID)$path
  .p("[BATCH4] read_excel(sheet='", BATCH4_SHEET_NAME, "') from: ", pth)
  dt <- as.data.table(readxl::read_excel(pth, sheet = BATCH4_SHEET_NAME))
  .p("[BATCH4] map dim: ", nrow(dt), " x ", ncol(dt))
  .p("[BATCH4] map colnames: ", paste(names(dt), collapse = " | "))
  .p("[BATCH4] map head:")
  print_df_head(dt, 6)
  dt
}

normalize_batch4_map <- function(raw_map) {
  dt <- copy(as.data.table(raw_map))

  c_ind      <- find_col(names(dt), c("Individual ID"))
  c_drug     <- find_col(names(dt), c("Drug"))
  c_tp       <- find_col(names(dt), c("Timepoint"))
  c_labelcap <- find_col(names(dt), c("Label on Cap"))
  c_labelvial<- find_col(names(dt), c("Label on Vial"))

  .p("[BATCH4] resolved columns:")
  .p("  individual_id=", ifelse(is.null(c_ind), "<NOT FOUND>", c_ind))
  .p("  drug         =", ifelse(is.null(c_drug), "<NOT FOUND>", c_drug))
  .p("  timepoint    =", ifelse(is.null(c_tp), "<NOT FOUND>", c_tp))
  .p("  label_on_cap =", ifelse(is.null(c_labelcap), "<NOT FOUND>", c_labelcap))
  .p("  label_on_vial=", ifelse(is.null(c_labelvial), "<NOT FOUND>", c_labelvial))

  if (is.null(c_ind) || is.null(c_drug) || is.null(c_tp)) {
    .p("[BATCH4][WARN] Missing required columns; returning 0 rows.")
    return(data.table(specimen_id=character(), individual_id=character(), drug=character(), timepoint=character()))
  }

  if (!is.null(c_labelcap)) {
    dt[, specimen_id := trimws(as.character(get(c_labelcap)))]
  } else if (!is.null(c_labelvial)) {
    dt[, specimen_id := trimws(as.character(get(c_labelvial)))]
  } else {
    dt[, specimen_id := ""]
  }

  dt[, individual_id := trimws(as.character(get(c_ind)))]
  dt[, drug          := trimws(as.character(get(c_drug)))]
  dt[, timepoint     := normalize_timepoint(get(c_tp))]

  dt[is.na(specimen_id) | !nzchar(specimen_id), specimen_id := "B4_unknown"]

  out <- dt[
    !is.na(individual_id) & nzchar(individual_id) &
      !is.na(drug) & nzchar(drug) &
      !is.na(timepoint) & nzchar(as.character(timepoint)),
    .(specimen_id, individual_id, drug, timepoint)
  ]

  .p("[BATCH4] normalized dim: ", nrow(out), " x ", ncol(out))
  .p("[BATCH4] normalized head:")
  print_df_head(out, 8)

  out
}

run_batch4 <- function(orig_samples, next_id) {
  build_from_standard_map(
    normalize_batch4_map(fetch_batch4_map()),
    orig_samples,
    BATCH4_KEY,
    next_id,
    require_base = FALSE
  )
}

# ============================================================
# Register all batches (order)
# ============================================================
ALL_BATCH_RUNNERS <- list(
  list(key=BATCH1_KEY, fn=run_batch1),
  list(key=BATCH2_KEY, fn=run_batch2),
  list(key=BATCH3_KEY, fn=run_batch3),
  list(key=BATCH4_KEY, fn=run_batch4)
)

# ============================================================
# MAIN
# ============================================================
args <- commandArgs(trailingOnly = TRUE)

.p("============================================================")
.p("=== 00b_sample_gen_treated.R (RUN ALL BATCHES) ===")
.p("Timestamp: ", .ts())
.p("Args: ", paste(args, collapse = " | "))
.p("============================================================")

if (length(args) < 1) {
  stop(
    "Usage: Rscript 00b_sample_gen_treated.R <prev_samples_path>\n",
    "Example: Rscript 00b_sample_gen_treated.R /tmp/mpnst_samples.csv",
    call. = FALSE
  )
}

prev_samples_path <- args[1]
.p("[ARGS] prev_samples_path=", prev_samples_path)
if (!file.exists(prev_samples_path)) stop("File not found: ", prev_samples_path, call. = FALSE)

.p("[READ] fread(samples) ...")
orig_samples <- fread(prev_samples_path)
.p("[READ] orig_samples dim: ", nrow(orig_samples), " x ", ncol(orig_samples))
print_df_head(orig_samples, 6)

# ------------------------------------------------------------
# Migrate already-existing treated other_ids to hyphen style
# so future runs match your desired format and don't skip forever
# ------------------------------------------------------------
orig_samples[, other_id_old := other_id]
orig_samples[, other_id := migrate_treated_other_id_prefix(other_id)]

n_changed <- sum(orig_samples$other_id != orig_samples$other_id_old, na.rm = TRUE)
if (n_changed > 0) .p("[MIGRATE] updated treated other_id prefix for ", n_changed, " row(s)")
orig_samples[, other_id_old := NULL]

# If migration created duplicates (same other_id + model_type), merge other_names and drop dup rows
orig_samples[, dedup_key := paste(other_id, model_type, sep = "||")]
dup_keys <- orig_samples[duplicated(dedup_key), unique(dedup_key)]

if (length(dup_keys) > 0) {
  .p("[MIGRATE] duplicates after migration: ", length(dup_keys), " key(s). Merging other_names + dropping dup rows.")

  merged <- orig_samples[dedup_key %in% dup_keys, .(
    merged_other_names = paste(
      sort(unique(trimws(unlist(strsplit(paste(na.omit(other_names), collapse = ", "), ",\\s*"))))),
      collapse = ", "
    )
  ), by = dedup_key]

  orig_samples[merged, on = "dedup_key", other_names := i.merged_other_names]
  orig_samples <- orig_samples[!duplicated(dedup_key)]
}

orig_samples[, dedup_key := NULL]

token <- Sys.getenv("SYNAPSE_AUTH_TOKEN")
.p("[ENV] SYNAPSE_AUTH_TOKEN present? ", nzchar(token))
if (!nzchar(token)) stop("Please set SYNAPSE_AUTH_TOKEN in your environment.", call. = FALSE)

.p("[SYN] synLogin() ...")
synLogin(authToken = token)
.p("[SYN] synLogin() OK")

max_id <- suppressWarnings(max(as.integer(orig_samples$improve_sample_id), na.rm = TRUE))
if (!is.finite(max_id)) max_id <- 0L
next_id <- max_id + 1L
.p("[IDS] max improve_sample_id=", max_id, " next_id=", next_id)

.p("============================================================")
.p("[RUN] Running ", length(ALL_BATCH_RUNNERS), " batch runner(s) in order...")
.p("============================================================")

for (item in ALL_BATCH_RUNNERS) {
  batch_key <- item$key
  runner    <- item$fn

  res <- runner(orig_samples = orig_samples, next_id = next_id)

  if (!is.null(res$new_samples_tbl) && nrow(res$new_samples_tbl) > 0) {
    orig_samples <- rbindlist(list(orig_samples, res$new_samples_tbl), use.names = TRUE, fill = TRUE)

    orig_samples[, dedup_key := paste(other_id, model_type, sep = "||")]
    before <- nrow(orig_samples)
    orig_samples <- orig_samples[!duplicated(dedup_key)]
    after <- nrow(orig_samples)
    if (after < before) .p("[", batch_key, "][WARN] post-append dedup removed ", before - after, " row(s)")
    orig_samples[, dedup_key := NULL]
  }

  next_id <- res$next_id

  .p("[RUN] Finished ", batch_key, " | added=", ifelse(is.null(res$new_samples_tbl), 0, nrow(res$new_samples_tbl)))
  write_batch_outputs(batch_key, res$debug_tbl, res$map_tbl)
}

.p("[WRITE] Writing updated samples -> ", prev_samples_path)
fwrite(orig_samples, prev_samples_path)
.p("[WRITE] OK")

.p("============================================================")
.p("[DONE] Completed all batches. Updated samples at: ", prev_samples_path)
.p("============================================================")
