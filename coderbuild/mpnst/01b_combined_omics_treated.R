# 01b_combined_omics_treated.R
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(synapser)
  library(dplyr)
  library(tidyr)
})

# -----------------------
# logging helpers
# -----------------------
.ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
.p  <- function(...) cat(.ts(), paste0(...), "\n")

# -----------------------
# batch inputs (Synapse)
# -----------------------
B1_SYN <- "syn64608366"  # salmon.merged.gene_tpm_corrected_id.csv
B2_SYN <- "syn64608370"  # salmon.merged.gene_tpm.tsv_corrected_ID.csv
B3_SYN <- "syn65887795"  # salmon.merged.gene_tpm.tsv
B4_SYN <- "syn68898091"  # salmon.merged.gene_tpm.tsv

# -----------------------
# helpers
# -----------------------
study_label <- function(type) {
  dplyr::case_when(
    type == "patient derived xenograft"     ~ "MPNST PDX",
    type == "tumor"                          ~ "MPNST Tumor",
    type == "xenograft derived organoid"     ~ "MPNST PDX MT",
    TRUE                                     ~ "MPNST"
  )
}

is_treated_microtissue <- function(other_id, other_names) {
  oid <- as.character(other_id)
  onm <- as.character(other_names)
  grepl("treated_microtissue$", oid) && !is.na(onm) && nzchar(trimws(onm))
}

split_rep_ids <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return(character(0))
  v <- unlist(strsplit(x, ","))
  v <- trimws(v)
  v <- v[nzchar(v)]
  unique(v)
}

guess_batch_from_rep <- function(rep_id) {
  r <- as.character(rep_id)
  if (grepl("^b1_", r)) return("b1")
  if (grepl("^b2_", r)) return("b2")
  if (grepl("^B3-", r) || grepl("^B3_", r)) return("b3")
  if (grepl("^B4-", r) || grepl("^B4\\.", r)) return("b4")
  return(NA_character_)
}

rep_to_colname <- function(rep_id, batch) {
  r <- as.character(rep_id)
  if (batch %in% c("b1", "b2")) return(r)
  if (batch == "b3") return(gsub("-", "_", r, fixed = FALSE))        # B3-10 -> B3_10
  if (batch == "b4") return(gsub("B4-", "B4.", r, fixed = TRUE))     # B4-10 -> B4.10
  return(r)
}

read_tpm_matrix <- function(syn_id) {
  .p("[SYN] synGet(", syn_id, ") ...")
  pth <- synGet(syn_id)$path
  .p("[READ] ", syn_id, " -> ", pth)
  dt <- fread(pth)
  stopifnot(all(c("gene_id", "gene_name") %in% names(dt)))
  dt
}

# returns a data.table with columns:
#   gene_id, transcriptomics, improve_sample_id, source, study
extract_treated_from_batch <- function(mat, treated_tbl, batch_key) {
  if (nrow(treated_tbl) == 0) return(data.table())

  tmp <- copy(treated_tbl)
  tmp[, rep_ids := lapply(other_names, split_rep_ids)]
  tmp[, rep_batch := vapply(rep_ids, function(v) {
    if (length(v) == 0) return(NA_character_)
    guess_batch_from_rep(v[1])
  }, character(1))]

  tmp <- tmp[rep_batch == batch_key]
  if (nrow(tmp) == 0) return(data.table())

  tmp[, rep_cols := lapply(rep_ids, function(v) rep_to_colname(v, batch_key))]

  .p("[", batch_key, "] treated samples in this batch: ", nrow(tmp))

  out_list <- vector("list", nrow(tmp))

  for (i in seq_len(nrow(tmp))) {
    samp_id    <- tmp$improve_sample_id[i]
    samp_model <- tmp$model_type[i]
    cols       <- tmp$rep_cols[[i]]
    cols       <- cols[cols %in% names(mat)]

    if (length(cols) == 0) {
      .p("[", batch_key, "][SKIP] improve_sample_id=", samp_id,
         " -> no replicate columns found in matrix")
      next
    }

    expr <- rowMeans(as.matrix(mat[, ..cols]), na.rm = TRUE)

    out_list[[i]] <- data.table(
      gene_id           = mat$gene_id,
      transcriptomics   = expr,
      improve_sample_id = samp_id,
      source            = "NF Data Portal",
      study             = study_label(samp_model)
    )
  }

  rbindlist(out_list, use.names = TRUE, fill = TRUE)
}

# -----------------------
# MAIN
# -----------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 01b_combined_omics_treated.R <PAT> <samples.csv> <genes.csv>", call. = FALSE)
}
PAT     <- args[1]
samples <- args[2]
genes   <- args[3]

.p("============================================================")
.p("=== 01b_combined_omics_treated.R ===")
.p("Timestamp: ", .ts())
.p("Args: ", paste(args, collapse = " | "))
.p("============================================================")

.p("[SYN] synLogin() ...")
synLogin(authToken = PAT)
.p("[SYN] synLogin() OK")

.p("[READ] samples: ", samples)
samples_df <- fread(samples)

treated <- samples_df[
  model_type == "xenograft derived organoid" &
    vapply(seq_len(nrow(samples_df)), function(i) {
      is_treated_microtissue(samples_df$other_id[i], samples_df$other_names[i])
    }, logical(1)),
  .(other_names, improve_sample_id, model_type)
]

.p("[SAMPLES] treated_microtissue rows: ", nrow(treated))
if (nrow(treated) == 0) stop("No treated microtissue samples found in samples.csv", call. = FALSE)

.p("[READ] genes: ", genes)
genes_df <- fread(genes)
stopifnot(all(c("other_id", "entrez_id") %in% names(genes_df)))

# matrices
m1 <- read_tpm_matrix(B1_SYN)
m2 <- read_tpm_matrix(B2_SYN)
m3 <- read_tpm_matrix(B3_SYN)
m4 <- read_tpm_matrix(B4_SYN)

# extract
t1 <- extract_treated_from_batch(m1, treated, "b1")
t2 <- extract_treated_from_batch(m2, treated, "b2")
t3 <- extract_treated_from_batch(m3, treated, "b3")
t4 <- extract_treated_from_batch(m4, treated, "b4")

treated_long <- rbindlist(list(t1, t2, t3, t4), use.names = TRUE, fill = TRUE)
.p("[OUT] treated_long rows: ", nrow(treated_long))
if (nrow(treated_long) == 0) stop("No treated data extracted (check replicate IDs vs column names).", call. = FALSE)

# ENSEMBL gene_id has version -> strip to match genes_df$other_id
treated_long[, other_id := tstrsplit(gsub("\\s+", "", gene_id), "\\.", keep = 1)]

treated_final <- treated_long %>%
  left_join(genes_df, by = "other_id") %>%
  dplyr::select(entrez_id, transcriptomics, improve_sample_id, source, study) %>%
  filter(!is.na(entrez_id), transcriptomics != 0) %>%
  distinct()

# enforce strict schema & order
treated_final <- as.data.table(treated_final)[, .(entrez_id, transcriptomics, improve_sample_id, source, study)]

treated_out_path <- file.path("/tmp", "mpnst_transcriptomics_treated_microtissue.csv")
.p("[WRITE] treated-only transcriptomics -> ", treated_out_path)
fwrite(treated_final, treated_out_path)
.p("[WRITE] OK")

# optional append into main transcriptomics (same strict schema)
main_path <- file.path("/tmp", "mpnst_transcriptomics.csv")
if (file.exists(main_path)) {
  .p("[APPEND] Found existing main transcriptomics: ", main_path)
  main_df <- fread(main_path)
  main_df <- main_df[, .(entrez_id, transcriptomics, improve_sample_id, source, study)]  # enforce schema

  combined <- rbindlist(list(main_df, treated_final), use.names = TRUE, fill = TRUE)

  # de-dup conservatively on (entrez_id, improve_sample_id, source, study, transcriptomics)
  combined[, dedup_key := paste(entrez_id, improve_sample_id, source, study, transcriptomics, sep = "||")]
  combined <- combined[!duplicated(dedup_key)]
  combined[, dedup_key := NULL]

  .p("[WRITE] overwriting main transcriptomics with appended treated rows -> ", main_path)
  fwrite(combined, main_path)
  .p("[WRITE] OK")
} else {
  .p("[APPEND] No existing /tmp/mpnst_transcriptomics.csv found. Not appending.")
}

.p("============================================================")
.p("[DONE] Wrote treated transcriptomics with strict schema.")
.p("============================================================")
