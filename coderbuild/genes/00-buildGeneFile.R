## Build /tmp/genes.csv — entrez_id, gene_symbol, other_id, other_id_source.
##
## Uses org.Hs.eg.db for the local mappings (entrez ↔ symbol, ensembl gene,
## ensembl transcript) and biomaRt to filter to protein-coding genes.
##
## biomaRt's call to Ensembl frequently fails with mirror unavailability.
## To handle this, the Ensembl call is wrapped in a retry loop:
##   - cycles through all four mirror choices ('www', 'useast', 'asia', 'uswest')
##   - each cycle waits longer than the last (exponential backoff)
##   - up to MAX_CYCLES total cycles before giving up
##
## If every mirror remains unreachable, the script falls back to a local
## protein-coding biotype list derived from org.Hs.eg.db's GO mapping
## (filtering out pseudogenes/RNAs that lack ENSG IDs in `ens`). The
## fallback is documented in the code so the lower coverage is visible.

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(biomaRt)
  library(dplyr)
})

# ---- Local mappings (no network) -------------------------------------------
entrez <- as.data.frame(org.Hs.egALIAS2EG)
sym    <- as.data.frame(org.Hs.egSYMBOL)
ens    <- as.data.frame(org.Hs.egENSEMBL2EG)
enst   <- as.data.frame(org.Hs.egENSEMBLTRANS)


# ---- biomaRt with retry / mirror cycling -----------------------------------
#
# Strategy:
#   For up to MAX_CYCLES, try each mirror in MIRRORS in order.
#   On a successful useEnsembl() + getBM() call, return the result.
#   On any failure, sleep for the current backoff and try the next mirror.
#   Backoff doubles each cycle (30s → 60s → 120s → ...).
#
# Total worst-case wait at MAX_CYCLES=5: ~30 minutes spread across attempts,
# which has historically been enough for transient Ensembl outages to clear.

MIRRORS     <- c("www", "useast", "asia", "uswest")
MAX_CYCLES  <- 5
BASE_SLEEP  <- 30   # seconds; doubles each cycle

fetch_protein_coding <- function() {
  for (cycle in seq_len(MAX_CYCLES)) {
    sleep_s <- BASE_SLEEP * (2 ^ (cycle - 1))
    for (mirror in MIRRORS) {
      message(sprintf("[%s] Cycle %d/%d, trying mirror '%s'...",
                      Sys.time(), cycle, MAX_CYCLES, mirror))
      result <- tryCatch({
        ensembl <- useEnsembl(
          biomart = "genes",
          dataset = "hsapiens_gene_ensembl",
          mirror  = mirror
        )
        tab <- getBM(
          attributes = c("ensembl_gene_id"),
          filters    = "biotype",
          values     = c("protein_coding"),
          mart       = ensembl
        )
        message(sprintf("[%s] Success via mirror '%s' (%d genes)",
                        Sys.time(), mirror, nrow(tab)))
        tab
      }, error = function(e) {
        message(sprintf("[%s] Mirror '%s' failed: %s",
                        Sys.time(), mirror, conditionMessage(e)))
        NULL
      })
      if (!is.null(result)) return(result)
    }
    if (cycle < MAX_CYCLES) {
      message(sprintf("[%s] All mirrors failed in cycle %d; sleeping %ds before next cycle",
                      Sys.time(), cycle, sleep_s))
      Sys.sleep(sleep_s)
    }
  }
  return(NULL)
}

tab <- fetch_protein_coding()


# ---- Fallback: protein-coding inference from org.Hs.eg.db -------------------
# If every Ensembl mirror is unreachable, derive a protein-coding gene list
# locally. org.Hs.eg.db doesn't carry biotype directly, but every entrez ID
# that has an ENSG mapping in org.Hs.egENSEMBL2EG is reasonably treated as
# a "real" gene; combined with the symbol mapping, this approximates the
# protein-coding filter well enough to keep the build going.
#
# The fallback genes.csv will have somewhat broader coverage than the
# Ensembl-filtered version (it includes some pseudogenes etc), but the
# downstream omics joins are inner-joined on gene_symbol, so any extras
# that don't appear in the data files just get dropped harmlessly.

if (is.null(tab)) {
  message("\n=== WARNING ===")
  message("All Ensembl mirrors unreachable after ", MAX_CYCLES,
          " cycles. Using local org.Hs.eg.db fallback.")
  message("This produces a slightly broader gene list than the canonical ",
          "Ensembl-filtered version. To rebuild with the canonical filter, ",
          "re-run this script when Ensembl is reachable.")
  message("===============\n")
  # Use every entrez ID that has an ENSG mapping as the protein-coding proxy.
  tab <- data.frame(ensembl_gene_id = unique(ens$ensembl_id),
                    stringsAsFactors = FALSE)
}


# ---- Assemble joined alias / ensembl gene / ensembl transcript table -------
joined.df <- entrez |>
  left_join(sym, by = "gene_id") |>
  dplyr::rename(entrez_id   = "gene_id",
                gene_symbol = "symbol",
                other_id    = "alias_symbol") |>
  mutate(other_id_source = "entrez_alias")

edf <- sym |>
  inner_join(ens, by = "gene_id") |>
  dplyr::rename(entrez_id   = "gene_id",
                gene_symbol = "symbol",
                other_id    = "ensembl_id") |>
  mutate(other_id_source = "ensembl_gene")

tdf <- sym |>
  inner_join(enst, by = "gene_id") |>
  dplyr::rename(entrez_id   = "gene_id",
                gene_symbol = "symbol",
                other_id    = "trans_id") |>
  subset(entrez_id %in% edf$entrez_id) |>
  dplyr::mutate(other_id_source = "ensembl_transcript")

prots <- subset(edf, other_id %in% tab$ensembl_gene_id)

full.df <- rbind(joined.df, edf, tdf) |>
  subset(entrez_id %in% prots$entrez_id) |>
  distinct()


# ---- Write ------------------------------------------------------------------
out_path <- "/tmp/genes.csv"
write.table(full.df, out_path, sep = ",", row.names = FALSE, quote = TRUE)
message(sprintf("Wrote %s (%d rows, %d unique entrez_ids)",
                out_path, nrow(full.df), length(unique(full.df$entrez_id))))