## Build /tmp/genes.csv — entrez_id, gene_symbol, other_id, other_id_source.
##
## Uses org.Hs.eg.db for the local mappings (entrez ↔ symbol, ensembl gene,
## ensembl transcript) and biomaRt to filter to protein-coding genes.
##
## biomaRt's call to Ensembl frequently fails with mirror unavailability.
## To handle this, the Ensembl call is wrapped in a retry loop:
##   - cycles through the three valid mirror choices ('www', 'useast', 'asia')
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

# ---- Ensembl "new site" compatibility (REQUIRED) ----------------------------
#
# In Aug 2026 Ensembl migrated all sites to the new ensembl.org. The new site
# performs Accept-header content negotiation on the legacy BioMart path:
#
#   Accept: */*                                     -> 308 redirect to the
#   Accept: text/xml                                   BioMart archive host
#   Accept: application/json                        -> 404 {"details":
#   Accept: application/json, text/xml, ... (httr)     "No supported new
#                                                       Ensembl equivalent
#                                                       for this URL"}
#
# httr -- which biomaRt uses -- sends `application/json` FIRST by default, so
# every biomaRt call gets a hard 404 and biomaRt reports it as the misleading
# "Ensembl site unresponsive". This is deterministic: without this override the
# retry loop below can never succeed, no matter how many cycles it runs.
#
# Forcing Accept: text/xml makes Ensembl issue the 308 to the BioMart archive,
# which is what biomaRt actually wants. Verified 2026-08-28: with this line the
# protein-coding query succeeds; without it, it fails 100% of the time.
httr::set_config(httr::add_headers(Accept = "text/xml"))

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
#
# Retries still matter even with the Accept fix above: during the Aug 2026
# migration the BioMart backends are load-balanced and only some of them
# answer -- roughly 1 request in 5 returns real MartRegistry XML, the rest
# return a "Service unavailable" page (served, unhelpfully, with HTTP 200).
# Across MAX_CYCLES x MIRRORS attempts that is very likely to land.
#
# NOTE: 'uswest' is NOT a valid biomaRt mirror -- biomaRt rejects it with
# "Invalid mirror. Select a mirror from [www, useast, asia]" and silently falls
# back to www, wasting a quarter of every cycle. Only the three real mirrors
# are listed here.
MIRRORS     <- c("www", "useast", "asia")
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


# ---- Fallback 1: query BioMart directly over HTTP ---------------------------
#
# This is a FALLBACK ONLY. biomaRt above stays the primary path so that when
# Ensembl and biomaRt are both healthy we use the supported, maintained API.
# This runs only when every biomaRt mirror attempt has already failed.
#
# WHY IT IS NEEDED
# ----------------
# As of 2026-09-01 biomaRt (2.62.1) cannot reach Ensembl for two independent
# reasons, neither of which is fixable from our side:
#
#   1. Ensembl's martservice no longer accepts POST, which is what biomaRt
#      sends:
#        POST https://www.ensembl.org/biomart/martservice  -> 405 Method Not Allowed
#        GET  https://www.ensembl.org/biomart/martservice?query=... -> 200
#
#   2. biomaRt's own health check, .test_ensembl(), performs a request against
#      BOTH www.ensembl.org and useast.ensembl.org and raises if EITHER fails.
#      useast currently returns 502 and asia returns 403, so a perfectly
#      healthy www is still reported as "Ensembl site unresponsive".
#
# Both are being reported upstream (Huber-group-EMBL/biomaRt). When they are
# fixed, biomaRt will simply succeed above and this code will never execute --
# so it is safe to leave in place, and it should NOT be promoted to the
# primary path.
#
# The query below is the exact equivalent of the getBM() call above: the
# hsapiens_gene_ensembl dataset, filtered to biotype = protein_coding,
# returning ensembl_gene_id. Verified 2026-09-01 to return 23,272 genes,
# matching what biomaRt returned when it last worked (2026-08-29).
#
# GET is safe here despite BioMart queries normally using POST: this query
# encodes to ~500 characters, far below the ~8 KB practical URL limit. A query
# with 30 attributes still only reaches ~1,700.
fetch_protein_coding_direct <- function() {
  if (!requireNamespace("httr2", quietly = TRUE)) return(NULL)
  query_xml <- paste0(
    '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>',
    '<Query virtualSchemaName="default" formatter="TSV" header="0" count="" ',
    'datasetConfigVersion="0.6">',
    '<Dataset name="hsapiens_gene_ensembl" interface="default">',
    '<Filter name="biotype" value="protein_coding"/>',
    '<Attribute name="ensembl_gene_id"/>',
    '</Dataset></Query>')

  for (attempt in seq_len(3)) {
    message(sprintf("[%s] Direct BioMart query, attempt %d/3...", Sys.time(), attempt))
    res <- tryCatch({
      req <- httr2::req_timeout(
        httr2::req_url_query(
          httr2::request("https://www.ensembl.org/biomart/martservice"),
          query = query_xml),
        300)
      body <- httr2::resp_body_string(httr2::req_perform(req))

      # BioMart signals failure in-band with a 200, so check the payload
      # rather than trusting the status code.
      if (grepl("Query ERROR|status.ensembl.org|<html", body, ignore.case = TRUE)) {
        stop("BioMart returned an error page rather than results")
      }
      ids <- trimws(strsplit(body, "\n")[[1]])
      ids <- ids[nzchar(ids)]
      if (length(ids) < 1000) {
        stop(sprintf("implausibly few genes returned (%d); treating as failure",
                     length(ids)))
      }
      data.frame(ensembl_gene_id = ids, stringsAsFactors = FALSE)
    }, error = function(e) {
      message(sprintf("[%s] Direct BioMart query failed: %s",
                      Sys.time(), conditionMessage(e)))
      NULL
    })
    if (!is.null(res)) {
      message(sprintf("[%s] Direct BioMart query succeeded (%d protein-coding genes)",
                      Sys.time(), nrow(res)))
      return(res)
    }
    if (attempt < 3) Sys.sleep(30)
  }
  NULL
}

if (is.null(tab)) {
  message("\n[fallback] biomaRt could not reach Ensembl; trying a direct ",
          "BioMart HTTP query before giving up.\n")
  tab <- fetch_protein_coding_direct()
}


# ---- Fallback 2: protein-coding inference from org.Hs.eg.db -----------------
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

# The fallback is OPT-IN. It used to run automatically, which meant a transient
# Ensembl problem would silently produce a genes.csv with a broader,
# non-canonical gene list -- and because the build then SUCCEEDED, that
# degraded file could flow all the way through validation to a published
# Figshare release with nobody noticing. Every downstream omics file joins
# against genes.csv, so this is not a local defect: it changes the whole
# release.
#
# Default is now to fail loudly. Set ALLOW_GENE_FALLBACK=1 to deliberately
# accept the degraded list (e.g. for a local debug build that is not going to
# be published).
if (is.null(tab)) {
  allow_fallback <- nzchar(Sys.getenv("ALLOW_GENE_FALLBACK"))
  if (!allow_fallback) {
    stop(paste0(
      "\n",
      "=====================================================================\n",
      " Ensembl unreachable -- REFUSING to build a degraded genes.csv\n",
      "=====================================================================\n",
      " All Ensembl mirrors failed after ", MAX_CYCLES, " cycles, AND the\n",
      " direct BioMart HTTP query fallback also failed.\n",
      "\n",
      " A local org.Hs.eg.db fallback exists, but it produces a BROADER,\n",
      " non-canonical gene list (it includes pseudogenes and other biotypes\n",
      " that the Ensembl protein-coding filter excludes). Every omics file\n",
      " joins against genes.csv, so publishing a release built on the\n",
      " fallback would silently change every dataset.\n",
      "\n",
      " This build has therefore been stopped rather than continuing quietly.\n",
      "\n",
      " Options:\n",
      "   - Wait for Ensembl and re-run (check https://status.ensembl.org)\n",
      "   - For a local, NON-published debug build only, re-run with:\n",
      "         ALLOW_GENE_FALLBACK=1\n",
      "\n",
      " NOTE: if every attempt failed with a 404 mentioning \"No supported new\n",
      " Ensembl equivalent\", the Accept-header override at the top of this\n",
      " script is missing or was undone -- see the comment there.\n",
      "=====================================================================\n"
    ))
  }
  message("\n=== WARNING ===")
  message("All Ensembl mirrors unreachable after ", MAX_CYCLES,
          " cycles. ALLOW_GENE_FALLBACK is set, so using the local ",
          "org.Hs.eg.db fallback.")
  message("This produces a BROADER gene list than the canonical ",
          "Ensembl-filtered version. DO NOT PUBLISH a release built this way; ",
          "re-run when Ensembl is reachable.")
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