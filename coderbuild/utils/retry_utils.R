## ---------------------------------------------------------------------------
## Retry helpers for flaky network fetches.
##
## WHY THIS EXISTS
## ---------------
## A full build makes hundreds of requests to a dozen third-party hosts over
## roughly a day. Any one of them being briefly unavailable used to fail the
## whole step, and because build_all.py restarts a failed step from the top,
## a single transient error could throw away many hours of work.
##
## Observed in real runs:
##   2026-09-02  Zenodo returned 504 Gateway Timeout on PSet_gCSI2019.rds,
##               after CTRPv2 and GDSC had already downloaded (one PSet is
##               1.5 GB). The drugs step had run for ~17 hours across three
##               attempts before giving up. The URL was fine minutes later.
##   2026-08-30  cog.sanger.ac.uk failed twice with an OpenSSL "unexpected
##               eof" before succeeding on the third attempt.
##
## These are transient upstream problems, not build defects, so the build
## should ride them out rather than fail.
##
## NOTE: this only covers fetches that were previously unprotected.
## 02-broadSangerOmics.R already has robust_download_httr2() for its own
## downloads, and the Python builders have their own retry logic.
## ---------------------------------------------------------------------------

## Waits between retries, in seconds: 1, 3, 10, then 15 minutes.
## Five attempts total (the initial one plus four retries), ~29 minutes of
## waiting if every attempt fails. Long enough to ride out a real upstream
## incident rather than only a momentary blip.
RETRY_SLEEPS <- c(60, 180, 600, 900)

## Run `expr` (a zero-argument function), retrying on error.
##
## IMPORTANT: if every attempt fails this STOPS. It never returns a partial or
## empty result, and callers must not catch the error to "carry on without
## this dataset" -- a dataset must never be silently dropped from a release.
## Exhausting the retries is meant to fail the step, and therefore the build,
## so the problem is visible rather than shipped.
##
##   what   - human-readable label used in log lines
##   sleeps - vector of waits between attempts; length(sleeps)+1 attempts
with_retries <- function(expr, what = "operation", sleeps = RETRY_SLEEPS) {
  tries <- length(sleeps) + 1L
  last_err <- NULL
  for (attempt in seq_len(tries)) {
    result <- tryCatch(
      list(ok = TRUE, value = expr()),
      error = function(e) list(ok = FALSE, msg = conditionMessage(e)))

    if (isTRUE(result$ok)) {
      if (attempt > 1) {
        message(sprintf("[%s] %s: succeeded on attempt %d/%d",
                        Sys.time(), what, attempt, tries))
      }
      return(result$value)
    }

    last_err <- result$msg
    if (attempt < tries) {
      sleep_s <- sleeps[attempt]
      message(sprintf("[%s] %s: attempt %d/%d failed (%s); retrying in %d min",
                      Sys.time(), what, attempt, tries,
                      substr(gsub("[\r\n]+", " ", last_err), 1, 120),
                      round(sleep_s / 60)))
      Sys.sleep(sleep_s)
    }
  }
  stop(sprintf(paste0(
    "%s failed after %d attempts over ~%d minutes. Last error: %s\n",
    "Failing the build rather than continuing without this data."),
    what, tries, round(sum(sleeps) / 60), last_err))
}


## Download a PharmacoGx PSet with retries, caching it so that a later attempt
## -- including a completely fresh container started by build_all.py's own
## retry loop -- does not re-download it.
##
## PSETS are large (CTRPv2 is 228 MB, one of them is 1.5 GB). downloadPSet()
## skips a file that is already present in saveDir, so pointing saveDir at a
## directory under /tmp (which build_all.py bind-mounts to the host's local/)
## makes the download survive across container runs. A subdirectory is used
## deliberately: build_all.py's intermediate cleanup globs `local/*.*` and
## skips directories, so the cache is not swept away mid-build.
PSET_CACHE_DIR <- Sys.getenv("PSET_CACHE_DIR", unset = "/tmp/pset_cache")

## Download a PSet, retrying, and never let a partial file poison the cache.
##
## PharmacoGx::downloadPSet only downloads when the file is ABSENT:
##
##   if (!file.exists(file.path(saveDir, pSetFileName))) { download.file(...) }
##   pSet <- readRDS(file.path(saveDir, pSetFileName))
##
## So an interrupted download leaves a truncated .rds that downloadPSet then
## reuses forever: every later attempt skips the download and fails instantly
## inside readRDS() with "error reading from connection". That is not a network
## error, and no amount of retrying -- or waiting for the remote host to come
## back -- can clear it. The cache poisons itself permanently.
##
## This happened on 2026-09-09. A Zenodo outage truncated CCLE_2015.rds at
## 2.4GB, and the drugs step then failed 15 times (5 retries x 3 step attempts)
## with an identical error, while the other five PSets downloaded fine.
##
## The fix is to delete the cached file whenever an attempt fails, so the next
## attempt performs a real download. We deliberately do NOT try to validate the
## file instead: the truncated CCLE file is still a readable gzip stream that
## R's gzfile() walks to the end without complaint, so a cheap integrity check
## reports it healthy. Only readRDS() detects it, and that would mean
## deserialising several GB just to decide whether to keep a file we are about
## to re-download anyway.
## Download any URL to `dest`, resuming an interrupted transfer.
##
## Generic counterpart to the PSet fetcher below, for the plain file downloads
## in the omics steps. Those are the largest transfers in the whole pipeline --
## Sanger ships rnaseq_all at ~897MB and WES_pureCN_CNV_genes at ~935MB -- and
## httr2's req_retry() restarts the whole request rather than resuming it. On a
## link that drops part-way, restarting a 935MB file simply fails again at a
## similar point: build v24 died exactly that way on a 112MB file, breaking
## 12MB in, three times.
##
## wget -c issues a Range request and continues from what is already on disk.
## Bytes accumulate in "<dest>.part", which is renamed onto `dest` only once the
## size matches the server's, so a partial file never appears at the
## destination. The ".part" is deliberately kept between attempts -- deleting it
## would mean a flaky link could never accumulate a large file.
download_file_resumable <- function(url, dest, sleeps = RETRY_SLEEPS,
                                    what = NULL) {
  what <- if (is.null(what)) basename(dest) else what
  part <- paste0(dest, ".part")

  expected <- suppressWarnings(tryCatch({
    out <- system2("wget", c("--spider", "-S", "--timeout=30", "--tries=2",
                             shQuote(url)), stdout = TRUE, stderr = TRUE)
    hits <- grep("Content-Length:", out, value = TRUE, ignore.case = TRUE)
    if (length(hits)) {
      as.numeric(sub(".*[Cc]ontent-[Ll]ength:[[:space:]]*([0-9]+).*", "\\1",
                     hits[length(hits)]))
    } else NA_real_
  }, error = function(e) NA_real_))

  with_retries(
    function() {
      status <- suppressWarnings(tryCatch(
        system2("wget", c("-c", "--tries=10", "--waitretry=15",
                          "--read-timeout=60", "--timeout=30",
                          "--progress=dot:giga",
                          "-O", shQuote(part), shQuote(url))),
        error = function(e) 1L))
      ok <- identical(as.integer(status), 0L)

      if (!file.exists(part)) {
        stop(sprintf("download of %s produced no file", what))
      }
      got <- file.info(part)$size

      if (!is.na(expected)) {
        if (got < expected) {
          stop(sprintf("%s incomplete (%s of %s bytes); next attempt resumes",
                       what, format(got, big.mark = ","),
                       format(expected, big.mark = ",")))
        }
        if (got > expected) {
          unlink(part)
          stop(sprintf("%s is oversized (%s vs %s bytes); discarded", what,
                       format(got, big.mark = ","),
                       format(expected, big.mark = ",")))
        }
      } else if (!ok) {
        stop(sprintf("%s transfer did not complete (%s bytes so far, remote size unknown); next attempt resumes",
                     what, format(got, big.mark = ",")))
      }

      if (!file.rename(part, dest)) {
        stop(sprintf("could not move %s into place", what))
      }
      message(sprintf("[%s] %s: downloaded %s bytes", Sys.time(), what,
                      format(got, big.mark = ",")))
      invisible(dest)
    },
    what   = paste0("download('", what, "')"),
    sleeps = sleeps)
}


## Look up a PSet's download URL from the PharmacoGx catalogue.
pset_download_url <- function(pset_name) {
  tbl <- tryCatch(PharmacoGx::availablePSets(), error = function(e) NULL)
  if (is.null(tbl)) return(NA_character_)
  i <- match(pset_name, as.character(tbl[, "PSet Name"]))
  if (is.na(i)) return(NA_character_)
  as.character(tbl[i, "Download"])
}

## Fetch the publisher's size and MD5 for a PSet file.
##
## PSets are hosted on Zenodo, whose REST API publishes both. Returns a list
## with NA entries when the URL is not a Zenodo record or the API is
## unreachable, so callers degrade to size-only checking rather than failing.
pset_remote_meta <- function(url) {
  none <- list(size = NA_real_, md5 = NA_character_)
  if (is.na(url) || !nzchar(url)) return(none)
  rec <- sub(".*/records?/([0-9]+)/files/.*", "\\1", url)
  if (identical(rec, url) || !grepl("^[0-9]+$", rec)) return(none)
  fname <- basename(sub("\\?.*$", "", url))
  api <- paste0("https://zenodo.org/api/records/", rec)
  txt <- suppressWarnings(tryCatch(
    system2("wget", c("-qO-", "--timeout=30", "--tries=3", shQuote(api)), stdout = TRUE),
    error = function(e) character(0)))
  if (!length(txt)) return(none)
  j <- tryCatch(jsonlite::fromJSON(paste(txt, collapse = "")), error = function(e) NULL)
  if (is.null(j) || is.null(j$files) || !nrow(j$files)) return(none)
  i <- match(fname, j$files$key)
  if (is.na(i)) return(none)
  list(size = as.numeric(j$files$size[i]),
       md5  = sub("^md5:", "", as.character(j$files$checksum[i])))
}

## Fetch with wget -c, which RESUMES an interrupted transfer via an HTTP range
## request instead of starting over.
##
## The CCLE PSet is 2.8GB and the connection can drop between 400MB and 1.1GB,
## so a non-resuming download never reaches the end: each retry throws away
## everything it had and dies at roughly the same place. Verified 2026-09-09
## that Zenodo answers range requests with 206 PARTIAL_CONTENT.
fetch_pset_resumable <- function(url, dest) {
  if (is.na(url) || !nzchar(url)) return(FALSE)
  status <- suppressWarnings(tryCatch(
    system2("wget", c("-c", "--tries=10", "--waitretry=15", "--read-timeout=60",
                      "--timeout=30", "--progress=dot:giga",
                      "-O", shQuote(dest), shQuote(url))),
    error = function(e) 1L))
  identical(as.integer(status), 0L)
}

## Load a cached PSet the way PharmacoGx::downloadPSet would, minus its final
## saveRDS() -- that re-save is only a caching optimisation, and it is what
## makes a verified file stop matching its published checksum.
##
## An unreadable file is discarded so the caller can fetch a clean copy. This
## is the last line of defence for files cached before verification existed.
pset_load_cached <- function(cached, pset_name) {
  pSet <- tryCatch(readRDS(cached), error = function(e) {
    message(sprintf("[%s] cached PSet '%s' is unreadable (%s); discarding it.",
                    Sys.time(), pset_name, conditionMessage(e)))
    unlink(cached)
    stop(e)
  })
  updateObject(pSet)
}

## Attempts with PharmacoGx's own downloader before the fallback kicks in.
## Short on purpose: long enough to ride out a blip, short enough that a real
## outage reaches the resumable path in minutes rather than half an hour.
PRIMARY_SLEEPS <- utils::head(RETRY_SLEEPS, 2)

## Download a PSet.
##
## Primary path: PharmacoGx::downloadPSet, so normal builds behave exactly as
## upstream intends and we do not diverge from the library without cause.
##
## Fallback path (only after the primary has exhausted PRIMARY_SLEEPS): fetch
## with wget -c into "<name>.rds.part" and rename to "<name>.rds" only once
## size and MD5 match what Zenodo publishes. That separation is what makes
## resuming safe -- a ".part" file is by construction an interrupted download,
## so resuming it is always correct, while "<name>.rds" is by construction
## complete, so it is never resumed onto.
##
## The fallback exists because downloadPSet cannot recover from two situations
## seen in real builds:
##
##   - It only downloads when the file is ABSENT, so an interrupted transfer
##     leaves a truncated .rds that it then reuses forever; every later attempt
##     fails instantly in readRDS() with "error reading from connection". A
##     Zenodo outage on 2026-09-09 truncated CCLE_2015.rds at 2.4GB of 2.8GB
##     and the drugs step then failed 15 times identically.
##   - It does not resume. CCLE is 2.8GB and the connection can drop between
##     400MB and 1.1GB, so every attempt discarded its progress and died in
##     roughly the same place.
download_pset_retry <- function(pset_name, sleeps = RETRY_SLEEPS) {
  if (!dir.exists(PSET_CACHE_DIR)) {
    dir.create(PSET_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
  }
  cached <- file.path(PSET_CACHE_DIR, paste0(pset_name, ".rds"))
  part   <- paste0(cached, ".part")

  ## ---- Primary: PharmacoGx's own downloader -------------------------------
  primary <- tryCatch(
    with_retries(
      function() {
        ## Deliberately does NOT delete the cached file on failure. If the
        ## file is good and the error was transient, deleting it would throw
        ## away a multi-GB download for nothing. If the file IS truncated,
        ## downloadPSet reuses it blindly and the remaining primary attempts
        ## fail instantly on the same bytes -- which is harmless and reaches
        ## the fallback sooner. Either way the fallback sorts it out.
        PharmacoGx::downloadPSet(pset_name,
                                 saveDir = PSET_CACHE_DIR,
                                 timeout = 10000)
      },
      what   = paste0("downloadPSet('", pset_name, "')"),
      sleeps = PRIMARY_SLEEPS),
    error = function(e) NULL)

  if (!is.null(primary)) return(primary)

  ## ---- Fallback: resumable, checksum-verified fetch ------------------------
  message(sprintf("[%s] PSet '%s': falling back to a resumable, checksum-verified download.",
                  Sys.time(), pset_name))
  ## Note we do not delete `cached` here. A file left by the primary may be
  ## perfectly good -- the primary can fail for reasons unrelated to it -- and
  ## discarding a multi-GB download on suspicion is expensive. pset_load_cached()
  ## below is the gate: if the file is unreadable it is removed there, and the
  ## next attempt downloads a clean, checksum-verified copy. Resuming onto it is
  ## impossible either way, because the fallback only ever resumes ".part".

  with_retries(
    function() {
      if (!file.exists(cached)) {
        url <- pset_download_url(pset_name)
        if (is.na(url) || !nzchar(url)) {
          stop(sprintf("no download URL for PSet '%s'", pset_name))
        }
        meta <- pset_remote_meta(url)
        complete <- fetch_pset_resumable(url, part)

        if (!file.exists(part)) {
          stop(sprintf("download of PSet '%s' produced no file", pset_name))
        }
        got <- file.info(part)$size

        ## A ".part" file is promoted to the final path ONLY on positive
        ## evidence that it is complete. Without that evidence it stays a
        ## ".part", so the next attempt RESUMES it.
        ##
        ## This is what failed on 2026-09-10 for CCLE_2019. Zenodo was timing
        ## out, so the metadata lookup timed out too and both checks below were
        ## skipped as NA. The partial file was then renamed into place,
        ## readRDS rejected it, and the file was deleted -- destroying the
        ## downloaded bytes on every cycle. Five attempts over 29 minutes each
        ## restarted from zero and the step could never finish.
        ##
        ## wget's exit status is the evidence when the API is unreachable: it
        ## is 0 only when the whole file was retrieved.
        if (!is.na(meta$size)) {
          if (got < meta$size) {
            stop(sprintf("PSet '%s' incomplete (%s of %s bytes); next attempt resumes",
                         pset_name, format(got, big.mark = ","),
                         format(meta$size, big.mark = ",")))
          }
          if (got > meta$size) {
            unlink(part)
            stop(sprintf("PSet '%s' is oversized (%s vs %s bytes); discarded",
                         pset_name, format(got, big.mark = ","),
                         format(meta$size, big.mark = ",")))
          }
        } else if (!complete) {
          ## No published size AND wget did not finish: assume partial. Keeping
          ## it is safe -- the worst case is one wasted resume attempt, whereas
          ## discarding it guarantees the download can never accumulate.
          stop(sprintf("PSet '%s' transfer did not complete (%s bytes so far, remote size unknown); next attempt resumes",
                       pset_name, format(got, big.mark = ",")))
        }

        if (!is.na(meta$md5)) {
          local_md5 <- unname(tools::md5sum(part))
          if (!identical(local_md5, meta$md5)) {
            unlink(part)
            stop(sprintf("PSet '%s' failed checksum (got %s, published %s); discarded",
                         pset_name, local_md5, meta$md5))
          }
          message(sprintf("[%s] PSet '%s' verified against published MD5 %s",
                          Sys.time(), pset_name, meta$md5))
        } else if (!is.na(meta$size)) {
          message(sprintf("[%s] PSet '%s': no published checksum; verified against published size %s bytes.",
                          Sys.time(), pset_name, format(meta$size, big.mark = ",")))
        } else {
          message(sprintf("[%s] PSet '%s': no published size or checksum available; accepted on a complete transfer, readRDS is the remaining check.",
                          Sys.time(), pset_name))
        }

        if (!file.rename(part, cached)) {
          stop(sprintf("could not move verified PSet '%s' into place", pset_name))
        }
      }
      pset_load_cached(cached, pset_name)
    },
    what   = paste0("downloadPSet('", pset_name, "') [resumable fallback]"),
    sleeps = sleeps)
}


## ---------------------------------------------------------------------------
## Per-PSet completion markers.
##
## Retrying the download is not enough on its own. build_all.py restarts a
## failed step from the very beginning, so before these markers existed, a
## failure on the 5th of 6 PSets meant re-downloading and re-processing the
## first four -- which is how one Zenodo 504 turned into a 17-hour failure
## across three attempts.
##
## A marker is written only after a PSet's results have been merged into the
## cumulative output, so a half-finished PSet is never marked done. Markers
## live under /tmp (build_all.py's bind mount to the host's local/) alongside
## the cumulative output they describe, so the two are cleared together and
## cannot disagree: wiping local/ resets both, and keeping local/ keeps both.
##
## Net effect: a rerun re-does only the PSet that actually failed.
## ---------------------------------------------------------------------------
PSET_STATE_DIR <- Sys.getenv("PSET_STATE_DIR",
                             unset = file.path(PSET_CACHE_DIR, "completed"))

pset_is_done <- function(tag) {
  file.exists(file.path(PSET_STATE_DIR, paste0(make.names(tag), ".done")))
}

mark_pset_done <- function(tag) {
  if (!dir.exists(PSET_STATE_DIR)) {
    dir.create(PSET_STATE_DIR, recursive = TRUE, showWarnings = FALSE)
  }
  file.create(file.path(PSET_STATE_DIR, paste0(make.names(tag), ".done")))
  invisible(TRUE)
}


## Canonical path for a cached per-study dose-response table (04a).
##
## Both the writer and the reader must agree on this name. They previously did
## not -- the writer used paste0(tolower(studyName),'DoseResponse') while the
## reader looked for paste0(cel,'doseResponse'), so the resume check never
## matched and every attempt recomputed everything. Routing both through one
## function makes that class of mismatch impossible.
##
## Kept under /tmp so it survives between container attempts, same rationale
## as the PSet cache above.
DOSE_RESPONSE_CACHE_DIR <- Sys.getenv("DOSE_RESPONSE_CACHE_DIR",
                                      unset = "/tmp/dose_response_cache")

dose_response_cache_path <- function(study_name) {
  if (!dir.exists(DOSE_RESPONSE_CACHE_DIR)) {
    dir.create(DOSE_RESPONSE_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
  }
  file.path(DOSE_RESPONSE_CACHE_DIR,
            paste0(tolower(study_name), "_doseResponse.tsv"))
}
