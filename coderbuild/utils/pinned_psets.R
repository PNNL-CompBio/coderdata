## ---------------------------------------------------------------------------
## Pinned PharmacoGx PSet list for the broad_sanger build.
##
## WHY THIS EXISTS
## ---------------
## PharmacoGx::availablePSets() is a LIVE list served by ORCESTRA. New PSets
## appear without notice, and the build previously consumed whatever that call
## returned. That is unsafe for a published dataset in two distinct ways:
##
##   1. A new PSet with an INCOMPATIBLE schema breaks the build.
##      In Aug 2026 a second CTRPv2 PSet ("CTRPv2.1_2016") appeared. Its
##      sensitivityInfo() has neither `treatmentid` nor `drugid` -- it uses
##      `Compound.Name` / `Broad.Compound.ID` / `master_cpd_id` instead. The
##      drugs step found 0 chemicals for it and then crashed in
##      update_dataframe_and_write_tsv with
##      "TypeError: 'NoneType' object is not iterable".
##
##   2. A new PSet with a COMPATIBLE schema silently changes the published
##      data between releases, with no diff and no reviewer.
##
## Both the drugs step (03-createDrugFile.R) and the experiments step
## (04a-drugResponseData.R) resolve PSets through this file, so they can never
## disagree about which PSets are in a release. A mismatch there would produce
## experiments referencing drugs that do not exist in the drugs file, which
## fails schema validation late and confusingly.
##
## ADDING A NEW PSet
## -----------------
## Adding a PSet changes the published dataset, so it is a deliberate,
## reviewed change -- not something that should happen because upstream
## published something new. To add one:
##   1. Confirm its sensitivityInfo() exposes `treatmentid` or `drugid`
##      (or add an explicit column mapping in 03-createDrugFile.R).
##   2. Confirm the experiments step can extract dose-response from it.
##   3. Add the exact `PSet Name` below and note the release it entered in.
##
## CURRENT PINS
## ------------
## These are exactly the PSets used to build coderdata v2.3 (verified against
## the v2.3 build log), so v2.4 keeps the same PharmacoGx provenance.
##
## KNOWN AVAILABLE BUT DELIBERATELY EXCLUDED:
##   CTRPv2.1_2016 -- DO NOT ADD. It is not new data; it is the SAME CTRP v2
##                    data as CTRPv2_2015 with different curation. Verified
##                    2026-08-28 by loading both PSets:
##
##                      sensitivityInfo rows : 395263  vs 395263  (identical)
##                      unique samples       :    887  vs    887  (identical)
##                      unique experiment_id :    907  vs    907  (identical
##                                                       SETS -- identical==TRUE)
##                      unique treatments    :    544  vs    545
##
##                    The only real difference is naming convention.
##                    CTRPv2_2015 is PharmacoGx-standardised (`treatmentid`,
##                    `sampleid`, INN names like "adavosertib", hyphenated
##                    codes like "azd-1480"). CTRPv2.1_2016 keeps raw CTRP
##                    columns (`Compound.Name`, `master_ccl_id`) and Broad
##                    code names ("abt-199", "azd1480"). Case-insensitively
##                    481 of ~545 treatment names match outright; most of the
##                    remainder are the same compound spelled differently
##                    (azd-1480 / azd1480, azd-7545 / azd7545).
##
##                    Including it would therefore DUPLICATE all 395,263
##                    measurements and mint a second improve_drug_id for
##                    compounds we already have under a different spelling --
##                    a data-quality regression, not an enrichment.
## ---------------------------------------------------------------------------

PINNED_PSETS <- c(
  "CTRPv2_2015",
  "CCLE_2015",
  "CCLE_2019",
  "FIMM_2016",
  "gCSI_2019",
  "GDSC_2020(v1-8.2)",
  "GDSC_2020(v2-8.2)",
  "PRISM_2020"
)


## Return the PSet names for one dataset, restricted to the pinned list.
##
## Loudly reports any PSet that upstream now offers but which is not pinned,
## so a new upstream release is visible in the build log rather than silently
## included or silently dropped.
psets_for_dataset <- function(all.dsets, cel) {
  available <- all.dsets |>
    subset(`Dataset Name` == cel) |>
    (\(d) d[["PSet Name"]])() |>
    unlist() |>
    unique()

  keep <- available[available %in% PINNED_PSETS]
  skipped <- setdiff(available, PINNED_PSETS)

  if (length(skipped) > 0) {
    message(sprintf(
      paste0("[pinned_psets] NOTE: dataset '%s' has %d PSet(s) available ",
             "upstream that are NOT pinned and are being SKIPPED: %s\n",
             "  If this PSet should be part of the release, add it to ",
             "PINNED_PSETS in coderbuild/utils/pinned_psets.R after verifying ",
             "its schema in BOTH the drugs and experiments steps."),
      cel, length(skipped), paste(skipped, collapse = ", ")))
  }

  if (length(keep) == 0) {
    stop(sprintf(
      paste0("No pinned PSet found for dataset '%s'. Available upstream: %s. ",
             "Either upstream renamed/removed a PSet, or PINNED_PSETS in ",
             "coderbuild/utils/pinned_psets.R is out of date."),
      cel, if (length(available)) paste(available, collapse = ", ") else "(none)"))
  }

  keep
}
