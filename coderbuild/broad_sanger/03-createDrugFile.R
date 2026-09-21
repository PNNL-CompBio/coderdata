###Here is a script that runs through all the data files one by one.
library('reticulate')
use_python("/opt/venv/bin/python3", required = TRUE)
library('tidyr')
library(dplyr)
library(data.table)
#this is a helper file that loads the data
source_python("pubchem_retrieval.py")

library('PharmacoGx')

## Pinned PSet list shared with 04a-drugResponseData.R
source("pinned_psets.R")
## Retry/caching helpers for flaky third-party downloads
source("retry_utils.R")

all.dsets<-PharmacoGx::availablePSets()



#' getCellLineData - gets cell line dose response data
getDepMapDrugData<-function(cell.lines=c('CTRPv2','FIMM','gCSI','PRISM','GDSC','CCLE'),efile=''){

    if(efile!=''){
        existing_ids=readr::read_tsv(efile)
    }else{
        existing_ids=NULL
        }
    output_file_path <- '/tmp/depmap_sanger_drugs.tsv'
    ignore_file_path <- '/tmp/ignore_chems.txt'
    for(cel in cell.lines){

        ## PSets are pinned (see coderbuild/utils/pinned_psets.R) rather than
        ## taken straight from availablePSets(), so a new upstream PSet cannot
        ## silently change or break the release.
        files<-psets_for_dataset(all.dsets, cel)

        for(f in files){
            print(f)
            if(f=='GDSC_2020(v2-8.2)')
                cel='GDSCv2'
            if(f=='GDSC_2020(v1-8.2)')
                cel='GDSCv1'

            ## Skip PSets already merged into the cumulative output on an
            ## earlier attempt, so a failure late in the list does not force a
            ## re-run of everything before it. The marker is written only after
            ## the merge below succeeds.
            done_tag <- paste0("drugs-", f)
            if (pset_is_done(done_tag)) {
                message(sprintf(
                  "[03-createDrugFile] PSet '%s' already processed in a previous attempt; skipping.", f))
                next
            }

            dset<<-download_pset_retry(f)


            mapping <- sensitivityInfo(dset)##get the dataset dose response data

            ##map the rownames
            if(!'exp_id'%in%names(mapping))
                mapping<-mapping%>%
                    tibble::rownames_to_column('exp_id')

            ##fix up the NSC ids
            if('NSC'%in%names(mapping))
                mapping<-mapping|>
                    dplyr::select(-treatmentid)|>
                    dplyr::mutate(treatmentid=paste0('NSC-',NSC))

            ##move drug to treatment id
            if("drugid"%in%names(mapping))
                mapping<-dplyr::rename(mapping,treatmentid='drugid')

            ##query to build the drug ids
            chem_list <- unique(mapping$treatmentid)
            if (!is.null(existing_ids)) {
                chem_list <- setdiff(chem_list, existing_ids$chem_name)
            }

            # testing: only keep first 10 per PSet
            # -------
            # chem_list <- sample(chem_list, min(10, length(chem_list)))
            # print(paste('Testing mode: using', length(chem_list), 'chemicals for dataset', cel))
            # -------

            print(paste('Found',length(chem_list),'chemicals for dataset',cel))

            ## Guard: never hand an empty/NULL chem_list to the Python writer.
            ## update_dataframe_and_write_tsv() iterates it and fails with the
            ## opaque "TypeError: 'NoneType' object is not iterable", which is
            ## how the Aug 2026 CTRPv2.1_2016 schema change surfaced.
            ##
            ## Zero is also a legitimate outcome -- e.g. every compound in this
            ## PSet was already covered by a previous PSet or by efile -- so
            ## this skips rather than errors. A PSet whose treatment column is
            ## missing entirely is reported loudly, since that means an
            ## unrecognised schema rather than "nothing new to do".
            ## An unrecognised schema is a HARD FAILURE, not a skip. PSets are
            ## pinned (coderbuild/utils/pinned_psets.R), so if a pinned PSet no
            ## longer exposes a treatment column it means upstream changed its
            ## shape. Continuing would silently omit that PSet's drugs from the
            ## release, which is exactly the kind of quiet data loss the pinning
            ## exists to prevent -- so stop and make it visible.
            if (is.null(mapping$treatmentid)) {
                stop(sprintf(
                  paste0("PSet '%s': sensitivityInfo() has no 'treatmentid' or ",
                         "'drugid' column (columns: %s).\n",
                         "A pinned PSet has changed schema. Refusing to continue, ",
                         "because skipping it would drop its drugs from the release.\n",
                         "Either add an explicit column mapping for this PSet or ",
                         "remove it from PINNED_PSETS after confirming that is correct."),
                  f, paste(colnames(mapping), collapse = ", ")))
            }

            ## Zero new chemicals is a legitimate outcome and NOT data loss:
            ## every compound in this PSet was already fetched via an earlier
            ## PSet in this run or via efile, so its drugs are already in the
            ## cumulative output. Nothing is dropped by moving on.
            if (length(chem_list) == 0) {
                message(sprintf(
                  "[03-createDrugFile] PSet '%s': all chemicals already present; nothing new to fetch.", f))
                mark_pset_done(done_tag)
                next
            }


            #build prev_drug_filepaths for this iteration (include original and accumulated) ---
            prev_list <- c()
            if (!is.na(efile) && nzchar(efile)) prev_list <- c(prev_list, efile)
            if (file.exists(output_file_path)) prev_list <- c(prev_list, output_file_path)
            prev_drug_filepaths <- if (length(prev_list) > 0) paste(prev_list, collapse = ",") else NULL

            # write this PSet's results to a temp file
            temp_out <- paste0(output_file_path, ".part.tsv")
            update_dataframe_and_write_tsv(
                unique_names           = chem_list,
                output_filename        = temp_out,
                ignore_chems           = ignore_file_path,
                batch_size             = 50L,
                isname                 = TRUE,
                prev_drug_filepaths    = prev_drug_filepaths,
                restrict_to_raw_names  = chem_list
            )

            # --- merge temp_out into cumulative output_file_path ---
            if (file.exists(temp_out)) {
                if (file.exists(output_file_path)) {
                    agg <- fread(output_file_path, sep="\t", header=TRUE)
                    new_part <- fread(temp_out, sep="\t", header=TRUE)
                    combined <- unique(rbindlist(list(agg, new_part), use.names=TRUE, fill=TRUE))
                } else {
                    combined <- fread(temp_out, sep="\t", header=TRUE)
                }
                fwrite(combined, output_file_path, sep="\t")
                file.remove(temp_out)
                ## Merge succeeded -- record this PSet as done so a later
                ## attempt skips it entirely.
                mark_pset_done(done_tag)
            } else {
                warning("Expected temporary output not found: ", temp_out)
            }

            ## The downloaded PSet is deliberately NOT deleted. It lives in
            ## PSET_CACHE_DIR under /tmp (the bind mount to the host's local/),
            ## so a retried attempt reuses it instead of re-downloading. These
            ## are large -- CTRPv2 is 228 MB and one PSet is 1.5 GB -- and
            ## re-fetching them is what made a single transient Zenodo 504
            ## expensive. Clear local/ between releases to reclaim the space.

        }
    }
}



main<-function(){
	args = commandArgs(trailingOnly=TRUE)
	if(length(args)<2){
	  print('Usage: Rscript 03-createDrugFile.R [datasets] [existing file]')
# 	  exit()
	  }
#	sfile = args[1]
        dsets<-unlist(strsplit(args[1],split=','))
        if(length(args)==2)
            efile=args[2]
        else
            efile=''
       dl1<-getDepMapDrugData(dsets,efile)


}

main()
