# coderdata cell line omics data capture
# pulls down static versions of depmap and sanger cell line omics measurements
# concatenates each into sinle file
library(readr)
library(tidyr)
library(dplyr)
library(rio)
library(httr2)
library(data.table)

## Shared retry/resume helpers (download_file_resumable, with_retries).
source("retry_utils.R")

Sys.setenv(VROOM_CONNECTION_SIZE=100000000)



## Resumable, size-verified download.
##
## Kept under the original name so callers are unchanged, but the transport is
## no longer httr2. req_retry() restarts the whole request instead of resuming
## it, and the files pulled here are the largest in the pipeline -- Sanger ships
## rnaseq_all at ~897MB and WES_pureCN_CNV_genes at ~935MB. On a link that drops
## part-way, restarting a 935MB file just fails again at a similar point: build
## v24 died exactly that way on a 112MB file, breaking 12MB in, three times.
##
## download_file_resumable() (coderbuild/utils/retry_utils.R) uses wget -c, so
## each attempt continues from what is already on disk, and only promotes the
## file once its size matches the server's.
##
## req_verbose() is also gone: full curl logging across ~2GB of transfers buried
## the real errors in the build log.
robust_download_httr2 <- function(url, dest,
                                  max_tries    = 5,
                                  timeout_secs = 1500) {
  download_file_resumable(url, dest, what = basename(dest))
}





# Helper to download a ZIP and extract it safely
download_and_extract_zip_httr2 <- function(url, dest_zip, extract_dir, max_tries = 5, timeout_secs = 1500) {
  robust_download_httr2(url, dest_zip, max_tries = max_tries, timeout_secs = timeout_secs)
  if (!file.exists(dest_zip)) stop(sprintf("Download failed, %s missing", dest_zip))
  tryCatch({
    utils::unzip(dest_zip, exdir = extract_dir)
  }, error = function(e) {
    file.remove(dest_zip)
    stop(sprintf("Failed to unzip %s: %s", dest_zip, e$message))
  })
}



##### DEPMAP FILES
## ---------------------------------------------------------------------------
## DEPMAP FILES ARE NO LONGER PROGRAMMATICALLY RETRIEVABLE.
##
## DepMap put a Cloudflare Turnstile ("verify you are a person") challenge in
## front of https://depmap.org/portal/api/download/files. That endpoint now
## returns an HTML challenge page with HTTP 200, so read_csv() silently parses
## the HTML and every downstream column lookup fails. It cannot be scripted
## around, and DepMap explicitly asks that the portal not be scraped.
##
## We therefore keep a STATIC COPY of the required DepMap files on Synapse, in
## the "DepMap Raw" folder: https://www.synapse.org/Synapse:syn75028495
##
## build_omics.sh runs coderbuild/utils/fetch_depmap.py BEFORE this script,
## which downloads them from Synapse into $DEPMAP_DIR. This script only reads
## from disk -- it never contacts the DepMap portal.
##
## THE SYNAPSE COPY MUST BE MANUALLY UPDATED FOR EVERY NEW DEPMAP RELEASE:
##   1. Download the files from https://depmap.org/portal/data_page/?tab=allData
##      (accept the terms, then use the download button on each file)
##   2. Upload them to syn75028495 as NEW VERSIONS of the existing entities
##   3. Bump DEPMAP_RELEASE in coderbuild/utils/fetch_depmap.py
##   4. Re-run the build with --depmap-ready
##
## Currently pinned to: DepMap Public 26Q1
## ---------------------------------------------------------------------------
## DEPMAP_DIR is exported by build_omics.sh from fetch_depmap.py, so the release
## name is defined in exactly one place. The fallback keeps this script runnable
## standalone for debugging.
DEPMAP_DIR <- Sys.getenv("DEPMAP_DIR", unset = "/opt/depmap_data")

depmap_filenames = list(
  copy_number     = file.path(DEPMAP_DIR, "PortalOmicsCNGeneLog2.csv"),
  transcriptomics = file.path(DEPMAP_DIR, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"),
  mutations       = file.path(DEPMAP_DIR, "OmicsSomaticMutations.csv")
)

missing_depmap <- unlist(depmap_filenames)[!file.exists(unlist(depmap_filenames))]
if (length(missing_depmap) > 0) {
  stop(sprintf(paste0(
    "Required DepMap file(s) not found:\n  %s\n",
    "They should have been downloaded from Synapse by fetch_depmap.py before\n",
    "this script ran. DepMap files can no longer be downloaded from the DepMap\n",
    "portal (Cloudflare challenge), so they are read from the Synapse folder\n",
    "syn75028495. To refresh them for a new release, download from\n",
    "https://depmap.org/portal/data_page/?tab=allData and re-upload to Synapse."),
    paste(missing_depmap, collapse = "\n  ")))
}
message("Using DepMap files from Synapse copy in ", DEPMAP_DIR)

## fetch_depmap.py has already downloaded these into a container-local scratch
## directory, so there is nothing to download here -- we read them in place.
##
## Each file is deleted as soon as it has been read: they are transient build
## inputs totalling ~1.3 GB, and this step is memory- and disk-sensitive (the
## copy_number melt alone is large). Freeing each input before processing the
## next keeps the container's footprint down. The container is `--rm` anyway,
## so this is belt-and-braces, but it matters while the step is running.
depmap_consume <- function(path) {
  if (!is.null(path) && file.exists(path)) {
    unlink(path)
    message("Removed consumed DepMap input: ", path)
  }
}
##### SANGER FILES
sanger_filenames=list(transcriptomics='https://cog.sanger.ac.uk/cmp/download/rnaseq_all_20220624.zip',
               copy_number='https://cog.sanger.ac.uk/cmp/download/WES_pureCN_CNV_genes_latest.csv.gz',
               mutations='https://cog.sanger.ac.uk/cmp/download/mutations_all_20230202.zip')


###### VARIANT SCHEMA HARMONIZATION
variant_schema =list(`3'UTR`=c("3'UTR",'THREE_PRIME_UTR','3prime_UTR_variant','3prime_UTR_ess_splice','3_prime_UTR_variant'),
                     `5'Flank`=c("FIVE_PRIME_FLANK","5'Flank",'upstream','upstream_gene_variant'),
                     `5'UTR`=c("5'UTR",'5prime_UTR_variant','5prime_UTR_variant','5prime_UTR_ess_splice','5_prime_UTR_variant'),
                     Undetermined=c('COULD_NOT_DETERMINE','protein_altering_variant'),
                     De_novo_Start_InFrame=c('DE_NOVO_START_IN_FRAME','De_novo_Start_InFrame'),
                     De_novo_Start_OutOfFrame=c('DE_NOVO_START_OUT_FRAME','De_novo_Start_OutOfFrame'),
                     Frame_Shift_Del=c('FRAME_SHIFT_DEL','Frame_Shift_Del','frameshift','frameshift_variant'),
                     Frame_Shift_Ins=c('FRAME_SHIFT_INS','Frame_Shift_Ins'),
                     IGR=c('IGR','nc_variant','intergenic_variant','downstream_gene_variant'),
                     In_Frame_Del=c('IN_FRAME_DEL','In_Frame_Del','inframe','inframe_deletion'),
                     In_Frame_Ins=c('IN_FRAME_INS','In_Frame_Ins','inframe_insertion'),
                     Intron=c('INTRON','Intron','intronic','intron','intron_variant'),
                     Missense_Mutation=c('Missense_Mutation','MISSENSE','missense','missense_variant'),
                     Nonsense_Mutation=c('Nonsense_Mutation','NONSENSE','nonsense','stop_gained'),
                     Nonstop_Mutation=c('Nonstop_Mutation','NONSTOP','stop_lost'),
                     RNA=c('RNA','non_coding_transcript_exon_variant','non_coding_transcript_variant'),
                     Start_Codon_SNP=c('START_CODON_SNP','Start_Codon_SNP'),
                     Start_Codon_Del=c('Start_Codon_Del','START_CODON_DEL','start_lost'),
                     Start_Codon_Ins=c('Start_Codon_Ins','START_CODON_INS'),
                     Stop_Codon_Del=c('Stop_Codon_Del','stop_lost'),
                     Stop_Codon_Ins=c('Stop_Codon_Ins'),
                     Silent=c('Silent','SILENT','silent','synonymous_variant'),
                     Splice_Site=c('Splice_Site','SPLICE_SITE','splice_region','splice_donor_variant','splice_acceptor_variant','splice_region_variant','splice_polypyrimidine_tract_variant','splice_donor_region_variant','splice_donor_5th_base_variant'),
                     Translation_Start_Site=c('Translation_Start_Site','start_lost'))

depmap_vtab<-do.call('rbind',sapply(names(variant_schema),function(x) cbind(rep(x,length(variant_schema[[x]])),unlist(variant_schema[[x]]))))
colnames(depmap_vtab)<-c('variant_classification','VariantInfo')

sanger_vtab<-do.call('rbind',sapply(names(variant_schema),function(x) cbind(rep(x,length(variant_schema[[x]])),unlist(variant_schema[[x]]))))
colnames(sanger_vtab)<-c('variant_classification','effect')


mirnaFixing<-function(mirlist){

  ##first let's get the prefix off
  trimmed<-lapply(mirlist,function(x)
    stringr::str_replace(x,'hsa-','')|>
      stringr::str_replace_all('-','')|>
      toupper())

  newmap<-data.frame(old=mirlist,gene_symbol=unlist(trimmed))
  return(newmap)
}


#### pull down ad process each of the sanger files
sanger_files<-function(fi,value){

  options(timeout=10000)

#  newres<-lapply(dt,function(value){

 #   fi=sanger_filenames[[value]]
    # fname=paste0('/tmp/sanger_',value,'.csv.gz')
    #print(paste('now reading',fi,'to store as',fname))
    ##now every data type is parsed slightly differently, so we need to change our formatting
    ##and mapping to get it into a unified 3 column schema
    if(value=='copy_number'){
      #read in file
      local_cn <- file.path(tempdir(), "sanger_copy_number.csv.gz")
      robust_download_httr2(fi, local_cn)
      exp_file <- data.table::fread(local_cn,
        select = c("model_id","symbol","gatk_mean_log2_copy_ratio","source","data_type","cn_category"),
        data.table = FALSE)
      # exp_file <- readr::read_csv(fi) ##already in long form <3 <3 <3
      # file.remove(fi)
      smap<-sanger_samples|>
          subset(other_id_source=='Sanger')|>
          subset(other_id%in%exp_file$model_id)|>
          dplyr::select(improve_sample_id,other_id)|>
          distinct()

      gvals <- intersect(exp_file$symbol,genes$gene_symbol)
      gmap<-genes|>
            subset(gene_symbol%in%gvals)|>
            distinct()

      print('wide to long')

      res<-exp_file|>
        dplyr::select(other_id='model_id',gene_symbol='symbol',gatk_mean_log2_copy_ratio,source,data_type,cn_category)|>
        mutate(copy_number=2^gatk_mean_log2_copy_ratio,.keep='all')|>
        distinct()|>
        left_join(gmap)|>
        dplyr::select(other_id,copy_number,entrez_id,Sanger='cn_category')|>
        left_join(smap)|>
          distinct()
       rm(exp_file); gc()
       if (file.exists(local_cn)) file.remove(local_cn)

        print('copy call')

        ##rename SANGER value
        # Amplification -> amp
        # Deletion -> deep del
        # Loss -> het loss
        # Gain -> gain
        # Neutral -> diploid
        #
        res$Sanger=sapply(res$Sanger,function(x) ifelse(x=='Amplification','amp',ifelse(x=='Deletion','deep del',ifelse(x=='Loss','het loss',ifelse(x=='Gain','gain','diploid')))))
      ##calibrate the copy call
      res<-res|> ##deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
        dplyr::mutate(IMPROVE=ifelse(copy_number<0.5210507,'deep del',
                                       ifelse(copy_number<0.7311832,'het loss',
                                              ifelse(copy_number<1.214125,'diploid',
                                                     ifelse(copy_number<1.422233,'gain','amp')))))|>
        dplyr::distinct()|>
        mutate(study='Sanger')

      full<-res|>
        tidyr::pivot_longer(cols=c(IMPROVE,Sanger),
                            names_to='source',
                            values_to='copy_call')
      rm(res); gc()
#      full<-lres

    }else if(value=='methylation'){ ###IF DATA REPRESENT RRBS###
      exp_file <- readr::read_csv(fi)[-c(1:2),]
      #the gene names are not unique and in some weird format, i willt ry to keep both in the metadata
      res = exp_file|>
        tidyr::pivot_longer(cols=c(2:ncol(exp_file)),
                            names_to='gene_region',values_to='methylation',
                            values_transform=list(methylation=as.numeric))|>
        dplyr::distinct()

      res<-res|>
        tidyr::separate(gene_region,into=c('gene_symbol','num','start','end'),sep='_')|>
        dplyr::left_join(genes)|>
        dplyr::distinct()

      colnames(res)[1]<-'other_id'
      vars=c('methylation','start','end')


    }else if(value=='mutations'){ ####IF DATA REPRESENTS MUTATIONS#####
      # res=download.file(fi,'/tmp/tmp.zip')
      # filist<-unzip('/tmp/tmp.zip',exdir='/tmp')
      # fi= "/tmp/mutations_all_20230202.csv"
      zip_path <- file.path(tempdir(), "sanger_mutations_20230202.zip")
      download_and_extract_zip_httr2(fi, zip_path, "/tmp")
      csv_path <- file.path("/tmp", "mutations_all_20230202.csv")
      if (!file.exists(csv_path)) stop("Expected mutations CSV not found after unzip")

    
      exp_file <- readr::read_csv(csv_path) |>
        dplyr::select(gene_symbol, other_id = 'model_id', effect, mutation = 'cdna_mutation', source) |>
        distinct()


      file.remove(csv_path)
      if (file.exists(zip_path)) file.remove(zip_path)

      smap<-sanger_samples|>
        dplyr::select(improve_sample_id,other_id)|>distinct()

      res<-exp_file|>
          left_join(genes)|>
        left_join(smap)|>
          mutate(study='Sanger')|>
          dplyr::select(-c(other_id,gene_symbol))|>
          left_join(as.data.frame(sanger_vtab))

      ##now many variants are missing???
      missing<-res|>
          select(effect,variant_classification)|>
          distinct()|>
          subset(is.na(variant_classification))
      print(missing)

###TODO double check to see if any variants are missing
      res<-res|>dplyr::select(-effect)|>
          subset(!is.na(improve_sample_id))|>
          distinct()

      rm(exp_file)

#      res$variant_classification=unlist(lapply(res$effect,function(x) names(variant_schema)[grep(x,variant_schema)]))
#      res<-res|>dplyr::select(-effect)
      #write_csv(res,file=fname)
      #rm(res)
                                        #return(fi)
      print(head(res))
      return(res)
    }else if(value=='transcriptomics'){ #if gene expression
      # res=download.file(fi,'/tmp/tmp.zip')
      # filist<-unzip('/tmp/tmp.zip',exdir='/tmp')

      zip_path <- "/tmp/rnaseq_all_20220624.zip"
      download_and_extract_zip_httr2(fi, zip_path, "/tmp")
      csv_path <- file.path("/tmp", "rnaseq_tpm_20220624.csv")
      if (!file.exists(csv_path)) stop("Expected transcriptomics CSV not found after unzip")
      exp_file <- readr::read_csv(csv_path)

      ##the rows have metadata
      samps<-t(exp_file[1:3,])
      colnames(samps)<-samps[1,]
      samps<-samps[-c(1:2),]|>as.data.frame()|>
        tibble::rownames_to_column('other_id')|>
        left_join(sanger_samples)|>
        dplyr::rename(source='data_source',study='dataset_name')

      missing<-subset(samps,is.na(improve_sample_id))|>
        dplyr::select(-c(other_id,improve_sample_id))|>
        dplyr::rename(other_id='model_name')

      ##the data starts at row 5
      dat<-exp_file[-c(1:4),-1]
      colnames(dat)[1]<-'gene_symbol'

      rm(exp_file); gc()

      ddat<-apply(dat,2,unlist)|>
          as.data.frame()|>
#      ddat<-ddat|>
        left_join(genes)|>
        dplyr::select(-gene_symbol)

      ## gc() here matters: without it R still holds dat when pivot_longer
      ## allocates the long frame, which is the peak of the whole step.
      ## Measured on a 20k gene x 800 sample matrix (16M output rows):
      ## peak RSS 1,351MB -> 1,219MB, run time 8.2s -> 8.1s. The other
      ## branches in this function already pair rm() with gc().
      rm(dat); gc()
      ddat$entrez_id<-as.numeric(ddat$entrez_id)
      res = tidyr::pivot_longer(data=as.data.frame(ddat),cols=c(1:(ncol(ddat)-1)),
                                names_to='other_id',values_to='transcriptomics',
                                values_transform=list(expression=as.numeric))|>
          distinct()

      rm(ddat); gc()
      smap<-samps|>
        dplyr::select(improve_sample_id,other_id,study,source)|>distinct()

      full<-res|>
        left_join(smap)
      rm(res)

      file.remove(csv_path)
      if (file.exists(zip_path)) file.remove(zip_path)


    }else if(value=='miRNA'){ #if mirna expression
      exp_file <- readr::read_csv(fi)

      res = tidyr::pivot_longer(data=exp_file,cols=c(2:ncol(exp_file)),
                                names_to='gene_symbol',values_to='miRNA',values_transform=list(mirnas=as.numeric))
      file.remove(fi)
      gmod<-genes
      gmod$gene_symbol<-tolower(gmod$gene_symbol)

      res$gene_symbol<-tolower(res$gene_symbol)

      res<-res|>
        dplyr::left_join(gmod)|>
        dplyr::distinct()
      #dplyr::mutate(gene_symbol=tolower(gene_symbol))

      missed<-subset(res,is.na(entrez_id))
      notmissed<-subset(res,!is.na(entrez_id))
      print(paste("matched",length(unique(notmissed$gene_symbol)),'miRNAs and missed',length(unique(missed$gene_symbol))))

      ##mirfix
      mirmap<-mirnaFixing(unique(missed$gene_symbol))
      fixed<-missed|>
        dplyr::rename(old='gene_symbol')|>
        dplyr::select(-entrez_id)|>
        left_join(mirmap)|>
        left_join(genes)
      notmissed<-fixed|>
        subset(!is.na(entrez_id))|>
        dplyr::select(-old)|>
        rbind(notmissed)

      missed<-subset(fixed,is.na(entrez_id))
      print(paste("after second pass, matched",length(unique(notmissed$gene_symbol)),'miRNAs and missed',length(unique(missed$gene_symbol))))

      colnames(res)[1]<-'other_id'
      vars=c('miRNA')
      full<-res

    }else if(value=='proteomics'){
      # res=download.file(fi,'/tmp/tmp.zip')
      # filist<-unzip('/tmp/tmp.zip',exdir='/tmp')

      zip_path <- "/tmp/Proteomics_20221214.zip"
      download_and_extract_zip_httr2(fi, zip_path, "/tmp")
      tsv_path <- file.path("/tmp", "Protein_matrix_averaged_zscore_20221214.tsv")
      if (!file.exists(tsv_path)) stop("Expected proteomics TSV not found after unzip")

      exp_file <- readr::read_tsv(tsv_path,skip=1)[-1,-1]
      colnames(exp_file)[1]<-'other_id'

      file.remove(tsv_path)
      if (file.exists(zip_path)) file.remove(zip_path)

      smap<-sanger_samples|>
        dplyr::select(improve_sample_id,other_id)|>distinct()

      res<-exp_file|>
        tidyr::pivot_longer(cols=(-c(other_id)),names_to='gene_symbol',values_to='proteomics')|>
        subset(!is.na(proteomics))|>
        left_join(genes)|>
        dplyr::select(-gene_symbol)|>
        left_join(smap)|>
          mutate(source='Sanger',study='Sanger')
      rm(exp_file)
      full<-res
      rm(res)

    }

    ##do the last join with samples
    #full<-res#|>
    #  left_join(samples)
    missed<-full|>subset(is.na(improve_sample_id))|>
      dplyr::select(improve_sample_id,other_id)|>
      distinct()
    print(paste('missing',nrow(missed),' sample identifiers'))
    print(missed)

    full<-full|>
        subset(!is.na(improve_sample_id))|>
        subset(!is.na(entrez_id))|>
      dplyr::select(-other_id)

    #write_csv(full,file=gzfile(fname))
#      rm(full)
      print(paste('Sanger',value,':'))
      print(head(full))
      return(full)
  #  })
  #    names(newres)<-names(sanger_filenames)

}


#### pull dodwn and process each depmap file per the 2023 q4 schema
depmap_files<-function(fi,value){

   ##runs through entire list of omcis types and files and returns  a list
#  newres<-lapply(values,function(value){

#    fi=depmap_filenames[[value]]
#    fname=paste0('/tmp/depmap_',value,'.csv.gz')
#    print(paste('now reading',fi,'to store as',fname))
    ##now every data type is parsed slightly differently, so we need to change our formatting
    ##and mapping to get it into a unified 3 column schema
    if(value=='copy_number'){
      local_path <- fi
      ## Use data.table::melt for memory-efficient wide->long conversion
      exp_dt <- data.table::fread(local_path, data.table = TRUE)
      depmap_consume(local_path)  # read into memory; drop the on-disk copy
      id_col <- colnames(exp_dt)[1]
      data.table::setnames(exp_dt, id_col, "other_id")
      gene_cols <- setdiff(colnames(exp_dt), "other_id")

      print('Long to wide (data.table melt)')
      res <- data.table::melt(exp_dt, id.vars = "other_id",
                              measure.vars = gene_cols,
                              variable.name = "gene_entrez",
                              value.name = "copy_number",
                              variable.factor = FALSE)
      rm(exp_dt, gene_cols); gc()
      ## (the DepMap input was already removed by depmap_consume above)

      res <- as.data.frame(res) |>
        dplyr::mutate(
          copy_number = suppressWarnings(as.numeric(copy_number)),
          ## PortalOmicsCNGeneLog2 stores log2(absolute copies); convert to ratio
          copy_number = 2^copy_number / 2
        ) |>
        dplyr::distinct()

      print('String manipulations')
      res <- res |>
        dplyr::mutate(
          gene_symbol = trimws(stringr::str_extract(gene_entrez, "^[^(]+")),
          entrez_id   = stringr::str_extract(gene_entrez, "(?<=\\()\\d+(?=\\))")
        ) |>
        dplyr::select(-gene_entrez)

      print('join with gene')
      res<-res|>
          dplyr::select(-entrez_id)|>
          left_join(genes)|>
          dplyr::select(other_id,entrez_id,copy_number)|>
          subset(!is.na(entrez_id))|>
          distinct()
      ##these are messing things up
    #  res$entrez_id<-stringr::str_replace(res$entrez_id,'\\)','')
    #  res$entrez_id<-stringr::str_replace(res$entrez_id,'\\(','')

      print('Adding copy calls')
      res<-res|> ##deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
        dplyr::mutate(copy_call=ifelse(copy_number<0.5210507,'deep del',
                                       ifelse(copy_number<0.7311832,'het loss',
                                              ifelse(copy_number<1.214125,'diploid',
                                                     ifelse(copy_number<1.422233,'gain','amp')))))|>
          dplyr::distinct()|>
          subset(!is.na(copy_number))


      vars=c('copy_number','copy_call')


    }else if(value=='methylation'){ ###IF DATA REPRESENT RRBS###
      exp_file <- readr::read_csv(fi)[-c(1:2),]
      #the gene names are not unique and in some weird format, i willt ry to keep both in the metadata
      res = exp_file|>
        tidyr::pivot_longer(cols=c(2:ncol(exp_file)),
                                names_to='gene_region',values_to='methylation',
                            values_transform=list(methylation=as.numeric))|>
        dplyr::distinct()

      res<-res|>
        tidyr::separate(gene_region,into=c('gene_symbol','num','start','end'),sep='_')|>
        dplyr::left_join(genes)|>
          dplyr::distinct()|>
          subset(!is.na(entrez_id))

      colnames(res)[1]<-'other_id'
      vars=c('methylation','start','end')


      }else if(value=='mutations'){ ####IF DATA REPRESENTS MUTATIONS#####

        local_mut <- fi

        exp_file <- readr::read_csv(local_mut)|>
          dplyr::filter(IsDefaultEntryForModel %in% c(TRUE, "Yes"))|>
          dplyr::select(EntrezGeneID,HgncName,other_id='ModelID',VariantInfo,mutation='DNAChange')|>
          distinct()|>
          dplyr::mutate(VariantInfo = sub("&.*", "", VariantInfo))
        depmap_consume(local_mut)  # read into memory; drop the on-disk copy

        res<-exp_file|>
          mutate(entrez_id=as.numeric(EntrezGeneID))|>
            filter(entrez_id %in% genes$entrez_id) |>
            left_join(as.data.frame(depmap_vtab))

              ##now many variants are missing???
        missing<-res|>
            select(VariantInfo,variant_classification)|>
            distinct()|>
            subset(is.na(variant_classification))
        print(missing)

        res<-res|>
            dplyr::select(-c(EntrezGeneID,VariantInfo))|>
            distinct()|>
          subset(!is.na(entrez_id)) ##removes thos with unknonw entrez

        rm(exp_file)

        smap<-depmap_samples|>
            dplyr::select(improve_sample_id,other_id)|>
            distinct()|>
            subset(other_id%in%res$other_id)

        full<-res|>  ###since we're already in ENTREZ we skip the mapping below
          dplyr::left_join(smap)|>
          #dplyr::rename(entrez_id=Entrez_id,mutations=Genome_Change,variant_classification=Variant_Classification)|>
          dplyr::select(entrez_id,improve_sample_id,mutation,variant_classification)|>
           dplyr::mutate(source='Broad',study='DepMap')|>
          dplyr::distinct()

          #write_csv(full,file=fname)
                                        #return(fi)
        print(head(full))
        return(full)
      }else if(value=='transcriptomics'){ #if gene expression
        # exp_file <- readr::read_csv(fi)
        local_tx <- fi
        ## 26Q1 format has metadata columns (SequencingID, ModelConditionID,
        ## IsDefaultEntryForMC, IsDefaultEntryForModel) before the gene columns.
        ## Filter to one canonical run per model, then keep only ModelID + gene cols.
        exp_file <- readr::read_csv(local_tx, show_col_types = FALSE)|>
          dplyr::filter(IsDefaultEntryForModel %in% c(TRUE, "Yes"))|>
          dplyr::rename(other_id = ModelID)|>
          dplyr::select(other_id, dplyr::matches("\\(\\d+\\)"))
        depmap_consume(local_tx)  # read into memory; drop the on-disk copy

        print("wide to long")
        res = tidyr::pivot_longer(data=exp_file,cols=c(2:ncol(exp_file)),
                                  names_to='gene_entrez',values_to='transcriptomics',
                                  values_transform=list(transcriptomics=as.numeric))|>
                                  dplyr::mutate(transcriptomics = 2^transcriptomics - 1)
        colnames(res)[1]<-'other_id'

        print('fixing gene names')
        res <- res |>
          dplyr::mutate(
            gene_symbol = trimws(stringr::str_extract(gene_entrez, "^[^(]+")),
            entrez_par  = stringr::str_extract(gene_entrez, "(?<=\\()\\d+(?=\\))")
          ) |>
          dplyr::select(-gene_entrez)

      print('join with gene')
      res<-res|>
          dplyr::select(-entrez_par)|>
          left_join(genes)|>
          dplyr::select(other_id,entrez_id,transcriptomics)|>
          subset(!is.na(entrez_id))|>
          distinct()

    #mutate(entrez_id=stringr::str_replace_all(entrez_par,'\\)|\\(',''))|>
    #      dplyr::select(-c(entrez_par,gene))|>
    #      distinct()
        rm(exp_file)
        vars=c('transcriptomics')

      }else if(value=='miRNA'){ #if mirna expression
        exp_file <- readr::read_csv(fi)

        res = tidyr::pivot_longer(data=exp_file,cols=c(2:ncol(exp_file)),
                                  names_to='gene_symbol',values_to='miRNA',values_transform=list(mirnas=as.numeric))

        rm(exp_file)
        gmod<-genes
        gmod$gene_symbol<-tolower(gmod$gene_symbol)

        res$gene_symbol<-tolower(res$gene_symbol)

        res<-res|>
          dplyr::left_join(gmod)|>
          dplyr::distinct()
          #dplyr::mutate(gene_symbol=tolower(gene_symbol))

        missed<-subset(res,is.na(entrez_id))
        notmissed<-subset(res,!is.na(entrez_id))
        print(paste("matched",length(unique(notmissed$gene_symbol)),'miRNAs and missed',length(unique(missed$gene_symbol))))

        ##mirfix
        mirmap<-mirnaFixing(unique(missed$gene_symbol))
        fixed<-missed|>
          dplyr::rename(old='gene_symbol')|>
          dplyr::select(-entrez_id)|>
          left_join(mirmap)|>
          left_join(genes)
        notmissed<-fixed|>
          subset(!is.na(entrez_id))|>
          dplyr::select(-old)|>
          rbind(notmissed)

        missed<-subset(fixed,is.na(entrez_id))
        print(paste("after second pass, matched",length(unique(notmissed$gene_symbol)),'miRNAs and missed',length(unique(missed$gene_symbol))))

        colnames(res)[1]<-'other_id'
        vars=c('miRNA')

      }

      ##do the last join with samples
      smap<-depmap_samples|>
          subset(other_id%in%res$other_id)|>
          distinct()

      print('joining with samples')
      full<-res|>
          dplyr::left_join(smap)
      rm(res)

      missed<-full|>subset(is.na(improve_sample_id))|>
          dplyr::select(improve_sample_id,other_id)|>
          distinct()
      print(paste('missing',nrow(missed),'identifiers'))
      print(missed)

      full<-full|>dplyr::select(c('entrez_id','improve_sample_id',vars))|>
          subset(entrez_id%in%genes$entrez_id)|>
          subset(!is.na(improve_sample_id))|>
          dplyr::distinct()|>
          dplyr::mutate(source='Broad',study='DepMap')

      #write_csv(full,file=gzfile(fname))
      #rm(full)
      print(paste('Depmap',value,':'))
      print(head(full))
      return(full)


}


main<-function(){
    args = commandArgs(trailingOnly=TRUE)
    if(length(args)!=2){
        stop('Usage: Rscript 02-depmap-sanger-omics.R [genefile] [samplefile]')
    }
    gfile = args[1]
    sfile = args[2]
###first reqad in all gene information so we can map appropriately
    allgenes = read_csv(gfile)

    genes <<- allgenes|>
        dplyr::select(gene_symbol,entrez_id)|>
        dplyr::distinct()

    ##here are the improve sample id indices
    depmap_samples <<- read_csv(sfile,
                               quote='"')|>
        dplyr::select(other_id,improve_sample_id)|>
        unique()

    sanger_samples <<- read_csv(sfile,
                                quote='"')|>
        dplyr::select(other_id,improve_sample_id,other_id_source)|>
        unique()

    alltypes<-c('mutations','transcriptomics','copy_number')

    lapply(alltypes,function(dt){
        print(dt)
        temps<-sanger_files(sanger_filenames[[dt]],dt)|>tidyr::drop_na()|>dplyr::distinct()
        readr::write_csv(temps,file=paste0('/tmp/sanger_',dt,'.csv.gz'))
        rm(temps); gc()
        tempd<-depmap_files(depmap_filenames[[dt]],dt)|>tidyr::drop_na()|>dplyr::distinct()
        readr::write_csv(tempd,file=paste0('/tmp/broad_',dt,'.csv.gz'))
        rm(tempd); gc()
    })

}

main()

