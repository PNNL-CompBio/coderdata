##this file creates a sample database table from Priyanka's depmap file
##and cellosaurus. it ONLY handles cell lines. ADditional samples must be added
##subsequent to running this code.

##author sara gosline
## email sara.gosline@pnnl.gov

library(curl)
library(dplyr)
library(readr)
library(XML)

## Retry helpers for flaky third-party downloads
source("retry_utils.R")

##the only thing that Priyanka has here is TRP identifiers, so collecting those
# tab<-read.table('DepMap_Argonne_Mapping.csv',sep=',',header=T)%>%
#   dplyr::select(Argonne_ID,DepMap_ID)%>%
#   distinct()%>%
#   tidyr::separate(Argonne_ID,into=c('id_source','other_id'),sep='\\.')%>%
#   subset(id_source=='CTRP')%>%
#   dplyr::rename(CTRP='other_id',ModelID='DepMap_ID')%>%
#   dplyr::select(-id_source)
##we miss these
###[1] "Panc-05-04"              "L3-3"                    "Panc-02-03"              "Panc-10-05"              "G-292-clone-A141B1"      "Hs-729"                  "Panc-03-27"              "Panc-04-03"              "Hs-688-A-T"
###[10] "Panc-08-13"              "PE-CA-PJ34-clone-C12"    "Hep-3B2-1-7"             "PE-CA-PJ41-clone-D2"     "PE-CA-PJ49"              "Ishikawa-Heraklio-02-ER"

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
## build_samples.sh runs coderbuild/utils/fetch_depmap.py BEFORE this script,
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
## Files required by coderdata (4 total, ~1.29 GB):
##   Model.csv                                          (samples - this script)
##   OmicsSomaticMutations.csv                          (omics)
##   OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv (omics)
##   PortalOmicsCNGeneLog2.csv                          (omics)
## ---------------------------------------------------------------------------
## DEPMAP_DIR is exported by build_samples.sh from fetch_depmap.py, so the
## release name is defined in exactly one place. The fallback keeps this script
## runnable standalone for debugging.
DEPMAP_DIR <- Sys.getenv("DEPMAP_DIR", unset = "/opt/depmap_data")
depmap_model_file <- file.path(DEPMAP_DIR, "Model.csv")

if (!file.exists(depmap_model_file)) {
  stop(sprintf(paste0(
    "Required DepMap file not found: %s\n",
    "It should have been downloaded from Synapse by fetch_depmap.py before\n",
    "this script ran. DepMap files can no longer be downloaded from the DepMap\n",
    "portal (Cloudflare challenge), so they are read from the Synapse folder\n",
    "syn75028495. To refresh them for a new release, download from\n",
    "https://depmap.org/portal/data_page/?tab=allData and re-upload to Synapse."),
    depmap_model_file))
}

message("Using DepMap file from Synapse copy: ", depmap_model_file)
depmap_models<-readr::read_csv(depmap_model_file)#|>
#  dplyr::rename(DepMap_ID='ModelID')
## Transient container-local input -- drop it once it is in memory.
unlink(depmap_model_file)

## cog.sanger.ac.uk intermittently drops the connection mid-transfer
## (OpenSSL "unexpected eof"); this failed twice before succeeding on the
## third attempt on 2026-08-30.
sanger_models<-with_retries(
  function() readr::read_csv("https://cog.sanger.ac.uk/cmp/download/model_list_20230923.csv",
                             show_col_types = FALSE),
  what = "download Sanger model_list")

print(paste("Downloaded",nrow(depmap_models),'dep map identifiers and',nrow(sanger_models),'sanger models'))

##query for cellosaurus automagically to get loadest version
url='https://ftp.expasy.org/databases/cellosaurus/cellosaurus.xml'
## 672 MB download; a partial file would poison every later run because the
## existence check below would treat it as complete. Download to a temp path
## and only move it into place once it has fully arrived.
if(!file.exists('/tmp/cell.xml')) {
  with_retries(function() {
    tmp_cello <- '/tmp/cell.xml.partial'
    curl_download(url, tmp_cello, quiet=TRUE)
    if (file.info(tmp_cello)$size < 1e8)
      stop(sprintf('cellosaurus download truncated (%s bytes)', file.info(tmp_cello)$size))
    file.rename(tmp_cello, '/tmp/cell.xml')
  }, what = "download cellosaurus.xml")
}
cello<-XML::xmlParse('/tmp/cell.xml')
cdf<-XML::xmlToList(cello)

print('Got all cellosaurus ids')
### next we parse through cellosaurus to get as many samples as we deem relevant
##ok, this command seems to have gotten file in appropriate state
cell.lines<-lapply(cdf$`cell-line-list`, function(x) unlist(x))

##now we need toe xtract columns
options(show.error.messages=TRUE)
full.res<-do.call(rbind,lapply(cell.lines,function(x){
    ##create a data frame for each cell lines
    x<-unlist(x)
                                        #should only be one acession
    acc<-x[grep('accession.text',names(x),fixed=T)]%>%unlist()

    cn<-x[grep('name-list.name.text',names(x),fixed=T)]%>%unlist()
                                        #these will fail if no key found
    spec<-x[grep("species-list.xref.label",names(x),fixed=T)]%>%unlist()
    #print(acc)
    #print(cn)
    #print(spec)
                                        #dis<-x[grep("disease-list.cv-term.text",names(x),fixed=T)]%>%unlist()
    data.frame(accession=cn,
               RRID=rep(acc,length(cn)),
               species=rep(spec,length(cn)))
                                        #          disease=rep(dis,length(cn)))
}))%>%
    subset(species=='Homo sapiens (Human)')


print(paste('Got',nrow(full.res),'human cellosaurus samples'))
######
#now we join the sagner table, depmap table, the cellosaurus table, and the CTRP identifiers

allmod<-sanger_models|>
  dplyr::rename(ModelID='BROAD_ID')|>
  left_join(depmap_models)|>
  dplyr::select(Sanger_ID='model_id',Depmap_ID='ModelID',PatientID,CellLineName,StrippedCellLineName,COSMIC_ID,CCLE_ID,RRID,cancer_type,ModelType)|>
  subset(!is.na(RRID))|> ##only take those we can map to cellosaurus
  distinct()

##we missed a handful here, so let's add back
missedmod<-subset(depmap_models,ModelID%in%setdiff(depmap_models$ModelID,allmod$Depmap_ID))|>
  subset(!is.na(RRID))|>
  dplyr::select(Depmap_ID='ModelID',CCLE_ID='CCLEName',COSMIC_ID='COSMICID',CellLineName,StrippedCellLineName,Sanger_ID='SangerModelID',RRID,cancer_type='OncotreeSubtype',ModelType)|>
  mutate(PatientID=NA)

allmod<-rbind(allmod,missedmod)
joined.df<-allmod%>%
  left_join(full.res)

##now lengethn the table to have the appopriate columns
full.df<-joined.df%>%
  dplyr::rename(DepMap='Depmap_ID',Sanger='Sanger_ID',##old file formatSanger='SangerModelID'
                CCLE='CCLE_ID',#old file CCLE='CCLE_ID'
                Cellosaurus='RRID',
                 common_name='CellLineName',
               # cancer_type='OncotreeSubtype',
                other_names='accession')%>%
  mutate(COSMIC=as.character(COSMIC_ID))|>#,WTSI=as.character(WTSIMasterCellID))%>%
  dplyr::select(-c(COSMIC_ID))|>#,WTSIMasterCellID))%>%
  distinct()


#we add idnetifiers for everything that has a cellosaurus id and those that dont
#has_id<-subset(full.df,Cellosaurus!="")
#no_id<-subset(full.df,is.na(Cellosaurus))

samp_ids<-data.frame(Cellosaurus=unique(full.df$Cellosaurus))
samp_ids$improve_sample_id<-seq(1,length(unique(full.df$Cellosaurus)))

#now we need to add in the missing ones
#if(nrow(no_id)>0){
#  extra_ids<-data.frame(DepMap=no_id$DepMap)
#  extra_ids$improve_sample_id<-seq(max(samp_ids$improve_sample_id)+1,
#                                  max(samp_ids$improve_sample_id)+nrow(extra_ids))
#
#}

#full.df<-rbind(left_join(has_id,samp_ids),
#               left_join(no_id,extra_ids))

long.df<-full.df%>%
  left_join(samp_ids)|>
  dplyr::select(-c(StrippedCellLineName))|>
  tidyr::pivot_longer(cols=c(PatientID,DepMap,Sanger,CCLE,COSMIC,Cellosaurus),names_to='other_id_source',
                      values_to='other_id')%>%
  mutate(model_type=dplyr::case_when(
    ModelType=='Organoid' ~ 'patient derived organoid',
    TRUE                  ~ 'cell line'
  ))%>%
  dplyr::select(-ModelType)%>%
  subset(!is.na(other_id))%>%
  subset(other_id!="")


readr::write_csv(long.df,'/tmp/broad_sanger_samples.csv',quote='needed')

