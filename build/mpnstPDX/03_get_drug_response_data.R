# Load required libraries
library(data.table)
# library(biomaRt)# biomart issues still exist
library(synapser)
library(dplyr)
library(stringr)
# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if a token was provided
if (length(args) == 0) {
  stop("No token or sample file provided. Usage: Rscript my_script.R <PAT> [samples] [drugs]", call. = FALSE)
}

# Set your personal access token
PAT <- args[1]
patients <- args[2]
drugfile <- args[3]

# Log in to Synapse
synLogin(authToken = PAT)


# Read the sample mapping CSV and genes.csv
samples_df <- fread(patients)|>
    dplyr::select(improve_sample_id,common_name,model_type)|>
    distinct()#"mpnst/synapse_NF-MPNST_samples.csv")
print(head(samples_df))

pdx_samps<-subset(samples_df,model_type=='patient derived xenograft')
org_samps<-subset(samples_df,model_type=='organoid')

##now get the manifest from synapse
manifest<-synapser::synTableQuery("select * from syn53503360")$asDataFrame()|>
                                                             as.data.frame()|>
                                                             dplyr::rename(common_name='Sample')


##PDX contain list of files
pdx<-manifest|>
    dplyr::select(common_name,PDX_Drug_Data)|>
    subset(!PDX_Drug_Data%in%c("NA",NA))|>
    left_join(pdx_samps)|>
    distinct()

print(pdx)


# Modify the extract_date_hour function to return a named vector
extract_date_hour <- function(experiment_id) {
  pattern <- "(\\d{6})_?(\\d{2,3})?"
  matches <- str_match(experiment_id, pattern)
  date <- matches[, 2]
  hour <- matches[, 3]
  date[is.na(date)] <- NA  # Replace with NA instead of blank
  hour[is.na(hour)] <- 48  # Replace with 48 instead of blank (default)
  return(list(date = date, hour = hour))
}



##define functions

##first function to get children from parentId

##now loop through manifest to get all the files
#mts_fold <- data.table(mts)[,strsplit(as.character(MicroTissueDrugFolder),","), by = .(improve_sample_id,common_name)]



##do the drug matching
drug_df<-fread(drugfile)|>
    dplyr::select('improve_drug_id','chem_name')|>
    distinct()

##update drug name PD901 since it's mussing
##now loop through manifest to get all the files
pdx_fold <- data.table(pdx)[,strsplit(as.character(PDX_Drug_Data),","), by = .(common_name)]|>
    dplyr::rename(id='V1')|>
    subset(!is.na(id))

pdx_meta<-do.call(rbind,lapply(pdx_fold$id, function(x) synapser::synGetAnnotations(x)|>
                                           as.data.frame()|>
                                          dplyr::select('experimentalCondition')|>
                                          dplyr::mutate(id=x)))|>left_join(pdx_fold)|>
   # tidyr::separate_rows(experimentalCondition,sep=';')|>
  #  mutate(chem_name=tolower(experimentalCondition))|>
   # left_join(drug_df)|>
    left_join(pdx_samps)|>
    dplyr::select(improve_sample_id,id)|>
    distinct()|>
    subset(!is.na(id))
pdx_meta$parentId=unlist(lapply(pdx_meta$id,function(x) synGet(x)$parentId))

##the older pdx data is in separate files. the newer is not.
#we need to reformat the older to look like the newer
oldfolders=c('syn22018363','syn22024460','syn22024428','syn22024429','syn22024437','syn22024438')

old_meta<-subset(pdx_meta,parentId%in%oldfolders)
old_data<-do.call(rbind,lapply(unique(old_meta$parentId),function(x){
    ids<-subset(old_meta,parentId==x)|>
        subset(!is.na(id))

  do.call(rbind,lapply(ids$id,function(y){
      tab<-readr::read_csv(synapser::synGet(y)$path)
      print(head(tab))
      tab<-dplyr::select(tab,c('specimen_id','compound_name','dose','dose_unit',
                             'experimental_time_point','experimental_time_point_unit',
                             'assay_type','assay_value','assay_units'))|>
        mutate(id=x)|>
        mutate(chem_name=tolower(compound_name))

   #   tab$single_or_combo=sapply(tab$chem_name,function(z) ifelse(length(grep('\\+',z))>0,'combo','single'))
      tab$chem_name=gsub('n/a','control',tab$chem_name)|>
        tidyr::replace_na('control')

      tab$chem_name=sapply(tab$chem_name,function(z) ifelse(z=='doxorubinsin','doxorubicin',z))
    #  tab<-tab|>left_join(drug_df)
      #print(head(tab))
      return(tab)
       }))
}))|>
  left_join(unique(select(old_meta,id=parentId,improve_sample_id)))|>
  dplyr::select(experiment=id,model_id=improve_sample_id,specimen_id,treatment=chem_name,time=experimental_time_point,volume=assay_value)|>distinct()



new_meta<-subset(pdx_meta,!parentId%in%oldfolders)

##now combine each of the old pdx files into single files
#each file has all experiments in it
new_data<-do.call(rbind,lapply(unique(new_meta$id), function(x){
    fpath=synapser::synGet(x)$path
    if(length(grep('xls',fpath))>0){
        tab<-readxl::read_xlsx(fpath)
    }else{
        tab<-readr::read_csv(fpath)
    }
    print(head(tab))
    tab<-dplyr::select(tab,c('specimen_id','compound_name','dose','dose_unit',
                             'experimental_time_point','experimental_time_point_unit',
                             'assay_type','assay_value','assay_units'))|>
        mutate(id=x)

   # tab$single_or_combo=sapply(tab$compound_name,function(x) ifelse(length(grep('\\+',x))>0,'combo','single'))
    tab$compound_name=gsub('N/A','control',tab$compound_name)|>tidyr::replace_na('control')
    tab<-tab|>
      mutate(chem_name=tolower(compound_name))#|>
   #   left_join(drug_df)
    #print(head(tab))
    return(tab)}))|>
    left_join(pdx_meta)|>
    dplyr::select(experiment=id,model_id=improve_sample_id,specimen_id,treatment=chem_name,time=experimental_time_point,volume=assay_value)|>distinct()

##maybe tweak the data frame a bit depending on curve fitting script
pdx_data<-rbind(old_data,new_data)

#single_pdx<-subset(pdx_data,single_or_combo=='single')
#combo_pdx<-subset(pdx_data,single_or_combo=='combo')
#print(head(pdx_data))
fwrite(pdx_data,'/tmp/curve_data.tsv',sep='\t')

##TODO: create new curve fitting script in python
pycmd = '/opt/venv/bin/python fit_pdx_curve.py --input /tmp/curve_data.tsv --output /tmp/experiments'
print('running curve fitting')
#system(pycmd)

##now read in data again, separate out by single/combo, then map to drug id

##mmve file name
#file.rename('/tmp/experiments.0','/tmp/mpnstpdx_experiments.tsv')


