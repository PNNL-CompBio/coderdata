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
drugs <- args[3]

# Log in to Synapse
synLogin(authToken = PAT)


# Read the sample mapping CSV and genes.csv
samples_df <- fread(patients)|>
    dplyr::select(improve_sample_id,common_name,model_type)|>
                                        distinct()#"mpnst/synapse_NF-MPNST_samples.csv")

pdx_samps<-subset(samples_df,model_type=='Patient derived xenograft')
org_samps<-subset(samples_df,model_type=='organoid')

##now get the manifest from synapse
manifest<-synapser::synTableQuery("select * from syn53503360")$asDataFrame()|>
                                                             dplyr::rename(common_name='Sample')


##PDX contain list of files
pdx<-manifest|>
    dplyr::select(common_name,PDX_Drug_Data)|>
    left_join(pdx_samps)|>
    distinct()|>
    subset(!is.na(PDX_Drug_Data))


##MTS contain lists of directories
mts<-manifest|>
    dplyr::select(common_name,MicroTissueDrugFolder)|>
    left_join(org_samps)|>
    distinct()|>
    subset(!is.na(MicroTissueDrugFolder))


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
getDrugDataByParent<-function(parid,sampleId){
    qtab<-synTableQuery(paste('select id,name,experimentalCondition,parentId from syn21993642 where parentId=\'',parid,'\''))$asDataFrame()|>
        subset(!is.na(experimentalCondition))|>dplyr::select(id,name,experimentalCondition)
    ##now we need to parse the metadatda table get the info

    res<-do.call(rbind,lapply(qtab$id,function(x){
        sname <- subset(qtab,id==x)
        print(sname)
        sname <-extract_date_hour(sname$name)

        #print(sname)
        data <- fread(synGet(x)$path)|>
            filter(response_type=='percent viability')|>
            mutate(improve_sample_id=sampleId,
                   DOSE=exp(dosage),
                   GROWTH=response /100,
                   source = "NF DATA PORTAL",
                   CELL = improve_sample_id,
                   chem_name = compound_name,
                   study = paste0('MT ',sname$date,' exp'),
                   time = sname$hour) %>%
            select(improve_sample_id,DOSE,GROWTH,source,CELL,chem_name,study,time)

    return(data)
    }))
    return(res)
}

##now loop through manifest to get all the files
mts_fold <- data.table(mts)[,strsplit(as.character(MicroTissueDrugFolder),","), by = .(improve_sample_id,common_name)]

alldrugs<-do.call(rbind,lapply(mts_fold$V1,function(x){
    samp<-subset(mts_fold,V1==x)
    res<-getDrugDataByParent(x,samp$improve_sample_id)
    return(res)
}))

##do the drug matching
drug_df<-fread(drugfile)

alldrugs<-alldrugs|>left_join(drug_df)

fwrite(alldrugs,'curve_data.tsv')


##then run the curve fitting

