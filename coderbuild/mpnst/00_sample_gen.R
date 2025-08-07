# This script generate a new sample table based on previous dataset's sample file (taking the max improve_sample_id)
# Load required libraries
library(data.table)
library(synapser)
library(dplyr)

##adding a command line argument
args = commandArgs(trailingOnly=TRUE)
if(length(args) > 1 ){
    stop("Up to one argument is allowed. This is the filepath to the previously run samples file.")
}

if (length(args) == 0 || is.na(args[1]) || args[1] == "" || !file.exists(args[1])) {
    orig_samples <- ""
} else {
    orig_samples <- fread(args[1])
}

# Check if Synapse token is available from the environment
synapse_token <- Sys.getenv("SYNAPSE_AUTH_TOKEN")
if (synapse_token == "") {
    stop("Error: SYNAPSE_AUTH_TOKEN environment variable is not set.")
}

synapser::synLogin(authToken=synapse_token)
manifest<-synapser::synTableQuery("select * from syn53503360")$asDataFrame()|>
                                                             as.data.frame()

#Drop contaminated sample JH-2-009
manifest <- manifest %>% 
  filter(Sample != "JH-2-009")


###sample file has a strict schema
## - improve_sample_id
## - other_id
## - other_id_source
## - common_name:
## - cancer_type:
## - other_names:
## - species:
## -  model_type:

##first create samples for the original tumors
tumorTable<-manifest|>
    dplyr::select(common_name='Sample')|>
    dplyr::mutate(other_id_source='NF Data Portal',other_names='',cancer_type="Malignant peripheral nerve sheath tumor",species='Homo sapiens (Human)',model_type='tumor')|>
    tidyr::unite(col='other_id',c('common_name','model_type'),sep=' ',remove=FALSE)

##then create samples for the PDX
sampTable<-manifest|>
    dplyr::select(common_name='Sample',MicroTissueDrugFolder)|>
    dplyr::mutate(other_id_source='NF Data Portal',other_names='',cancer_type="Malignant peripheral nerve sheath tumor",species='Homo sapiens (Human)',model_type='patient derived xenograft')|>
    tidyr::unite(col='other_id',c('common_name','model_type'),sep=' ',remove=FALSE)


##third, generate a sample for the MTs if they were generated
pdxmt<-subset(sampTable,!is.na(MicroTissueDrugFolder))
pdxmt$model_type=rep('xenograft derived organoid',nrow(pdxmt))
print(pdxmt)

main<-rbind(sampTable,pdxmt)|>
    dplyr::select(-MicroTissueDrugFolder)|>
    rbind(tumorTable)

# If there is no previous samples file - start at 1, else, continue where the previous one left off.
if (identical(orig_samples, "")) {
    max_id <- 1  
} else {
    max_id <- max(orig_samples$improve_sample_id, na.rm = TRUE)
}

main$improve_sample_id <- seq(from = max_id + 1, length.out = nrow(main))

fwrite(main,'/tmp/mpnst_samples.csv')
