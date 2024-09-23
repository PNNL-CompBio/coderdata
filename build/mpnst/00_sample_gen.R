# This script generate a new sample table based on pervious beatAML improved sample ID
# It will take the maximum value of beatAML improved sample ID and continue from ID count from there
# Load required libraries
library(data.table)
library(synapser)
library(dplyr)

##adding a command line argument
args = commandArgs(trailingOnly=TRUE)
if(length(args) > 1 ){
    stop("Up to one argument is allowed. This is the filepath to the previously run samples file.")
}


if (file.size(args[1]) == 0) {
    orig_samples <- ""
} else {
    orig_samples <- fread(args[1])
}


# Check if Synapse token is available from the environment
synapse_token <- Sys.getenv("SYNAPSE_AUTH_TOKEN")
if (synapse_token == "") {
    stop("Error: SYNAPSE_AUTH_TOKEN environment variable is not set.")
}

synapser::synLogin(authToken=args[2])
manifest<-synapser::synTableQuery("select * from syn53503360")$asDataFrame()|>
                                                             as.data.frame()


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
    dplyr::mutate(other_id_source='NF Data Portal',other_names='',cancer_type="Malignant peripheral nerve sheath tumor",species='Human',model_type='tumor')|>
    tidyr::unite(col='other_id',c('common_name','model_type'),sep=' ',remove=FALSE)

##then create samples for the PDX
sampTable<-manifest|>
    dplyr::select(common_name='Sample',MicroTissueDrugFolder)|>
    dplyr::mutate(other_id_source='NF Data Portal',other_names='',cancer_type="Malignant peripheral nerve sheath tumor",species='Human',model_type='patient derived xenograft')|>
    tidyr::unite(col='other_id',c('common_name','model_type'),sep=' ',remove=FALSE)


##third, generate a sample for the MTs if they were generated
pdxmt<-subset(sampTable,!is.na(MicroTissueDrugFolder))
pdxmt$model_type=rep('organoid',nrow(pdxmt))
print(pdxmt)

main<-rbind(sampTable,pdxmt)|>
    dplyr::select(-MicroTissueDrugFolder)|>
    rbind(tumorTable)

#main <- fread("mpnst/NF_MPNST_samples.csv")
#previous_aml <- fread(args[1])#"beatAML/beataml_samples.csv")

# If there is no previous samples file - start at 1, else, continue where the previous one left off.
if (identical(orig_samples, "")) {
    max_id <- 1  
} else {
    max_id <- max(orig_samples$improve_sample_id, na.rm = TRUE)
}


main$improve_sample_id <- seq(from = max_id + 1, length.out = nrow(main))

#synapse_main <- fread("mpnst/synapse_NF-MPNST_samples.csv")
# Step 1: Create a dictionary from 'main'
#id_dict <- setNames(main$improve_sample_id, main$other_id)

# Step 2: Update 'ID' in 'synapse_main'
#synapse_main$ID <- id_dict[synapse_main$Sample]

# Handling NA values if any mismatch occurs (Optional based on your data integrity)
# If there are NAs generated, you might need to check for unmatched keys
# synapse_main$ID[is.na(synapse_main$ID)] <- -1  # Assign a placeholder like -1 for unmatched rows

# Step 3: Save the updated 'synapse_main'
#fwrite(synapse_main, "mpnst/synapse_NF-MPNST_samples.csv")
#fwrite(main, "mpnst/NF_MPNST_samples.csv") # updated sample file
fwrite(main,'/tmp/mpnst_samples.csv')

