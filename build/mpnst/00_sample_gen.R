# This script generate a new sample table based on pervious beatAML improved sample ID
# It will take the maximum value of beatAML improved sample ID and continue from ID count from there
# Load required libraries
library(data.table)
library(synapser)
library(dplyr)

##adding a command line argument
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=2){
    stop("Need a sample file and synapse token as argument. Rscript 00_sample_gen.R [samplefile] [synapse token]")

}

orig_samples<-fread(args[1])

synapser::synLogin(authToken=args[2])
manifest<-synapser::synTableQuery("select * from syn53503360")$asDataFrame()


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
    dplyr::mutate(other_id_source='NF Data Portal',cancer_type="Malignant peripheral nerve sheath tumor",species='Human',model_type='Tumor')|>
    tidyr::unite(col='other_id',c('common_name','model_type'),sep=' ',remove=FALSE)

##then create samples for the PDX
sampTable<-manifest|>
    dplyr::select(c(common_name='Sample',MicroTissueDrugFolder))|>
    dplyr::mutate(other_id_source='NF Data Portal',cancer_type="Malignant peripheral nerve sheath tumor",species='Human',model_type='Patient derived xenograft')|>
    tidyr::unite(col='other_id',c('common_name','model_type'),sep=' ',remove=FALSE)


##third, generate a sample for the MTs if they were generated
pdxmt<-subset(sampTable,!is.na(MicroTissueDrugFolder))
pdxmt$model_type=rep('organoid',nrow(pdxmt))
#print(pdxmt)

main<-rbind(sampTable,pdxmt)|>
    dplyr::select(-MicroTissueDrugFolder)|>
    rbind(tumorTable)

#main <- fread("mpnst/NF_MPNST_samples.csv")
#previous_aml <- fread(args[1])#"beatAML/beataml_samples.csv")
max_id <- max(orig_samples$improve_sample_id)
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
fwrite(main,'/tmp/MPNST_samples.csv')

