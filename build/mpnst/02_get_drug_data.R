# Load required libraries
library(data.table)
# library(biomaRt)# biomart issues still exist
library(dplyr)
library(stringr)
library(synapser)


# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)


# Check the number of arguments provided
if (length(args) < 1) {
  stop("At least one argument is required. Usage: Rscript 02_get_drug_data.R <newdrugfile> [olddrugfile]", call. = FALSE)
}


# Assign arguments
newdrugfile <- args[1]  # Path to the new drug file
olddrugfiles <- ifelse(length(args) >= 2 && args[2] != "", args[2], NA)

# Read SYNAPSE_AUTH_TOKEN from the environment
synapse_token <- Sys.getenv("SYNAPSE_AUTH_TOKEN")
if (synapse_token == "") {
  stop("Error: SYNAPSE_AUTH_TOKEN environment variable is not set.")
}

synLogin(authToken = synapse_token)

##now get the manifest from synapse
manifest<-synapser::synTableQuery("select * from syn53503360")$asDataFrame()|>
                                                             as.data.frame()|>
                                                             dplyr::rename(common_name='Sample')


##PDX contain list of files
pdx<-manifest|>
    dplyr::select(common_name,PDX_Drug_Data)|>
    distinct()|>
    subset(!is.na(PDX_Drug_Data))


##MTS contain lists of directories
mts<-manifest|>
    dplyr::select(common_name,MicroTissueDrugFolder)|>
    subset(!is.na(MicroTissueDrugFolder))



##define functions

##first function to get children from parentId
getDrugsByParent<-function(parid){
    qtab<-synTableQuery(paste('select id,name,experimentalCondition,parentId from syn21993642 where parentId=\'',parid,'\''))$asDataFrame()|>
        subset(!is.na(experimentalCondition))|>dplyr::select(id,name,experimentalCondition)
    ##now we need to parse the metadatda table get the info

    return(unique(qtab$experimentalCondition))

}

##now loop through manifest to get all the files
mts_fold <- data.table(mts)[,strsplit(as.character(MicroTissueDrugFolder),","), by = .(common_name)]

alldrugs<-unique(unlist(lapply(mts_fold$V1,function(x){
    samp<-subset(mts_fold,V1==x)
    res<-getDrugsByParent(x)
    return(res)
})))


alldrugs[which(alldrugs=='PD901')]<-'PD-0325901'

print(paste(alldrugs,collapse=','))



## new code:


# Handle old drugs
if (!is.na(olddrugfiles)) {
  # Read and combine old drug files
  olddrug_list <- lapply(unique(unlist(strsplit(olddrugfiles, split = ','))), function(x) {
    if (file.exists(x)) {
      return(fread(x, header = TRUE, sep = '\t', quote = ''))
    } else {
      warning(paste("Old drug file does not exist:", x))
      return(NULL)
    }
  })
  
  # Remove NULL entries and ensure uniqueness
  olddrug_list <- Filter(Negate(is.null), olddrug_list)
  
  if (length(olddrug_list) > 0) {
    olddrugs <- unique(rbindlist(olddrug_list, use.names = TRUE, fill = TRUE))
    print(paste('Read in', nrow(olddrugs), 'old drugs'))
  } else {
    olddrugs <- data.frame(
      improve_drug_id = integer(),
      chem_name = character(),
      pubchem_id = character(),
      canSMILES = character(),
      # isoSMILES = character(),
      InChIKey = character(),
      formula = character(),
      weight = numeric(),
      stringsAsFactors = FALSE
    )
    print("Old drug files not valid. Created empty olddrugs dataframe.")
  }
} else {
  # Create an empty dataframe with specified columns
  olddrugs <- data.frame(
    improve_drug_id = integer(),
    chem_name = character(),
    pubchem_id = character(),
    canSMILES = character(),
    # isoSMILES = character(),
    InChIKey = character(),
    formula = character(),
    weight = numeric(),
    stringsAsFactors = FALSE
  )
  print("No old drug file provided. Created empty olddrugs dataframe.")
}

# Write the initial drug file (old drugs)
write.table(olddrugs, file = newdrugfile, sep = '\t', row.names = FALSE, quote = FALSE,col.names=T)


# Define the ignore file path
ignore_file_path <- '/tmp/mpnst_ignore_chems.txt'


# ##copy old drug to new drug
# olddrugs<-do.call(rbind,lapply(unique(unlist(strsplit(olddrugfiles,split=','))),function(x) read.table(x,header=T,sep='\t',quote='',comment.char='')))
# olddrugs<-unique(olddrugs)

# print(paste('Read in ',nrow(olddrugs),'old drugs'))
#                                         #file.copy(olddrugfile,newdrugfile)
# write.table(olddrugs,file=newdrugfile,sep='\t',row.names=F,quote=FALSE,col.names=T)


##now load reticulate down here

library(reticulate)

use_python("/opt/venv/bin/python3", required = TRUE)
source_python("pubchem_retrieval.py")

update_dataframe_and_write_tsv(unique_names=alldrugs,output_filename=newdrugfile,ignore_chems=ignore_file_path)


tab<-read.table(newdrugfile,sep='\t',header=T,quote="",fill=TRUE)

newdrugs<-tab|>
    subset(chem_name%in%tolower(alldrugs))

tab<-tab|>
    subset(improve_drug_id%in%newdrugs$improve_drug_id)

write.table(tab,file=newdrugfile,sep='\t',row.names=FALSE,quote=FALSE)

print(paste("Final drug table written to", newdrugfile))


##now call the python drug script


