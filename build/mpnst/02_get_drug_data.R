# Load required libraries
library(data.table)
# library(biomaRt)# biomart issues still exist
library(dplyr)
library(stringr)
library(reticulate)

use_python("/opt/venv/bin/python3", required = TRUE)
source_python("pubchem_retrieval.py")

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if a token was provided
if (length(args) == 0) {
  stop("No token or sample file provided. Usage: Rscript my_script.R <PAT> [olddrugfile] [newdrugfile]", call. = FALSE)
}

# Set your personal access token
PAT <- args[1]
olddrugfiles <- args[2]
newdrugfile <- args[3]
# Log in to Synapse
library(synapser)
synLogin(authToken = PAT)


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


##copy old drug to new drug
olddrugs<-do.call(rbind,lapply(unique(unlist(strsplit(olddrugfiles,split=','))),function(x) read.table(x,header=T,sep='\t',quote='',comment.char='')))
olddrugs<-unique(olddrugs)

print(paste('Read in ',nrow(olddrugs),'old drugs'))
                                        #file.copy(olddrugfile,newdrugfile)
write.table(olddrugs,file=newdrugfile,sep='\t',row.names=F,quote=FALSE,col.names=T)
output_file_path <- newdrugfile
ignore_file_path <- '/tmp/mpnst_ignore_chems.txt'

update_dataframe_and_write_tsv(unique_names=alldrugs,output_filename=output_file_path,ignore_chems=ignore_file_path)


tab<-read.table(newdrugfile,sep='\t',header=T,quote="",comment.char="")

newdrugs<-tab|>
    subset(chem_name%in%tolower(alldrugs))

tab<-tab|>
    subset(improve_drug_id%in%newdrugs$improve_drug_id)

write.table(tab,file=newdrugfile,sep='\t',row.names=FALSE,quote=FALSE)


##now call the python drug script


