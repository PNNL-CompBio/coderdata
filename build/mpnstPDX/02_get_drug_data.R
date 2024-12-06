# Load required libraries
library(data.table)
# library(biomaRt)# biomart issues still exist
library(dplyr)
library(stringr)
library(reticulate)
library(synapser)
library(tidyr)


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





##define functions

#print(pdx)
##now loop through manifest to get all the files
pdx_fold <- data.table(pdx)[,strsplit(as.character(PDX_Drug_Data),","), by = .(common_name)]|>
    subset(!is.na(V1))|>
    subset(V1!='NA')|>
    dplyr::rename(id='V1')

#print(pdx_fold)
###this is not all of themju
pdx_meta<-do.call(rbind,lapply(pdx_fold$id, function(x) synapser::synGetAnnotations(x)|>
                                          as.data.frame()|>
                                          dplyr::select('experimentalCondition')|>
                                          dplyr::mutate(id=x)))|>
    left_join(pdx_fold)|>
    tidyr::separate_rows(experimentalCondition,sep=';')|>
    mutate(chem_name=tolower(experimentalCondition))

#pdx_drug <- data.table(pdx_meta)[,strsplit(as.character(experimentalCondition),';'),by= .(common_name,id)]|>
#    mutate(drug=tolower(experimentalCondition))
#drugs<-sapply(pdx_meta$experimentalCondition,function(x) tolower(unlist(strsplit(x,split=';'))))|>
#    unlist()|>
#    unique()

drugs<-setdiff(pdx_meta$chem_name,'control')


print(paste(drugs,collapse=','))


##copy old drug to new drug
olddrugs<-do.call(rbind,lapply(unique(unlist(strsplit(olddrugfiles,split=','))),function(x) read.table(x,header=T,sep='\t',quote='',comment.char='')))
olddrugs<-unique(olddrugs)

print(paste('Read in ',nrow(olddrugs),'old drug files'))

fdrugs<-subset(olddrugs,chem_name%in%drugs)
if(nrow(fdrugs)>0){
    dids<-fdrugs$improve_drug_id
}else{
    dids<-c()
}
newdrugs<-subset(olddrugs,improve_drug_id%in%dids)

print(paste('Found',length(dids),'improved drug ids that exist, saving those'))


                                        #file.copy(olddrugfile,newdrugfile)
write.table(newdrugs,file=newdrugfile,sep='\t',row.names=F,quote=FALSE,col.names=T)
output_file_path <- newdrugfile
ignore_file_path <- '/tmp/mpnstpdx_ignore_chems.txt'


##now load reticulate down here



use_python("/opt/venv/bin/python3", required = TRUE)
source_python("pubchem_retrieval.py")

update_dataframe_and_write_tsv(unique_names=drugs,output_filename=output_file_path,ignore_chems=ignore_file_path)


tab<-read.table(newdrugfile,sep='\t',header=T,quote="",comment.char="")

newdrugs<-tab|>
    subset(chem_name%in%tolower(alldrugs))

tab<-tab|>
    subset(improve_drug_id%in%newdrugs$improve_drug_id)

write.table(tab,file=newdrugfile,sep='\t',row.names=FALSE,quote=FALSE)


##now call the python drug script


