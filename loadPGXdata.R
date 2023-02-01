##get all the pgx info here.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require('PharmacoGx')){
  BiocManager::install("PharmacoGx",force=TRUE)
  library('PharmacoGx')
}

require(webchem)

library(dplyr)
library(tidyr)

all.dsets<-PharmacoGx::availablePSets()

##load existing gene/sample/drug files

##initialize experimental files
doseRep<-data.frame(DRUG=c(),CELL=c(),DOSE=c(),RESPONSE=c(),GROWTH=c(),SOURCE=c(),STUDY=c())
gex<-data.frame()
copy_number<-data.frame()
muts<-data.frame()


buildDrugTable<-function(druglist){
  print(paste("Finding ids for",length(druglist),'drugs'))
  
  pubchem_id<-webchem::get_cid(druglist)%>%
    dplyr::rename(common_name='query',pubchem_id='cid')
  
  print(paste('found',nrow(pubchem_id),'drugs'))
  
  props<-webchem::pc_prop(cids$cid)%>%
    dplyr::select(pubchem_id='CID',formula='MolecularFormula',
                  weight='MolecularWeight',canSMILES='CanonicalSMILES',
                  isoSMILES='IsomericSMILES',InChIKey='InChIKey')
  
  ##now we join them together
  joined.df<-pubchem_id%>%left_join(props)
  
  ##read in additional drugs
  if(file.exists('drugs.csv'))
    existing.drugs<-
}
