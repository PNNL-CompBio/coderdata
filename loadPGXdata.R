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
doseRep<-data.frame(DRUG=c(),CELL=c(),DOSE=c(),RESPONSE=c(),
                    GROWTH=c(),SOURCE=c(),STUDY=c())
gex<-data.frame()
copy_number<-data.frame()
muts<-data.frame()

##load in gene/sample files
if(!exists('candle_genes'))
  candle_genes<<-read.csv('https://github.com/sgosline/candleDataProcessing/raw/main/data/genes.csv')
if(!exists('candle_samples'))
  candle_samples<<-read.csv('https://github.com/sgosline/candleDataProcessing/raw/main/data/samples.csv')
if(!exists('candle_drugs')){
  ##load in initial drug file
}
##we have to build the drug table with every new addition.
##need to return a drug.map that has all the drug information
##and also the 'drug_id' of the druglist that was used as input
buildDrugTable<-function(druglist){
  print(paste("Finding ids for",length(druglist),'drugs'))
  
  #drug.map<-data.frame()
  
  if(!exists('candle_drugs')){
    if(file.exists('drugs.csv')){
      candle_drugs<<-read.csv('drugs.csv',sep=',')
    
      new_drugs<-setdiff(druglist,candle_drugs$common_name)
      print(paste('of those drugs',length(new_drugs),'are not in database'))
    #  drug.map<-subset(candle_drugs,common_name%in%candle_drugs)%>%
    #    mutate(drug_id=common_name)
      druglist<-new_drugs
  }}
  
  
        
  if(length(druglist)>0){
    
    pubchem_id<-webchem::get_cid(druglist)%>%
      dplyr::rename(common_name='query',pubchem_id='cid')
    
    print(paste('found',nrow(pubchem_id),'new drugs'))
    
    props<-webchem::pc_prop(pubchem_id$pubchem_id)%>%
    mutate(pubchem_id=as.character(CID))%>%
        dplyr::select(pubchem_id,formula='MolecularFormula',
                  weight='MolecularWeight',canSMILES='CanonicalSMILES',
                    isoSMILES='IsomericSMILES',InChIKey='InChIKey')
    
    ##now we join them together
    joined.df<-pubchem_id%>%left_join(props)
    missing<-setdiff(druglist,joined.df$common_name)
    
    print(paste("still missing",length(missing),'drug names'))
    joined.df<-rbind(joined.df,
                     data.frame(common_name=missing,
                                pubchem_id=rep(NA,length(missing)),
                                formula=rep(NA,length(missing)),
                                weight=rep(NA,length(missing)),
                                canSMILES=rep(NA,length(missing)),
                                isoSMILES=rep(NA,length(missing)),
                                InChIKey=rep(NA,length(missing))))
    
    if(!exists('candle_drugs')){
      joined.df$candle_drug_id<-seq(1,nrow(joined.df))
      candle_drugs<<-joined.df
    }
    else{
      joined.df$candle_drug_id<-seq(max(candle_drugs$candle_drug_id)+1,
                                    max(candle_drugs$candle_drug_id)+nrow(joined.df))
      candle_drugs<<-rbind(candle_drugs,joined.df)
    }
  }
  return(candle_drugs)
}
                                               