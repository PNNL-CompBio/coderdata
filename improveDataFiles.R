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

options(timeout=10000)

all.dsets<-PharmacoGx::availablePSets()

##load existing gene/sample/drug files

##initialize experimental files

##load in gene/sample files
if(!exists('improve_genes'))
  improve_genes<<-read.csv('https://github.com/sgosline/candleDataProcessing/raw/main/data/genes.csv')
if(!exists('improve_samples'))
  improve_samples<<-read.csv('https://github.com/sgosline/candleDataProcessing/raw/main/data/samples.csv')
if(!exists('improve_drugs')){
  ##load in initial drug file
}




#' buidlDrugTable - appends drug table by querying
#' pubChem for any available matches, and then adds them to
#' the general list of all drugs
#' @return drug mapping table
buildDrugTable<-function(druglist){
  print(paste("Finding ids for",length(druglist),'drugs'))
  
  #drug.map<-data.frame()
  
  if(file.exists('data/drugs.csv')){
      improve_drugs<<-read.csv('data/drugs.tsv',sep='\t')
    
      new_drugs<-setdiff(tolower(druglist),tolower(improve_drugs$chem_name))
      print(paste('of those drugs',length(new_drugs),'are not in database'))
    #  drug.map<-subset(improve_drugs,common_name%in%improve_drugs)%>%
    #    mutate(drug_id=common_name)
      druglist<-new_drugs
  }
  
        
  if(length(druglist)>0){
    ##lets divude queries into 1000 beacuse i was getting impatient.
    qsize=500
    reps<-round(length(druglist)/qsize)
    for(i in 1:max(reps,1)){
      
      ##get a subset of the list
      newlist=druglist[(1+qsize*(i-1)):min(qsize*i,length(druglist))]
    
      #query pubchem for CIDs
      pubchem_id<-webchem::get_cid(newlist)%>%
        dplyr::rename(common_name='query',pubchem_id='cid')%>%
        subset(!is.na(pubchem_id))
    
      print(paste('found',nrow(pubchem_id),'new drugs'))
      print(pubchem_id)
      
      if(nrow(pubchem_id)==0){
        print("no new drugs found, returning existing")
        return(improve_drugs)
      }
      #now get chemical properties
      qres<-webchem::pc_prop(pubchem_id$pubchem_id,
                             properties=c('MolecularFormula','MolecularWeight',
                                          'CanonicalSMILES','IsomericSMILES','InChIKey'))
      #also get chemical synonyms
      qsyn<-webchem::pc_synonyms(pubchem_id$pubchem_id,from='cid',match='all')

      syntab<-do.call('rbind',lapply(names(qsyn),function(x)
        return(data.frame(chem_name=qsyn[[x]],pubchem_id=x))))
      
      ##join toegether properties and synonyms
      props<-qres%>%
        subset(!is.na(CID))%>%
        mutate(pubchem_id=as.character(CID))%>%
          dplyr::select(pubchem_id,formula='MolecularFormula',
                    weight='MolecularWeight',canSMILES='CanonicalSMILES',
                    isoSMILES='IsomericSMILES',InChIKey='InChIKey')%>%
        mutate(improve_chem_id=paste0('PC_',pubchem_id))%>%
        left_join(syntab)%>%
        dplyr::select(-pubchem_id)
      
      ##now we join them together
     # joined.df<-pubchem_id%>%
    #    left_join(props)
      
      ##now append the missing
      missing<-setdiff(tolower(newlist),tolower(props$chem_name))
    
      print(paste("still missing",length(missing),'drug names, creating ids'))
      
      joined.df<-rbind(props,
                       data.frame(chem_name=missing,
                                  improve_chem_id=rep(NA,length(missing)),
                                  formula=rep(NA,length(missing)),
                                  weight=rep(NA,length(missing)),
                                  canSMILES=rep(NA,length(missing)),
                                  isoSMILES=rep(NA,length(missing)),
                                  InChIKey=rep(NA,length(missing))))
    
      if(!exists('improve_drugs')){
        #joined.df$improve_drug_id<-seq(1,nrow(joined.df))
        improve_drugs<<-joined.df
      }
      else{
        #joined.df$improve_drug_id<-seq(max(improve_drugs$improve_drug_id)+1,
        #                              max(improve_drugs$improve_drug_id)+nrow(joined.df))
        improve_drugs<<-rbind(improve_drugs,joined.df)
      }
      
      ##now we add int he IMP_ ids
      idfix<-improve_drugs%>%
        separate(improve_chem_id,into=c('source','value'),sep='_')
        
      matched<-subset(idfix,!is.na(source))
      unmatched<-subset(idfix,is.na(source))
      
      ##get max impid
      imp_id<-subset(matched,source=='IMP')
      maxval<-max(as.numeric(imp_id$value),na.rm=T)
      if(is.na(maxval))
        maxval=0
      
      unmatched$source<-rep("IMP",nrow(unmatched))
      unmatched$value<-seq(maxval+1,maxval+nrow(unmatched))
      
      improve_drugs<<-rbind(matched,unmatched)%>%
        tidyr::unite(col='improve_chem_id','source','value',sep='_')
      ##write new table on every iteration in case it fails
      write.table(improve_drugs,file='data/drugs.tsv',sep='\t',quote=F,col.names=T,row.names=F)
    }
    
  }
  return(improve_drugs)
  
}



                                               