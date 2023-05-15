##get all the pgx info here.

require(webchem)
library(dplyr)
library(tidyr)

options(timeout=10000)

##load existing gene/sample/drug files

##load in gene/sample files
if(!exists('improve_genes'))
  improve_genes<<-read.csv('../data/genes.csv')
if(!exists('improve_samples'))
   improve_samples<<-read.csv('../data/samples.csv')
if(!exists('improve_drugs')){
  ##load in initial drug file
  improve_drugs<<-read.table('../data/drugs.tsv.gz',sep='\t',header=T,quote='',comment.char='')
}




#' buildDrugTable - This is a generic drug search function that
#' returns the improve_drug_id for the specific drug of interest
#' by either identifying it in the database, or querying it in
#' pubchem and then appending it
#' @return drug mapping table
buildDrugTable<-function(druglist){
  print(paste("Finding ids for",length(druglist),'drugs'))

  druglist<-unique(druglist)
  if(file.exists('../data/drugs.tsv.gz')){
      improve_drugs<<-read.table('../data/drugs.tsv.gz',sep='\t',header=T,comment.char = '',quote='')

      new_drugs<-setdiff(tolower(druglist),tolower(improve_drugs$chem_name))
      print(paste('of those drugs',length(new_drugs),'are not in database'))
      druglist<-new_drugs
  }

  if(length(druglist)>0){
    ##lets divude queries into 500 because it's easier
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
     # print(pubchem_id)


      if(nrow(pubchem_id)==0){ ## if we found no new drugs in pubchem, we need to create some
        print("no new drugs found in pubchem, creating new ones")
        #return(improve_drugs)
        props<-data.frame(formula=c(),
                          weight=c(),
                          canSMILES=c(),
                          isoSMILES=c(),
                          InChIKey=c(),
                          chem_name=c(),
                          improve_drug_id=c())
      }else{
        #now get chemical properties
        qres<-webchem::pc_prop(pubchem_id$pubchem_id,
                             properties=c('MolecularFormula','MolecularWeight',
                                              'CanonicalSMILES','IsomericSMILES','InChIKey'))

        #also get chemical synonyms
        qsyn<-webchem::pc_synonyms(pubchem_id$pubchem_id,from='cid',match='all')
        #keep all synonyms for future lookup
        syntab<-do.call('rbind',lapply(names(qsyn),function(x)
          return(data.frame(chem_name=qsyn[[x]],pubchem_id=x))))

        ##join toegether properties and synonyms
        props<-qres%>%
          subset(!is.na(CID))%>%
          mutate(pubchem_id=as.character(CID))%>%
            dplyr::select(pubchem_id,formula='MolecularFormula',
                      weight='MolecularWeight',canSMILES='CanonicalSMILES',
                      isoSMILES='IsomericSMILES',InChIKey='InChIKey')%>%
          left_join(syntab)%>%
          mutate(improve_drug_id=paste0('PC_',pubchem_id))%>%
          dplyr::select(-pubchem_id)
      }
      ##now append the missing
      missing<-setdiff(tolower(newlist),tolower(props$chem_name))

      print(paste("still missing",length(missing),'drug names, creating ids'))
      print(head(props))
      joined.df<-rbind(props,
                         data.frame(formula=rep(NA,length(missing)),
                                    weight=rep(NA,length(missing)),
                                    canSMILES=rep(NA,length(missing)),
                                    isoSMILES=rep(NA,length(missing)),
                                    InChIKey=rep(NA,length(missing)),
                                    chem_name=missing,
                                    improve_drug_id=rep(NA,length(missing))))

      if(!exists('improve_drugs')){
          improve_drugs<<-joined.df
      }
      else{
          improve_drugs<<-rbind(improve_drugs,joined.df)
      }


      ##let's add identifiers to the unmatched
      matched<-subset(improve_drugs,!is.na(improve_drug_id))
      unmatched<-subset(improve_drugs,is.na(improve_drug_id))

      if(nrow(unmatched)>0){
                     ##get max improve id
            orig<-matched%>%
              tidyr::separate(improve_drug_id,into=c('source', 'value'),sep='_')%>%
              subset(source=='IMP')

            maxval<-max(as.numeric(orig$value),na.rm=T)
            if(!is.finite(maxval))
		maxval<-0
            print(maxval)
            unmatched$improve_drug_id<-paste0('IMP_',seq(maxval+1,maxval+nrow(unmatched)))

            improve_drugs<<-rbind(matched,unmatched)
         }
      else{
           improve_drugs<<-matched
      }
      improve_drugs<<-unique(improve_drugs)
          ##write new table on every iteration that adds values in case it fails
      write.table(improve_drugs,file=gzfile('../data/drugs.tsv.gz'),sep='\t',quote=F,col.names=T,row.names=F)
    } ##end each rep of drugs

  }
  return(improve_drugs)

}



