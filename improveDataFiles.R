##get all the pgx info here.

require(webchem)
library(dplyr)
library(tidyr)
library(readr)
options(timeout=10000)
library(rfigshare)
##load existing gene/sample/drug files

##load in gene/sample files directly from figshare
if(!exists('improve_genes'))
  improve_genes<<-read.csv('../cell_line/genes.csv')
if(!exists('improve_samples'))
  # improve_samples<<-fs_download('40576103')#read.csv('../cell_line/samples.csv')
  improve_samples<<-read_csv('cell_lines/samples.csv')
if(!exists('improve_drugs')){
  ##load in initial drug file
  improve_drugs<<-read.table('cell_line/drugs.tsv.gz',sep='\t',header=T,quote='',comment.char='')
}
if(!exists('experiments')){
  experiments <<-readr::read_tsv('cell_line/experiments.tsv.gz')
}


#' lookupDrug
getDataForDrug<-function(drug_name){


}


#' lookupSample
getDataForSample<-function(sample_name){


}


#' getDrugData - this is a generic function that gets dose response data by drug or sample id
#' or both



#' buildDrugTable - This is a generic drug search function that
#' returns the improve_drug_id for the specific drug of interest
#' by either identifying it in the database, or querying it in
#' pubchem and then appending it
#' @return drug mapping table
buildDrugTable<-function(druglist){
  print(paste("Finding ids for",length(druglist),'drugs'))

  druglist<-unique(druglist)
  if(file.exists('../cell_line/drugs.tsv.gz')){
      improve_drugs<<-read.table('../cell_line/drugs.tsv.gz',sep='\t',header=T,comment.char = '',quote='')

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
        rbind(.,webchem::get_cid(newlist,domain='substance'))|>
        dplyr::rename(common_name='query',pubchem_id='cid')%>%
        subset(!is.na(pubchem_id))|>
        distinct()

      print(paste('found',length(unique(pubchem_id$pubchem_id)),'new drugs'))
     # print(pubchem_id)
      missed_ids <-setdiff(tolower(newlist),tolower(pubchem_id$common_name))
      print(paste('found',length(unique(pubchem_id$pubchem_id)),'new drugs, missed',length(missed_ids),'trying to find NSC'))

      nscs<-missed_ids[grep('nsc-',missed_ids)]
      #print(paste(length(nscs),'of those are NSC'))
      missed_nsc_id<-data.frame()### create an empty one for joining belo
      if(length(nscs)>0){

        ##create a new table that permutes the NSC ids a bit to deal with wobble
        nsc_val=do.call(rbind,lapply(nscs, function(x){
          data.frame(common_name=x,query=c(gsub('-','',x),gsub('-',' ',x)))
        }))##first create substute values

        ##do new query
        new_id=rbind(webchem::get_cid(nsc_val$query),  ##get CIDs from compound
            webchem::get_cid(nsc_val$query,domain='substance'))|> ##get cdis from substance
          dplyr::rename(pubchem_id='cid')|>
          subset(!is.na(pubchem_id))

        ##append the matches to the list
        if(nrow(new_id)>0){
          print(paste('found',length(unique(new_id$pubchem_id)),'missing NSC ids'))
          missed_nsc_id<-nsc_val|>
            dplyr::inner_join(new_id)|>
            dplyr::select(-query)

          pubchem_id=rbind(pubchem_id,missed_nsc_id)|>
            distinct()

        }
      }
      print(paste('now have',length(unique(pubchem_id$pubchem_id)),'matched drugs'))

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
                                              'CanonicalSMILES','IsomericSMILES','InChIKey',
                                          'IUPACName'))


        #also get chemical synonyms
        qsyn<-webchem::pc_synonyms(pubchem_id$pubchem_id,from='cid',match='all')

        #keep all synonyms for future lookup
        syntab<-do.call('rbind',lapply(names(qsyn),function(x)
          return(data.frame(common_name=qsyn[[x]],pubchem_id=x))))|>
          rbind(pubchem_id)|>
          dplyr::rename(chem_name='common_name')|>
          distinct()

        ##now add in the IUPAC names
        names_only<-qres|>
          dplyr::select(pubchem_id=CID,chem_name=IUPACName)|>
          rbind(syntab)|>
          distinct()

        ##now add in the names that were missed with our NSC munging
        if(nrow(missed_nsc_id)>0){
          names_only<-rbind(names_only,dplyr::rename(missed_nsc_id,chem_name='common_name'))
        }

        ##join toegether properties and synonyms
        props<-qres%>%
          subset(!is.na(CID))%>%
          mutate(pubchem_id=as.character(CID))%>%
            dplyr::select(pubchem_id,formula='MolecularFormula',
                      weight='MolecularWeight',canSMILES='CanonicalSMILES',
                      isoSMILES='IsomericSMILES',InChIKey='InChIKey')%>%
          left_join(names_only)%>%
          mutate(improve_drug_id=paste0('PC_',pubchem_id))%>%
          dplyr::select(-pubchem_id)
      }
      ##now append the missing
      missing<-setdiff(tolower(newlist),tolower(props$chem_name))

      print(paste("still missing",length(missing),'drug names, creating ids'))
      print(head(missing))
      print(tail(props))
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
      write.table(improve_drugs,file=gzfile('../data/drugs.tsv.gz'),sep='\t',
                  quote=F,col.names=T,row.names=F)
    } ##end each rep of drugs

  }
  return(improve_drugs)

}



