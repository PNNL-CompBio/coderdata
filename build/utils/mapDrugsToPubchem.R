##get all the pgx info here.

require(webchem)
library(dplyr)
library(tidyr)
library(readr)
options(timeout=10000)
##load existing gene/sample/drug files

#' buildDrugTable - This is a generic drug search function that
#' returns the improve_drug_id for the specific drug of interest
#' by either identifying it in the database, or querying it in
#' pubchem and then appending it
#' @return drug mapping table
buildDrugTable<-function(druglist,path_to_file='drugs.tsv.gz'){
  print(paste("Finding ids for",length(druglist),'drugs'))

  druglist<-unique(druglist)

  max_id=0
  if(file.exists(path_to_file)){
      improve_drugs<<-read.table(path_to_file,sep='\t',header=T,comment.char = '',quote='')

      res = improve_drugs|>
          tidyr::separate(improve_drug_id,sep='_',into=c('imp','id'))
      max_id=max(as.numeric(res$id))

      new_drugs<-setdiff(tolower(druglist),tolower(improve_drugs$chem_name))
      print(paste('of those drugs',length(new_drugs),'are not in database'))
      druglist<-new_drugs
  }##reduce druglist to those that are no in list

  if(length(druglist)>0){
    ##lets divude queries into 500 because it's easier
    qsize=500
    reps<-round(length(druglist)/qsize)
    for(i in 1:max(reps,1)){
      ##get a subset of the list
        newlist=druglist[(1+qsize*(i-1)):min(qsize*i,length(druglist))]

      #query pubchem for CIDs
        pubchem_id<-webchem::get_cid(newlist)%>%
            rbind(.,webchem::get_cid(newlist,domain='compound'))|>
            dplyr::rename(common_name='query',pubchem_id='cid')%>%
            subset(!is.na(pubchem_id))|>
            distinct()

        print(paste('found',length(unique(pubchem_id$pubchem_id)),'new drugs'))
     # print(pubchem_id)
        missed_ids <-setdiff(tolower(newlist),tolower(pubchem_id$common_name))


        ## print(paste('found',length(unique(pubchem_id$pubchem_id)),'new drugs, missed',length(missed_ids),'trying to find pubchem substance'))
        ## ##maybe we need to look for suvstances too?

        ## substance_id<-webchem::get_cid(missed_ids)%>%
        ##     rbind(.,webchem::get_cid(newlist,domain='substance'))|>
        ##     dplyr::rename(common_name='query',pubchem_id='cid')%>%
        ##     subset(!is.na(pubchem_id))|>
        ##     distinct()

        ## missed_ids <-setdiff(tolower(missed_ids),tolower(substance_id$common_name))

        ## if(nrow(substance_id)>0)
        ##     pubchem_id<-rbind(pubchem_id,substance_id)

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
            webchem::get_cid(nsc_val$query,domain='compound'))|> ##get cdis from substance
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

#      if(nrow(pubchem_id)==0){ ## if we found no new drugs in pubchem, we need to create some
#        print("no new drugs found in pubchem")
        #return(improve_drugs)
#        props<-data.frame(formula=c(),
#                          weight=c(),
#                          canSMILES=c(),
#                          isoSMILES=c(),
#                          InChIKey=c(),
#                          chem_name=c(),
#                          improve_drug_id=c())
#      }else{
        if(nrow(pubchem_id)>0){
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
                left_join(names_only)|>
                subset(!is.na(chem_name))
                                        #%>%
                                        #          mutate(improve_drug_id=paste0('PC_',pubchem_id))%>%
                                        #          dplyr::select(-pubchem_id)


            if(exists('improve_drugs')){
                ##are any of those props already in pubchemID?
                inpubchem<-subset(props,isoSMILES%in%improve_drugs$isoSMILES)
                print(paste(nrow(inpubchem),'drugs entries are already in the file by SMILE string'))

                print('Adding new improve ids for new drugs')
                iid_map<-improve_drugs|>
                    subset(isoSMILES%in%inpubchem$isoSMILES)|>
                    dplyr::select(isoSMILES,improve_drug_id)|>
                    distinct()
                new_i <- setdiff(props$isoSMILES,improve_drugs$isoSMILES)
                iid_map<-rbind(iid_map,data.frame(isoSMILES=new_i,improve_drug_id=paste0('SMI_',seq(max_id+1,max_id+length(new_i)))))

            }else{
                new_i <-unique(props$isoSMILES)
                iid_map<-data.frame(isoSMILES=new_i,improve_drug_id=paste0('SMI_',seq(max_id+1,max_id+length(new_i))))
            }


            ##now append the missing
            missing<-setdiff(tolower(newlist),tolower(props$chem_name))

            print(paste("still missing",length(missing),'drug names:'))
            print(head(missing))
                                        #    print(tail(props))
            joined.df <- props|>
                left_join(iid_map)|>
                distinct()

            if(!exists('improve_drugs')){
                improve_drugs<<-joined.df
            }
            else if(nrow(joined.df)>0){
                print(head(improve_drugs))
                print(head(joined.df))
                improve_drugs<<-rbind(improve_drugs,joined.df)
            }


            improve_drugs<<-unique(improve_drugs)
        }


          ##write new table on every iteration that adds values in case it fails
      write.table(improve_drugs,file=gzfile(path_to_file),sep='\t',
                  quote=F,col.names=T,row.names=F)
    } ##end each rep of drugs

  }
  return(improve_drugs)

}



