###Here is a script that runs through all the data files one by one.

#this is a helper file that loads the data
source("mapDrugsToPubchem.R")

if(!require('PharmacoGx')){
  BiocManager::install("PharmacoGx",force=TRUE)
  library('PharmacoGx')
}

all.dsets<-PharmacoGx::availablePSets()
##first define a generic dose response function

improve_samples<<-readr::read_csv('cell_line_samples.csv')
#' getDoseRespData
#' Generic function to get dose and response data from PGX object
#' out of dataset object, and store with dataset name
getDoseRespData<-function(dset,studyName){

  mapping <- sensitivityInfo(dset)##get the dataset dose response data

  ##map the rownames
  if(!'exp_id'%in%names(mapping))
    mapping<-mapping%>%
      tibble::rownames_to_column('exp_id')

  ##fix up the NSC ids
  if('NSC'%in%names(mapping))
    mapping<-mapping|>
      dplyr::select(-treatmentid)|>
      dplyr::mutate(treatmentid=paste0('NSC-',NSC))

  ##move drug to treatment id
  if("drugid"%in%names(mapping))
    mapping<-dplyr::rename(mapping,treatmentid='drugid')

  ##move cellid to sampleid
  if('cellid'%in%names(mapping))
    mapping<-dplyr::rename(mapping,sampleid='cellid')

  ##query to build the drug ids
  drug.map<-buildDrugTable(unique(mapping$treatmentid,'drugs.tsv.gz'))%>%
    dplyr::select(common_drug_name='chem_name',improve_drug_id)%>%
    distinct()

  #reduce drug ids to only one pubchem id, in case there are more!
  red.drug.map<-drug.map|>
    subset(tolower(common_drug_name)%in%tolower(unique(mapping$treatmentid)))|>
    tidyr::separate(improve_drug_id,into=c('pc','num'),sep='_')

  ##now for each PubChem/Improve duplicate, just get minimum id number
  minvals<-red.drug.map|>
    group_by(common_drug_name)|>
    summarize(minVal=min(num))

  ##now reduce drug map to those minvals
  new.drug.map<-red.drug.map|>
    subset(num%in%minvals$minVal)|>
    tidyr::unite('pc','num',col='improve_drug_id')

  ##first get the sample id
  samp.map<-mapping%>%
    dplyr::select(sampleid,exp_id,treatmentid)%>%distinct()

  comm.map<-samp.map%>%
    dplyr::rename(other_names='sampleid')%>%
    left_join(improve_samples)%>%
    subset(!is.na(improve_sample_id))

  print(paste('By common name, found',length(unique(comm.map$improve_sample_id)),
              'matches out of',length(unique(mapping$sampleid)),'for study',studyName))

  ldmap<-new.drug.map%>%
    mutate(chem_name=tolower(common_drug_name))
  ##then join the sample id
  full.map<-comm.map%>%
    dplyr::select(exp_id,improve_sample_id,treatmentid)%>%
    mutate(chem_name=tolower(treatmentid))%>%
    distinct()%>%
    left_join(ldmap)

  #all_ids<-unique(samps$other_id)
  #  dplyr::rename(other_id='sampleid')%>%##this maps it to the file

  ##now map the sample information to sample file. here, grep is hte most reliable


  ##get the raw data
  ##now
  alldat <- sensitivityRaw(dset)

  doseDat<-alldat[,,1]%>%
    as.data.frame()%>%
    tibble::rownames_to_column('exp_id')%>%
    tidyr::pivot_longer(cols=starts_with('dose'),names_to='doseNum',values_to='Dose')


  respDat<-alldat[,,2]%>%
    as.data.frame()%>%
    tibble::rownames_to_column('exp_id')%>%
    tidyr::pivot_longer(cols=starts_with('dose'),names_to='doseNum',values_to='Response')


  doseRep<-doseDat%>%
    dplyr::full_join(respDat,by=c('doseNum','exp_id'))%>%
    left_join(full.map)%>%
    dplyr::select(DRUG=improve_drug_id,CELL=improve_sample_id,DOSE=Dose,GROWTH=Response)%>%
    #dplyr::mutate(DOSE=-log10(Dose/1000))###curve fitting code requires -log10(M), these are mM
    #rename(GROWTH=RESPONSE)%>%
    mutate(SOURCE='pharmacoGX')%>%
    mutate(STUDY=studyName)

  print(head(doseRep))

  ##now map  drugs, samples, genes to ids in database files

  write.table(doseRep,file=paste0(tolower(studyName),'DoseResponse'),sep='\t',row.names=F,quote=F)
  doseRep

}

###################################################################
##now we iterate through the drug files one by one
####




#' getCellLineData - gets cell line dose response data
getCellLineDoseData<-function(cell.lines=c('CTRPv2','FIMM','gCSI','PRISM','GDSC','NCI60','CCLE')){
  ###first get cell lines
  all.dose.rep<-do.call(rbind,lapply(cell.lines,function(cel){

  # print(cel)
   files<-subset(all.dsets,`Dataset Name`==cel)%>%
      dplyr::select(`PSet Name`)%>%
      unlist()


    res<-do.call(rbind,lapply(files,function(f){
      print(f)
      if(f=='GDSC_2020(v2-8.2)')
        cel='GDSCv2'
      if(f=='GDSC_2020(v1-8.2)')
        cel='GDSCv1'
      tmpfile<-paste0(cel,'doseResponse')

      if(file.exists(tmpfile)){
        dres<-read.table(tmpfile,sep='\t',header=T)
      }else{
        dset<<-downloadPSet(f,saveDir='.',timeout=10000)

        url=subset(all.dsets,`PSet Name`==f)$Download
      #print(url)

        dres<-getDoseRespData(dset,cel)
       # print(dres)
      }
      return(dres)

    }))
      #print(res)
      return(res)
  }))

  #print(head(all.dose.rep))
  all.dose.rep%>%
    subset(!is.na(DOSE))
  #write.table(all.dose.rep,file='allDoseRepPreCalc.tsv',sep='\t')
}

if(FALSE){

  cl1<-c('CTRPv2','FIMM','GDSC')
  dl1<-getCellLineDoseData(cl1)

  cl2<-c('gCSI','PRISM','CCLE')
  dl2<-getCellLineDoseData(cl2)

  cl2<-c('NCI60') ###this is the biggest dataset by far, and has lots of drugs that require lookup
  dl2<-getCellLineDoseData(cl2)

}
