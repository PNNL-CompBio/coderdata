
source("loadPGXdata.R")


dset<-PharmacoGx::downloadPSet('BeatAML_2018')
mapping <- sensitivityInfo(dset)%>%
  tibble::rownames_to_column('exp_id')

##get the raw data

alldat <- sensitivityRaw(dset)


doseDat<-alldat[,,1]%>%
  as.data.frame()%>%
  tibble::rownames_to_column('exp_id')%>%
  tidyr::pivot_longer(cols=starts_with('dose'),names_to='doseNum',
                      values_to='Dose')


respDat<-alldat[,,2]%>%
  as.data.frame()%>%
  tibble::rownames_to_column('exp_id')%>%
  tidyr::pivot_longer(cols=starts_with('dose'),names_to='doseNum',values_to='Response')

 
doseRep<-doseDat%>%
  dplyr::full_join(respDat,by=c('doseNum','exp_id'))%>%
  left_join(mapping)%>%
  dplyr::select(DRUG=treatmentid,CELL=sampleid,DOSE=Dose,RESPONSE=Response)%>%
  mutate(GROWTH=100-RESPONSE,SOURCE='PharmacoGX',STUDY='BEATAML')

  
#  tidyr::separate(experiment,int=c('DRUG','CELL'),sep='_')%>%
#  dplyr::select(-doseNum)

print(head(doseRep))
write.table(doseRep,file='beatAMLDoseResponse.tsv',sep='\t',row.names=F,quote=F)


