
source('loadPGXdata.R')

dset<-PharmacoGx::downloadPSet('FIMM_2016')
mapping <- sensitivityInfo(dset)%>%
  tibble::rownames_to_column('exp_id')

##get the raw data

##now get the drug ids
drug.map<-buildDrugTable(unique(mapping$treatmentid))%>%
  dplyr::select(common_drug_name='common_name',candle_drug_id)%>%distinct()

##first get the sample id
samp.map<-mapping%>%
  dplyr::select(sampleid,exp_id,treatmentid)%>%distinct()

comm.map<-samp.map%>%
  dplyr::rename(other_names='sampleid')%>%
  left_join(candle_samples)%>%
  subset(!is.na(candle_sample_id))

print(paste('By common name, found',length(unique(comm.map$candle_sample_id)),
            'matches out of',length(unique(mapping$sampleid))))


##then get the sample id
full.map<-comm.map%>%
  dplyr::select(exp_id,candle_sample_id,common_drug_name='treatmentid')%>%
  distinct()%>%
  left_join(drug.map)

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
  dplyr::select(DRUG=candle_drug_id,CELL=candle_sample_id,DOSE=Dose,RESPONSE=Response)%>%
  mutate(GROWTH=100-RESPONSE)%>%
  mutate(SOURCE='pharmacoGX')%>%
  mutate(STUDY='FIMM')


print(head(doseRep))

##now map  drugs, samples, genes to ids in database files

write.table(doseRep,file='fimmDoseResponse',sep='\t',row.names=F,quote=F)

##there are no other data files available, so we can write empty files



