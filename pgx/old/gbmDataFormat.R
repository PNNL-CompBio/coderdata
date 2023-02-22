
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require('PharmacoGx')){
  BiocManager::install("PharmacoGx")
  library('PharmacoGx')
}

library(dplyr)
library(tidyr)
all.dsets<-PharmacoGx::availablePSets()


dset1<-PharmacoGx::downloadPSet('GBM_2021_scr3')
mapping <- sensitivityInfo(dset1)%>%
  tibble::rownames_to_column('exp_id')


alldat <- sensitivityRaw(dset1)

doseDat<-alldat[,,1]%>%
  as.data.frame()%>%
  tibble::rownames_to_column('exp_id')%>%
  tidyr::pivot_longer(cols=starts_with('dose'),names_to='doseNum',values_to='Dose')
#%>%
#  tidyr::separate(experiment,int=c('drug1','cellline','drug2','drug3'),sep='::')

respDat<-alldat[,,2]%>%
  as.data.frame()%>%
  tibble::rownames_to_column('exp_id')%>%
  tidyr::pivot_longer(cols=starts_with('dose'),names_to='doseNum',values_to='Response')
#%>%
# 

doseRep1<-doseDat%>%
  dplyr::full_join(respDat,by=c('doseNum','exp_id'))%>%
  subset(!is.na(Dose))%>%
  left_join(mapping)%>%
  dplyr::select(DRUG=treatmentid,CELL=sampleid,DOSE=Dose,
                RESPONSE=Response)%>%
  mutate(SOURCE='PharmacoGx', STUDY='GBM_scr3',GROWTH=100-RESPONSE)

###########################################################

dset2<-PharmacoGx::downloadPSet('GBM_2021_scr2')
mapping2 <- sensitivityInfo(dset2)%>%
  tibble::rownames_to_column('exp_id')


alldat <- sensitivityRaw(dset2)

doseDat<-alldat[,,1]%>%
  as.data.frame()%>%
  tibble::rownames_to_column('exp_id')%>%
  tidyr::pivot_longer(cols=starts_with('dose'),names_to='doseNum',values_to='Dose')
#%>%
#  tidyr::separate(experiment,int=c('drug1','cellline','drug2','drug3'),sep='::')

respDat<-alldat[,,2]%>%
  as.data.frame()%>%
  tibble::rownames_to_column('exp_id')%>%
  tidyr::pivot_longer(cols=starts_with('dose'),names_to='doseNum',values_to='Response')
#%>%
# 

doseRep2<-doseDat%>%
  dplyr::full_join(respDat,by=c('doseNum','exp_id'))%>%
  subset(!is.na(Dose))%>%
  left_join(mapping2)%>%
  dplyr::select(DRUG=treatmentid,CELL=sampleid,DOSE=Dose,
                RESPONSE=Response)%>%
  mutate(SOURCE='PharmacoGx', STUDY='GBM_scr2',GROWTH=100-RESPONSE)

cdoseRep<-rbind(doseRep,doseRep2)
print(head(cdoseRep))
write.table(cdoseRep,file='gbmDoseResponse.tsv',sep='\t',row.names=F,quote=F)

#####expression parsring

##dset1
expression<-molecularProfiles(dset1,'rna')%>%
  as.data.frame()%>%
  tibble::rownames_to_column('geneid')%>%
  tidyr::pivot_longer(cols=starts_with('GSM'),names_to='sample',values_to='tpm')

genemapping<-featureInfo(dset1,'rna')%>%
  as.data.frame()%>%
  tibble::rownames_to_column('geneid')%>%
  right_join(expression)

##this is the final table with all possible columns. now we need to extract those we need for hte schema
sampGeneMapping<-phenoInfo(dset1,'rna')%>%
  as.data.frame()%>%
  tibble::rownames_to_column('sample')%>%
  right_join(genemapping)


##dset2
molecularProfiles(dset12,'rna')
