
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if(!require('PharmacoGx')){
  BiocManager::install("PharmacoGx")
  library('PharmacoGx')
}

library(dplyr)
library(tidyr)
all.dsets<-PharmacoGx::availablePSets()

dset<-PharmacoGx::downloadPSet('Tavor_2020')

##get the raw data

alldat <- sensitivityRaw(dset)

mapping <- sensitivityInfo(dset)

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
  left_join(mapping)%>%
  dplyr::select(treatmentid,sampleid,Dose,Response)

#  tidyr::separate(experiment,int=c('DRUG','CELL'),sep='_')%>%
#  dplyr::select(-doseNum)

print(head(doseRep))
write.table(doseRep,file='tavorDoseResponse.tsv',sep='\t',row.names=F,quote=F)

#tavor also has expression data, so we need to extract counts

