
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if(!require('PharmacoGx')){
  BiocManager::install("PharmacoGx")
  library('PharmacoGx')
}

library(dplyr)
library(tidyr)
all.dsets<-PharmacoGx::availablePSets()


dset<-PharmacoGx::downloadPSet('FIMM_2016')
mapping <- sensitivityInfo(dset)%>%
  tibble::rownames_to_column('exp_id')

##get the raw data

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
  left_join(mapping)%>%
  dplyr::select(DRUG=treatmentid,CELL=sampleid,DOSE=Dose,RESPONSE=Response)%>%
  mutate(GROWTH=100-RESPONSE)%>%
  mutate(SOURCE='pharmacoGX')%>%
  mutate(STUDY='FIMM')

print(head(doseRep))
write.table(doseRep,file='fimmDoseResponse.tsv',sep='\t',row.names=F,quote=F)

