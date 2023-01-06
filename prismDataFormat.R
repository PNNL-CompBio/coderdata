
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if(!require('PharmacoGx')){
  BiocManager::install("PharmacoGx")
  library('PharmacoGx')
}

library(dplyr)
library(tidyr)
all.dsets<-PharmacoGx::availablePSets()

dset<-PharmacoGx::downloadPSet('PRISM_2020')

##get the raw data

alldat <- sensitivityRaw(dset)

doseDat<-alldat[,,1]%>%
  as.data.frame()%>%
  tibble::rownames_to_column('experiment')%>%
  tidyr::pivot_longer(cols=starts_with('dose'),names_to='doseNum',values_to='Dose')
#%>%
#  tidyr::separate(experiment,int=c('drug1','cellline','drug2','drug3'),sep='::')

respDat<-alldat[,,2]%>%
  as.data.frame()%>%
  tibble::rownames_to_column('experiment')%>%
  tidyr::pivot_longer(cols=starts_with('dose'),names_to='doseNum',values_to='Response')
#%>%
 # 
 
doseRep<-doseDat%>%
  dplyr::outer_join(respDat,by=c('doseNum','experiment'))%>%
  tidyr::separate(experiment,int=c('drug1','cellline','drug2','drug3'),sep='::')


