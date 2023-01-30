
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require('PharmacoGx')){
  BiocManager::install("PharmacoGx")
  library('PharmacoGx')
}

library(dplyr)
library(tidyr)
all.dsets<-PharmacoGx::availablePSets()

dset<-PharmacoGx::downloadPSet('PDTX_2019')
mapping <- sensitivityInfo(dset)%>%
  tibble::rownames_to_column('exp_id')
##get the raw data

alldat <- sensitivityRaw(dset)