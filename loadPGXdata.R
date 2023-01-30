##get all the pgx info here.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require('PharmacoGx')){
  BiocManager::install("PharmacoGx",force=TRUE)
  library('PharmacoGx')
}

library(dplyr)
library(tidyr)

all.dsets<-PharmacoGx::availablePSets()

##load existing gene/sample/drug files

##initialize experimental files
doseRep<-data.frame(DRUG=c(),CELL=c(),DOSE=c(),RESPONSE=c(),GROWTH=c(),SOURCE=c(),STUDY=c())
gex<-data.frame()
copy_number<-data.frame()
muts<-data.frame()
