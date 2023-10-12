

source('loadPGXdata.R')

download.file('https://zenodo.org/record/5827917/files/GBM_scr3.rds?download=1',destfile='tmp.rds')
dset<-readRDS('tmp.rds')%>%
  updateObject()


dset<-PharmacoGx::downloadPSet('CTRPv2_2015')

getDoseRespData(dset,'CTRPv2')


##now we want to get the rna expression

