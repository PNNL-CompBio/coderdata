
source('loadPGXdata.R')

dset<-PharmacoGx::downloadPSet('FIMM_2016')

getDoseRespData(dset,'FIMM')


