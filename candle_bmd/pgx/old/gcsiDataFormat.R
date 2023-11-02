

source('loadPGXdata.R')


dset<-PharmacoGx::downloadPSet('gCSI_2019')

getDoseRespData(dset,'gCSI')
##now we want to get the rna expression

