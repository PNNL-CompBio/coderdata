##this file creates a sample database table from Priyanka's depmap file

require(remotes)
if(!require(TidyComb)){
  BiocManager::install("synergyfinder")
  remotes::install_github("DrugComb/TidyComb")
}

library(curl)
library(dplyr)



tab<-read.table('DepMap_Argonne_Mapping.csv',sep=',',header=T)

##here are all the models in depmap
depmap_models<-read.table('Model.csv',sep=',',header=T)

##query for cellosaurus
url='https://ftp.expasy.org/databases/cellosaurus/cellosaurus.xml'
curl_download(url,'cell.xml',quiet=TRUE)#curl(url, "r", h)
cello<-TidyComb::ParseCell('cell.xml')
cdf<-XML::xmlToDataFrame(cello)
cols<-c('accession-list','name-list','species-list','disease-list')

##add column for cell line type
cdf.red<-cdf%>%
  dplyr::select(all_of(cols))%>%
  distinct()%>%
  subset(!is.na(`disease-list`))

##get humans cell lines only
cdf.red<-cdf.red[grep('Homo sapiens',cdf.red$`species-list`),]%>%
  dplyr::rename(RRID='accession-list')%>%
  dplyr::select(-c(`species-list`,`derived-from`))

##add column for disease 
joined.df<-depmap_models%>%
  full_join(cdf.red)

##now lengethn the table to have the appopriate columns
long.df<-joined.df%>%
  dplyr::rename(DepMap='ModelID',Sanger='SangerModelID',
                CCLE='CCLEName')%>%
  mutate(COSMIC=as.character(COSMICID),WTSI=as.character(WTSIMasterCellID))%>%
  dplyr::select(-c(COSMICID,WTSIMasterCellID))%>%
  tidyr::pivot_longer(cols=c(DepMap,Sanger,CCLE,COSMIC,WTSI),names_to='id_source',
                      values_to='id_values')

write.table(long.df,'samples.tsv',sep='\t',row.names=F,col.names=T)

