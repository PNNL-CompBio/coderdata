##this file creates a sample database table from Priyanka's depmap file
##and cellosaurus. it ONLY handles cell lines. ADditional samples must be added
##subsequent to running this code. 

##author sara gosline
## email sara.gosline@pnnl.gov

require(remotes)
if(!require(TidyComb)){
  BiocManager::install("synergyfinder")
  remotes::install_github("DrugComb/TidyComb")
}

library(curl)
library(dplyr)


##the only thing that Priyanka has here is TRP identifiers, so collecting those
tab<-read.table('DepMap_Argonne_Mapping.csv',sep=',',header=T)%>%
  dplyr::select(Argonne_ID,DepMap_ID)%>%
  distinct()%>%
  tidyr::separate(Argonne_ID,into=c('id_source','other_id'),sep='\\.')%>%
  subset(id_source=='CTRP')%>%
  dplyr::rename(CTRP='other_id',ModelID='DepMap_ID')%>%
  dplyr::select(-id_source)


##here are all the models in depmap, downloded on 1/31/2023
depmap_models<-read.table('Model.csv',sep=',',header=T)

##query for cellosaurus automagically
url='https://ftp.expasy.org/databases/cellosaurus/cellosaurus.xml'
curl_download(url,'cell.xml',quiet=TRUE)#curl(url, "r", h)
cello<-TidyComb::ParseCell('cell.xml')
cdf<-XML::xmlToDataFrame(cello)
cols<-c('accession-list','name-list','species-list','disease-list')

cdf.red<-cdf%>%
  dplyr::select(all_of(cols))%>%
  distinct()%>%
  subset(!is.na(`disease-list`)) #try to remove the non cancer stuff

##get humans cell lines only
cdf.red<-cdf.red[grep('Homo sapiens',cdf.red$`species-list`),]%>%
  dplyr::rename(RRID='accession-list')%>%
  dplyr::select(-c(`species-list`))

######
#now we join the depmap table, the cellosaurus table, and the CTRP identifiers

joined.df<-depmap_models%>%
  full_join(cdf.red)%>%
  left_join(tab)%>%
  dplyr::select(-c(PatientID,SourceType,GrowthPattern,
                     PrimaryOrMetastasis,MolecularSubtype,
                   CatalogNumber,PublicComments,SampleCollectionSite,
                   Sex,Age,SourceDetail,OncotreeLineage,OncotreeCode,
                   OncotreePrimaryDisease,DepmapModelType,StrippedCellLineName))

##now lengethn the table to have the appopriate columns
full.df<-joined.df%>%
  dplyr::rename(DepMap='ModelID',Sanger='SangerModelID',
                CCLE='CCLEName',Cellosaurus='RRID',
                common_name='CellLineName',cancer_type='OncotreeSubtype',
                other_names='name-list')%>%
  mutate(COSMIC=as.character(COSMICID),WTSI=as.character(WTSIMasterCellID))%>%
  dplyr::select(-c(COSMICID,WTSIMasterCellID))

##this is hte full table, we add identifiers HEREEEE
full.df<-full.df%>%
  mutate(candle_sample_id=seq(1:nrow(full.df)))

long.df<-full.df%>%
  tidyr::pivot_longer(cols=c(DepMap,Sanger,CCLE,COSMIC,WTSI,CTRP,Cellosaurus),names_to='id_source',
                      values_to='other_id')%>%
  mutate(model_type='cell line')


write.table(long.df,'samples.csv',sep=',',row.names=F,col.names=T)

