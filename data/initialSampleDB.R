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
tab<-read.table('data/DepMap_Argonne_Mapping.csv',sep=',',header=T)%>%
  dplyr::select(Argonne_ID,DepMap_ID)%>%
  distinct()%>%
  tidyr::separate(Argonne_ID,into=c('id_source','other_id'),sep='\\.')%>%
  subset(id_source=='CTRP')%>%
  dplyr::rename(CTRP='other_id',ModelID='DepMap_ID')%>%
  dplyr::select(-id_source)


##here are all the models in depmap, downloded on 1/31/2023
depmap_models<-read.table('data/Model.csv',sep=',',header=T)

##query for cellosaurus automagically
url='https://ftp.expasy.org/databases/cellosaurus/cellosaurus.xml'
curl_download(url,'cell.xml',quiet=TRUE)#curl(url, "r", h)
cello<-XML::xmlParse('cell.xml')
cdf<-XML::xmlToList(cello)



##these are the data we need
##ok this command seems to have gotten file in appropriate state
cell.lines<-lapply(cdf$`cell-line-list`, function(x) unlist(x))

##now we need toe xtract columns
options(show.error.messages=TRUE)
full.res<-do.call(rbind,lapply(cell.lines,function(x){
  ##create a data frame for each cell lines
  x<-unlist(x)
  #should only be one acession
  acc<-x[grep('accession.text',names(x),fixed=T)]%>%unlist()
  
  cn<-x[grep('name.text',names(x),fixed=T)]%>%unlist()
  #these will fail if no key found
  spec<-x[grep("species-list.cv-term.text",names(x),fixed=T)]%>%unlist()
  #dis<-x[grep("disease-list.cv-term.text",names(x),fixed=T)]%>%unlist()
  data.frame(accession=cn,
             RRID=rep(acc,length(cn)),
             species=rep(spec,length(cn)))
   #          disease=rep(dis,length(cn)))
}))%>%
  subset(species=='Homo sapiens')



######
#now we join the depmap table, the cellosaurus table, and the CTRP identifiers

joined.df<-depmap_models%>%
  full_join(full.res)%>%
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
                other_names='accession')%>%
  mutate(COSMIC=as.character(COSMICID),WTSI=as.character(WTSIMasterCellID))%>%
  dplyr::select(-c(COSMICID,WTSIMasterCellID))%>%
  distinct()

#we add idnetifiers for everything that has a cellosaurus id and those that dont
has_id<-subset(full.df,Cellosaurus!="")
no_id<-subset(full.df,Cellosaurus=="")

samp_ids<-data.frame(Cellosaurus=unique(has_id$Cellosaurus))
samp_ids$candle_sample_id<-seq(1,nrow(samp_ids))

#now we need to add in the missing ones
extra_ids<-data.frame(DepMap=no_id$DepMap)
extra_ids$candle_sample_id<-seq(max(samp_ids$candle_sample_id)+1,
                                max(samp_ids$candle_sample_id)+nrow(extra_ids))

full.df<-rbind(left_join(has_id,samp_ids),
               left_join(no_id,extra_ids))

long.df<-full.df%>%
  tidyr::pivot_longer(cols=c(DepMap,Sanger,CCLE,COSMIC,WTSI,CTRP,Cellosaurus),names_to='id_source',
                      values_to='other_id')%>%
  mutate(model_type='cell line')%>%
  subset(!is.na(other_id))


write.table(long.df,'data/samples.csv',sep=',',row.names=F,col.names=T)

