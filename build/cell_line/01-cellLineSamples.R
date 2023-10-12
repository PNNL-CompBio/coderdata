##this file creates a sample database table from Priyanka's depmap file
##and cellosaurus. it ONLY handles cell lines. ADditional samples must be added
##subsequent to running this code.

##author sara gosline
## email sara.gosline@pnnl.gov

library(curl)
library(dplyr)
library(readr)

##the only thing that Priyanka has here is TRP identifiers, so collecting those
# tab<-read.table('DepMap_Argonne_Mapping.csv',sep=',',header=T)%>%
#   dplyr::select(Argonne_ID,DepMap_ID)%>%
#   distinct()%>%
#   tidyr::separate(Argonne_ID,into=c('id_source','other_id'),sep='\\.')%>%
#   subset(id_source=='CTRP')%>%
#   dplyr::rename(CTRP='other_id',ModelID='DepMap_ID')%>%
#   dplyr::select(-id_source)
##we miss these
###[1] "Panc-05-04"              "L3-3"                    "Panc-02-03"              "Panc-10-05"              "G-292-clone-A141B1"      "Hs-729"                  "Panc-03-27"              "Panc-04-03"              "Hs-688-A-T"
###[10] "Panc-08-13"              "PE-CA-PJ34-clone-C12"    "Hep-3B2-1-7"             "PE-CA-PJ41-clone-D2"     "PE-CA-PJ49"              "Ishikawa-Heraklio-02-ER"

##here are all the models in depmap 23Q2, downloded on 10/11/2023
depmap_models<-readr::read_csv('https://figshare.com/ndownloader/files/40448834')#|>
#  dplyr::rename(DepMap_ID='ModelID')


##query for cellosaurus automagically to get loadest version
url='https://ftp.expasy.org/databases/cellosaurus/cellosaurus.xml'
if(!file.exists('cell.xml'))
  curl_download(url,'cell.xml',quiet=TRUE)#curl(url, "r", h)
cello<-XML::xmlParse('cell.xml')
cdf<-XML::xmlToList(cello)



### next we parse through cellosaurus to get as many samples as we deem relevant
##ok, this command seems to have gotten file in appropriate state
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
  subset(species=='Homo sapiens (Human)')


print(paste('Got',nrow(full.res),'human cellosaurus samples'))
######
#now we join the depmap table, the cellosaurus table, and the CTRP identifiers

joined.df<-depmap_models%>%
  full_join(full.res)%>%
  #left_join(tab)%>%
  ##this for 22q2 data
  #dplyr::select(-c(patient_id,sample_collection_site,primary_or_metastasis,primary_disease,
  #                 sex,age,source,depmap_public_comments,lineage,lineage_sub_subtype,Subtype,Cellosaurus_NCIt_id,lineage_subtype,
  #                 lineage_molecular_subtype,default_growth_pattern,stripped_cell_line_name,alias,source,sample_collection_site,
  #                 model_manipulation,model_manipulation_details,parent_depmap_id,Cellosaurus_issues))
  ##use this for 22q4 data
  dplyr::select(-c(PatientID,SourceType,GrowthPattern,
                     PrimaryOrMetastasis,LegacyMolecularSubtype,LegacySubSubtype,
                   EngineeredModel,TreatmentStatus,PlateCoating,AgeCategory,
                   CatalogNumber,PublicComments,SampleCollectionSite,OnboardedMedia,
                   Sex,Age,SourceDetail,OncotreeLineage,OncotreeCode,
                   OncotreePrimaryDisease,DepmapModelType,StrippedCellLineName))

##now lengethn the table to have the appopriate columns
full.df<-joined.df%>%
  dplyr::rename(DepMap='ModelID',Sanger='SangerModelID',##old file formatSanger='SangerModelID'
                CCLE='CCLEName',#old file CCLE='CCLEName'
                Cellosaurus='RRID',
                 common_name='CellLineName',
                cancer_type='OncotreeSubtype',
                other_names='accession')%>%
  mutate(COSMIC=as.character(COSMICID),WTSI=as.character(WTSIMasterCellID))%>%
  dplyr::select(-c(COSMICID,WTSIMasterCellID))%>%
  distinct()

#we add idnetifiers for everything that has a cellosaurus id and those that dont
has_id<-subset(full.df,Cellosaurus!="")
no_id<-subset(full.df,Cellosaurus=="")

samp_ids<-data.frame(Cellosaurus=unique(has_id$Cellosaurus))
samp_ids$improve_sample_id<-seq(1,nrow(samp_ids))

#now we need to add in the missing ones
if(nrow(no_id)>0){
  extra_ids<-data.frame(DepMap=no_id$DepMap)
  extra_ids$improve_sample_id<-seq(max(samp_ids$improve_sample_id)+1,
                                  max(samp_ids$improve_sample_id)+nrow(extra_ids))

}
full.df<-rbind(left_join(has_id,samp_ids),
               left_join(no_id,extra_ids))

long.df<-full.df%>%
  tidyr::pivot_longer(cols=c(DepMap,Sanger,CCLE,COSMIC,WTSI,Cellosaurus),names_to='id_source',
                      values_to='other_id')%>%
  mutate(model_type='cell line')%>%
  subset(!is.na(other_id))%>%
  subset(other_id!="")


write.table(long.df,'samples.csv',sep=',',row.names=F,col.names=T)

