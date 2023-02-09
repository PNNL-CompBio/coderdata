
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if(!require('PharmacoGx')){
  BiocManager::install("PharmacoGx")
  library('PharmacoGx')
}

library(dplyr)
library(tidyr)
source('loadPGXdata.R')
dset<-PharmacoGx::downloadPSet('PRISM_2020')
mapping <- sensitivityInfo(dset)%>%
  tibble::rownames_to_column('exp_id')
##get the raw data



##now get the drug ids
drug.map<-buildDrugTable(unique(mapping$treatmentid))%>%
  dplyr::select(common_name,candle_drug_id)%>%distinct()

##then get the sample id
full.map<-mapping%>%
  dplyr::rename(other_id='depmap_id')%>%
  left_join(subset(candle_samples,id_source=='DepMap'))%>%
  dplyr::select(exp_id,candle_sample_id,common_name='treatmentid')%>%
  distinct()%>%
  left_join(drug.map)



##now we want to select the drug candle ids and sample candle ids only

alldat <- sensitivityRaw(dset)

doseDat<-alldat[,,1]%>%
  as.data.frame()%>%
  tibble::rownames_to_column('exp_id')%>%
  tidyr::pivot_longer(cols=starts_with('dose'),names_to='doseNum',values_to='Dose')
#%>%
#  tidyr::separate(experiment,int=c('drug1','cellline','drug2','drug3'),sep='::')

respDat<-alldat[,,2]%>%
  as.data.frame()%>%
  tibble::rownames_to_column('exp_id')%>%
  tidyr::pivot_longer(cols=starts_with('dose'),names_to='doseNum',values_to='Response')
#%>%
 # 
 
doseRep<-doseDat%>%
  dplyr::full_join(respDat,by=c('doseNum','exp_id'))%>%
  subset(!is.na(Dose))%>%
  left_join(full.map)%>%
  dplyr::select(DRUG=candle_drug_id,CELL=candle_sample_id,DOSE=Dose,
                RESPONSE=Response)%>%
  mutate(SOURCE='PharmacoGx', STUDY='PRISM',GROWTH=100-RESPONSE)

print(head(doseRep))
write.table(doseRep,file='prismDoseResponse.tsv',sep='\t',row.names=F,quote=F)

