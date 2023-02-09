

files<-commandArgs(trailingOnly=TRUE)
sampFile<-files[1]
geneFile<-files[2]
drugFile<-files[3]

source('loadPGXdata.R')
samps<-read.csv(sampFile)
genes<-read.csv(geneFile)

if(file.exists(drugFile))
  drugs<-read.csv(drugFile)

dset<-PharmacoGx::downloadPSet('CTRPv2_2015')
mapping <- sensitivityInfo(dset)%>%
  tibble::rownames_to_column('exp_id')



##first lets try to map samples
samp_ids<-unique(mapping$sampleid)

matches<-lapply(samp_ids,function(x) grep(x,samps$other_names,fixed=TRUE))

#all_ids<-unique(samps$other_id)
#  dplyr::rename(other_id='sampleid')%>%##this maps it to the file

    ##now map the sample information to sample file. here, grep is hte most reliable
    
    
##get the raw data

alldat <- sensitivityRaw(dset)

doseDat<-alldat[,,1]%>%
  as.data.frame()%>%
  tibble::rownames_to_column('exp_id')%>%
  tidyr::pivot_longer(cols=starts_with('dose'),names_to='doseNum',values_to='Dose')


respDat<-alldat[,,2]%>%
  as.data.frame()%>%
  tibble::rownames_to_column('exp_id')%>%
  tidyr::pivot_longer(cols=starts_with('dose'),names_to='doseNum',values_to='Response')


doseRep<-doseDat%>%
  dplyr::full_join(respDat,by=c('doseNum','exp_id'))%>%
  left_join(mapping)%>%
  dplyr::select(DRUG=treatmentid,CELL=sampleid,DOSE=Dose,RESPONSE=Response)%>%
  mutate(GROWTH=100-RESPONSE)%>%
  mutate(SOURCE='pharmacoGX')%>%
  mutate(STUDY='FIMM')

drug.tab<-buildDrugTable(unique(doseRep$DRUG))

write.table(drug.tab,file='data/drugs.csv',sep=',')

print(head(doseRep))

##now map  drugs, samples, genes to ids in database files

#write.table(doseRep,file='fimmDoseResponse',sep='\t',row.names=F,quote=F)

##there are no other data files available, so we can write empty files



