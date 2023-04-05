


source('loadPGXdata.R')

download.file("https://zenodo.org/record/3905462/files/CCLE.rds?download=1",destfile='tmp.rds')
dset<-readRDS('tmp.rds')%>%
  updateObject()

getDoseRespData(dset,'CCLE')

##there are no other data files available, so we can write empty files


edat<-molecularProfiles(dset)$Kallisto_0.46.1.rnaseq.counts
sampmap<-colData(edat)%>%
  as.data.frame()%>%
  tibble::rownames_to_column('samples')%>%
  dplyr::select('samples','sampleid')%>%
  distinct()
##now we want to get the rna expression
geneExp<-assays(edat)$exprs%>%
  as.data.frame()%>%
  tibble::rownames_to_column('ensgene')%>%
  tidyr::separate(ensgene,into=c('gene','vers'),'\\.')%>%
  tidyr::pivot_longer(cols=seq(3,2+ncol(assays(edat)$exprs)),
                      names_to='samples',
                      values_to='counts')%>%
  left_join(sampmap)

cmap<-candle_samples%>%
  dplyr::select(other_id,candle_sample_id)%>%
  distinct()

##now we map samples
geneExp2<-geneExp%>%
  dplyr::rename(other_id='sampleid')%>%
  left_join(cmap)%>%
  subset(!is.na(candle_sample_id))%>%
  select(gene,counts,candle_sample_id)%>%
  distinct()%>%
  dplyr::rename(other_id='gene')%>%
  left_join(candle_genes)%>%
  dplyr::select(entrez_id,candle_sample_id,counts)%>%
  subset(!is.na(entrez_id))%>%
  mutate(source='PharmacoGX',study='CCLE')

#now we map gene names
write.table(geneExp2,file='ccleGeneExp.csv',sep=',',row.names=F,quote=F)
