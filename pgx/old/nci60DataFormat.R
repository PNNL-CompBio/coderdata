

source('loadPGXdata.R')


dset<-PharmacoGx::downloadPSet('NCI60_2021')

getDoseRespData(dset,'NCI60')

##now we want to get the rna expression
geneExp<-assays(molecularProfiles(dset)$rnaseq.iso)$exprs%>%
  as.data.frame()%>%
  tibble::rownames_to_column('gene')%>%
  tidyr::pivot_longer(cols=seq(2,1+ncol(assays(molecularProfiles(dset)$rnaseq.iso)$exprs)),
                      names_to='samples',
                      values_to='counts')

##now we map samples
geneExp<-geneExp%>%
  dplyr::rename(other_names='samples')%>%
  left_join(candle_samples)%>%
  subset(!is.na(candle_sample_id))%>%
  select(gene,counts,candle_sample_id)%>%
  distinct()%>%
  dplyr::rename(gene_symbol='gene')%>%
  left_join(candle_genes)%>%
  dplyr::select(entrez_id,candle_sample_id,counts)%>%
  subset(!is.na(entrez_id))%>%
  mutate(source='PharmacoGX',study='NCI60')

#now we map gene names
write.table(doseRep,file='nci60GeneExp.csv',sep=',',row.names=F,quote=F)

