##this uses the Hg38 database to seed the gene name mappings that we might want to use

if(!require('org.Hs.eg.db')){
  BiocManager::install('org.Hs.eg.db')
  library(org.Hs.eg.db)
}

library(dplyr)
##get entrez ids to symbol
entrez<-as.data.frame(org.Hs.egALIAS2EG)

##get entriz ids to ensembl
ens<-as.data.frame(org.Hs.egENSEMBL2EG)

##get transcript ids as well
enst<-as.data.frame(org.Hs.egENSEMBLTRANS)


joined.df<-entrez%>%full_join(ens)%>%
  dplyr::rename(entrez_id='gene_id',gene_symbol='alias_symbol',other_id='ensembl_id')%>%
  mutate(other_id_source='ensembl_gene')

tdf<-entrez|>
    full_join(enst)|>
    dplyr::rename(entrez_id='gene_id',gene_symbol='alias_symbol',other_id='trans_id')|>
    dplyr::mutate(other_id_source='ensembl_transcript')

joined.df<-rbind(joined.df,tdf)|>
    distinct()

#save to file and version
write.table(joined.df,'/tmp/genes.csv',sep=',',row.names=F,quote=T)

##store this file somewhere!


