##this uses the Hg38 database to seed the gene name mappings that we might want to use

#if(!require('org.Hs.eg.db')){
#  BiocManager::install('org.Hs.eg.db')
library(org.Hs.eg.db)
library(biomaRt)
#}

library(dplyr)
##get entrez ids to symbol
entrez<-as.data.frame(org.Hs.egALIAS2EG)

##get entriz ids to ensembl
ens<-as.data.frame(org.Hs.egENSEMBL2EG)

##get transcript ids as well
enst<-as.data.frame(org.Hs.egENSEMBLTRANS)

                                        #now we can filter by protein coding using biomart
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

tab <- getBM(attributes=c('ensembl_gene_id'),filters='biotype', values=c('protein_coding'),mart=ensembl)


joined.df<-entrez%>%full_join(ens)%>%
  dplyr::rename(entrez_id='gene_id',gene_symbol='alias_symbol',other_id='ensembl_id')%>%
    mutate(other_id_source='ensembl_gene')|>
    mutate(is_protein=other_id%in%tab$ensembl_gene_id)|>
    subset(is_protein)|>
    dplyr::select(-is_protein)

tdf<-entrez|>
    full_join(enst)|>
    dplyr::rename(entrez_id='gene_id',gene_symbol='alias_symbol',other_id='trans_id')|>
    subset(entrez_id%in%joined.df$entrez_id)|>
    subset(gene_symbol%in%joined.df$gene_symbol)|>
    dplyr::mutate(other_id_source='ensembl_transcript')

joined.df<-rbind(joined.df,tdf)|>
    distinct()

#save to file and version
write.table(joined.df,'/tmp/genes.csv',sep=',',row.names=F,quote=T)

##store this file somewhere!


