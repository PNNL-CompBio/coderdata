##this uses the Hg38 database to seed the gene name mappings that we might want to use

#if(!require('org.Hs.eg.db')){
#  BiocManager::install('org.Hs.eg.db')
#library(org.Hs.eg.db)
library(biomaRt)
#}

library(dplyr)
##get entrez ids to symbol
# entrez<-as.data.frame(org.Hs.egALIAS2EG)
#
# sym <- as.data.frame(org.Hs.egSYMBOL)
#
# ##get entriz ids to ensembl
# ens<-as.data.frame(org.Hs.egENSEMBL2EG)
#
# ##get transcript ids as well
# enst<-as.data.frame(org.Hs.egENSEMBLTRANS)

                                        #now we can filter by protein coding using biomart
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

tab <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol','entrezgene_id','ensembl_transcript_id'),
             filters = 'biotype',
             values = c('protein_coding'), mart = ensembl) ##only protein coding

new.df <- tab |>
  dplyr::rename(entrez_id = 'entrezgene_id', gene_sybmol = 'hgnc_symbol',#rename existing columns
                ensembl_transcript = 'ensembl_transcript_id', ensembl_gene = 'ensembl_gene_id') |>  #consistent with previous versions
  subset(!is.na(entrez_id)) |> # we dont want missing entrez ids
  tidyr::pivot_longer(cols = c('ensembl_gene','ensembl_transcript'), #compress these
                      names_to = 'other_id_source',
                      values_to = 'other_id')

# tab <- getBM(attributes=c('ensembl_gene_id'),filters='biotype', values=c('protein_coding'),mart=ensembl)
#
#
# joined.df<-entrez|>
#     left_join(sym)|>
#     dplyr::rename(entrez_id='gene_id',gene_symbol='symbol',other_id='alias_symbol',gene_symbol='symbol')%>%
#     mutate(other_id_source='entrez_alias')
#
# ##now get aliases from ensembl
# edf <- sym|>
#     inner_join(ens)|>
#     dplyr::rename(entrez_id='gene_id',gene_symbol='symbol',other_id='ensembl_id')%>%
#     mutate(other_id_source='ensembl_gene')
#
#
# tdf<-sym|>
#     inner_join(enst)|>
#     dplyr::rename(entrez_id='gene_id',gene_symbol='symbol',other_id='trans_id')|>
#     subset(entrez_id%in%edf$entrez_id)|>
# #    subset(gene_symbol%in%ed.df$gene_symbol)|>
#     dplyr::mutate(other_id_source='ensembl_transcript')
#
#
# prots<-subset(edf,other_id%in%tab$ensembl_gene_id)
#
# full.df<-rbind(joined.df,edf,tdf)|>
#     subset(entrez_id%in%prots$entrez_id)|>
#     distinct()

#save to file and version
write.table(new.df,'/tmp/genes.csv',sep = ',', row.names = F, quote = T)

##store this file somewhere!

