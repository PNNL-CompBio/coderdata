
library(readr)
library(tidyr)
genes = read_csv('genes.csv')|>
  dplyr::select(gene_symbol,entrez_id)|>
  dplyr::distinct()

samples = read_csv('samples.csv',quote='"')|>
  dplyr::select(other_id,improve_sample_id)

Sys.setenv(VROOM_CONNECTION_SIZE=10000000)

basename='https://ftp.mcs.anl.gov/pub/candle/public/improve/Data/Omics/Curated_CCLE_Multiomics_files/'
filenames=paste0(basename,c('CCLE_AID_expression_full.csv',
                            'CCLE_AID_RRBS_TSS_1kb_20180614.csv',
                            'CCLE_AID_gene_cn.csv',
                            'Mutation_AID_count.csv'))#,
                            #'CCLE_AID_RPPA_20180123.csv' #wont' work because no gene identifiers

newres<-lapply(filenames,function(fi){
  if(length(grep('cn',fi))>0||length(grep("Mutation",fi))>0){
    exp_file <- readr::read_csv(fi,skip=1)[-1,]
  }else if(length(grep('RRBS',fi))>0){
    exp_file <- readr::read_csv(fi,skip=3)

    }else{ #if gene expression
    exp_file <- readr::read_csv(fi,skip=2)
  }
  fname=stringr::str_replace(basename(fi),'AID','IID')

  #colnames(exp_file)<-str(colnames(exp_file))
  ## reshape and rename
  res = tidyr::pivot_longer(exp_file,cols=c(2:ncol(exp_file)),
                            names_to='gene_symbol',values_to='counts',values_transform=list(counts=as.numeric))|>
    dplyr::left_join(genes)|>
    dplyr::distinct()

  colnames(res)[1]<-'other_id'

  full<-res|>
    dplyr::left_join(samples)|>
    dplyr::select(entrez_id,improve_sample_id,counts)|>
    dplyr::distinct()|>
    dplyr::mutate(source='DepMap',study='CCLE')

  write_csv(full,file=fname)
  return(fi)

})

