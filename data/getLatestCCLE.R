
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
                            'Mutation_AID_binary2.csv'))
#                            'Mutation_AID_count.csv'))#,
                            #'CCLE_AID_RPPA_20180123.csv' #wont' work because no gene identifiers

newres<-lapply(filenames,function(fi){

  fname=stringr::str_replace(basename(fi),'AID','IID')

  if(length(grep('cn',fi))>0){
    exp_file <- readr::read_csv(fi,skip=1)[-1,]

  }else if(length(grep('RRBS',fi))>0){
    exp_file <- readr::read_csv(fi,skip=3)

    }else if(length(grep("Mutation",fi))>0){
      exp_file <- readr::read_csv(fi)[-1,]

      ##current file has 1024 patient samples, so let's read them in ~100 at a time
      npats<-length(grep('CVCL',colnames(exp_file)))
      reps=seq(1,ceiling(npats/100))
      res<-lapply(reps,function(x){
        print(x)
        cols<-seq(13+(x-1)*100,min(ncol(exp_file),12+(x)*100))
        dres<-exp_file[,c(3,12,cols)]
        ret <-tidyr::pivot_longer(dres,cols=c(3:ncol(dres)),names_to='other_id',values_to='mutation')
        print(head(ret))
        ret|>
          subset(mutation!=0)
      })
      res<-do.call(rbind,res)
      full<-res|>
        dplyr::left_join(samples)|>
        dplyr::rename(entrez_id=Entrez_Gene_Id,alteration=Genome_Change)|>
        dplyr::select(entrez_id,improve_sample_id,alteration)|>
        dplyr::distinct()|>
        dplyr::mutate(source='DepMap',study='CCLE')
      write_csv(full,file=fname)
      return(fi)
      }else{ #if gene expression
    exp_file <- readr::read_csv(fi,skip=2)
  }

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

  if(length(grep('cn',fi)>1)){
    full<-dplyr::rename(full,alteration=counts)
    full<-full|> ##deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
      dplyr::mutate(copy_call=ifelse(alteration<0.5210507,'deep del',
                               ifelse(alteration<0.7311832,'het loss',
                                       ifelse(alteration<1.214125,'diploid',
                                               ifelse(alteration<1.422233,'gain','amp')))))
  }

  write_csv(full,file=fname)
  return(fi)

})

