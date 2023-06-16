
library(readr)
library(tidyr)

###first reqad in all gene information so we can map appropriately
allgenes = read_csv("https://raw.githubusercontent.com/PNNL-CompBio/candleDataProcessing/main/data/genes.csv")
genes = allgenes|>
  dplyr::select(gene_symbol,entrez_id)|>
  dplyr::distinct()


##here are the improve sample id indices
samples = read_csv('https://raw.githubusercontent.com/PNNL-CompBio/candleDataProcessing/main/data/samples.csv',
                   quote='"')|>
  dplyr::select(other_id,improve_sample_id)

Sys.setenv(VROOM_CONNECTION_SIZE=10000000)
###PATH TO FILES
###THESE ARE FILES FROM THE DEPMAP repository from PRIAY
basename='https://ftp.mcs.anl.gov/pub/candle/public/improve/Data/Omics/Curated_CCLE_Multiomics_files/'
filenames=list(expression='CCLE_AID_expression_full.csv',
                            copy_number='CCLE_AID_gene_cn.csv',
                            mutations='Mutation_AID_binary.csv',
                            methylation='CCLE_AID_RRBS_TSS_1kb_20180614.csv',
                            proteins='CCLE_AID_RPPA_20180123.csv',
                            mirnas='CCLE_AID_miRNA_20180525.csv')#wont' work because no gene identifiers

###run through each file and rewrite
newres<-lapply(names(filenames),function(value){

  fi=paste0(basename,filenames[[value]])
  fname=paste0(value,'.csv')

  ##now every data type is parsed slightly differently, so we need to change our formatting
  ##and mapping to get it into a unified 3 column schema
  if(value=='copy_number'){
    exp_file <- readr::read_csv(fi,skip=1)[-1,]

    res = exp_file|>
      tidyr::pivot_longer(cols=c(2:ncol(exp_file)),
                          names_to='gene_symbol',values_to='alteration',
                          values_transform=list(alteration=as.numeric))|>
      dplyr::distinct()

    res<-res|> ##deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
      dplyr::mutate(copy_call=ifelse(alteration<0.5210507,'deep del',
                                     ifelse(alteration<0.7311832,'het loss',
                                            ifelse(alteration<1.214125,'diploid',
                                                   ifelse(alteration<1.422233,'gain','amp')))))

    colnames(res)[1]<-'other_id'
    vars=c('copy_number','alteration')


  }else if(value=='methylation'){ ###IF DATA REPRESENT RRBS###
    exp_file <- readr::read_csv(fi)[-c(1:2),]
    #the gene names are not unique and in some weird format, i willt ry to keep both in the metadata
    res = exp_file|>
      tidyr::pivot_longer(cols=c(2:ncol(exp_file)),
                              names_to='gene_region',values_to='methylation',
                          values_transform=list(methylation=as.numeric))|>
      dplyr::distinct()

    res<-res|>
      tidyr::separate(gene_region,into=c('gene_symbol','num','start','end'),sep='_')
    colnames(res)[1]<-'other_id'
    vars=c('methylation','start','end')


    }else if(value=='mutations'){ ####IF DATA REPRESENTS MUTATIONS#####
      exp_file <- readr::read_csv(fi)[-1,]

      ##current file has 1024 patient samples, so let's read them in ~100 at a time
      pat<-grep('CVCL',colnames(exp_file))
      npats<-length(pat)
      reps=seq(1,ceiling(npats/100))
      res<-lapply(reps,function(x){
        print(x)
        cols<-seq(min(pat)+(x-1)*100,min(ncol(exp_file),min(pat)-1+(x)*100)) ##this can keep changing!
        dres<-exp_file[,c('Entrez_id','Genome_Change',colnames(exp_file)[cols])]
        ret <-tidyr::pivot_longer(dres,cols=c(3:ncol(dres)),names_to='other_id',values_to='mutations')
        #print(head(ret))
        ret|>
          subset(mutations!=0)
      })

      res<-do.call(rbind,res)
      full<-res|>  ###since we're already in ENTREZ we skip the mapping below
        dplyr::left_join(samples)|>
        dplyr::rename(entrez_id=Entrez_id,alteration=Genome_Change)|>
        dplyr::select(entrez_id,improve_sample_id,alteration)|>
        dplyr::distinct()|>
        dplyr::mutate(source='DepMap',study='CCLE')

        write_csv(full,file=fname)
        return(fi)
    }
    else if(value=='expression'){ #if gene expression
      exp_file <- readr::read_csv(fi,skip=2)

      res = tidyr::pivot_longer(data=exp_file,cols=c(2:ncol(exp_file)),
                                names_to='gene_symbol',values_to='expression',
                                values_transform=list(expression=as.numeric))|>
        dplyr::left_join(genes)|>
        dplyr::distinct()
      colnames(res)[1]<-'other_id'
      vars=c('expression')

    }
    else if(value=='mirnas'){ #if mirna expression
      exp_file <- readr::read_csv(fi)

      res = tidyr::pivot_longer(data=exp_file,cols=c(2:ncol(exp_file)),
                                names_to='gene_symbol',values_to='mirnas',values_transform=list(mirnas=as.numeric))|>
        dplyr::left_join(genes)|>
        dplyr::distinct()|>
        mutate(gene_symbol=tolower(gene_symbol))

      colnames(res)[1]<-'other_id'
      vars=c('mirnas')

    }
    else if(value=='proteins'){ #if protein expression
      exp_file <- readr::read_csv(fi)
      ###WARNING: this will not match most genes!!! no gene names are provided!!!
      res = tidyr::pivot_longer(data=exp_file,cols=c(2:ncol(exp_file)),
                                names_to='gene_symbol',values_to=proteins,
                                values_transform=list(proteins=as.numeric))|>
        dplyr::left_join(genes)|>
        dplyr::distinct()|>
        mutate(gene_symbol=toupper(gene_symbol))

      colnames(res)[1]<-'other_id'
      vars=c('proteins')
    }

  ##do the last join with samples
  full<-res|>
    dplyr::left_join(samples)|>
    dplyr::select(c('entrez_id','improve_sample_id',vars))|>
    dplyr::distinct()|>
    dplyr::mutate(source='DepMap',study='CCLE')

  write_csv(full,file=fname)
  return(fi)

})

