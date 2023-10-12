
library(readr)
library(tidyr)
library(dplyr)
###first reqad in all gene information so we can map appropriately
allgenes = read_csv("https://raw.githubusercontent.com/PNNL-CompBio/candleDataProcessing/main/cell_line/genes.csv")
genes = allgenes|>
  dplyr::select(gene_symbol,entrez_id)|>
  dplyr::distinct()



##here are the improve sample id indices
samples = read_csv('https://raw.githubusercontent.com/PNNL-CompBio/candleDataProcessing/main/cell_line/samples.csv',
                   quote='"')|>
  dplyr::select(other_id,improve_sample_id)|>
  unique()

Sys.setenv(VROOM_CONNECTION_SIZE=100000000)
###PATH TO FILES
###THESE ARE FILES FROM THE DEPMAP repository from PRIAY
basename='https://ftp.mcs.anl.gov/pub/candle/public/improve/Data/Omics/Curated_CCLE_Multiomics_files/'
filenames=list(transcriptomics='CCLE_AID_expression_full.csv',
                            copy_number='CCLE_AID_gene_cn.csv',
                            methylation='CCLE_AID_RRBS_TSS_1kb_20180614.csv',
                           # proteins='CCLE_AID_RPPA_20180123.csv', #wont' work well because no gene identifiers
                            miRNA='CCLE_AID_miRNA_20180525.csv',
                            mutations='Mutation_AID_binary.csv')


getProteomics<-function(){
  #pull directly from gygi lab
  proteomics <- 'https://gygi.hms.harvard.edu/data/ccle/Table_S2_Protein_Quant_Normalized.xlsx'
  options(timeout=300)

  res<-download.file(proteomics,'prot.xlsx')
  pdat<-readxl::read_xlsx('prot.xlsx',sheet = 'Normalized Protein Expression')
  pdat[,7:ncol(pdat)]<-apply(pdat[,7:ncol(pdat)],2,as.numeric)
  pdat<-pdat|>
    dplyr::select(!starts_with('Ten'))
  plong<-pdat|>
    tidyr::pivot_longer(cols=7:ncol(pdat),names_to='cellLine',values_to='proteomics',values_drop_na=TRUE)

  pfilt<-plong|>
    dplyr::select(gene_symbol='Gene_Symbol',cellLine,proteomics)|>
    dplyr::distinct()|>
    tidyr::separate(cellLine,into=c('other_id','res'),sep='_Ten')

  res<-pfilt|>
    dplyr::left_join(samples)|>
    dplyr::left_join(genes)|>
    dplyr::select(improve_sample_id,entrez_id,proteomics)|>
    dplyr::distinct()
  write_csv(res,file=gzfile('proteomics.csv.gz'))
}


mirnaFixing<-function(mirlist){

  ##first let's get the prefix off
  trimmed<-lapply(mirlist,function(x)
    stringr::str_replace(x,'hsa-','')|>
      stringr::str_replace_all('-','')|>
      toupper())

  newmap<-data.frame(old=mirlist,gene_symbol=unlist(trimmed))
  return(newmap)
}


###run through each file and rewrite
newres<-lapply(names(filenames),function(value){

  fi=paste0(basename,filenames[[value]])
  fname=paste0(value,'.csv.gz')
  print(paste('now reading',fi,'to store as',fname))
  ##now every data type is parsed slightly differently, so we need to change our formatting
  ##and mapping to get it into a unified 3 column schema
  if(value=='copy_number'){
    exp_file <- readr::read_csv(fi,skip=1)[-1,]

    res = exp_file|>
      tidyr::pivot_longer(cols=c(2:ncol(exp_file)),
                          names_to='gene_symbol',values_to='copy_number',
                          values_transform=list(copy_number=as.numeric))|>
      dplyr::distinct()

    res<-res|> ##deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
      dplyr::mutate(copy_call=ifelse(copy_number<0.5210507,'deep del',
                                     ifelse(copy_number<0.7311832,'het loss',
                                            ifelse(copy_number<1.214125,'diploid',
                                                   ifelse(copy_number<1.422233,'gain','amp')))))|>
      dplyr::left_join(genes)|>
      dplyr::distinct()


    colnames(res)[1]<-'other_id'
    vars=c('copy_number','copy_call')


  }else if(value=='methylation'){ ###IF DATA REPRESENT RRBS###
    exp_file <- readr::read_csv(fi)[-c(1:2),]
    #the gene names are not unique and in some weird format, i willt ry to keep both in the metadata
    res = exp_file|>
      tidyr::pivot_longer(cols=c(2:ncol(exp_file)),
                              names_to='gene_region',values_to='methylation',
                          values_transform=list(methylation=as.numeric))|>
      dplyr::distinct()

    res<-res|>
      tidyr::separate(gene_region,into=c('gene_symbol','num','start','end'),sep='_')|>
      dplyr::left_join(genes)|>
      dplyr::distinct()

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
        dres<-exp_file[,c('Entrez_id','Genome_Change','Variant_Classification',colnames(exp_file)[cols])]
        ret <-tidyr::pivot_longer(dres,cols=c(4:ncol(dres)),names_to='other_id',values_to='num_muts')
        #print(head(ret))
        ret|>
          subset(num_muts!=0)
      })

      res<-do.call(rbind,res)
      full<-res|>  ###since we're already in ENTREZ we skip the mapping below
        dplyr::left_join(samples)|>
        dplyr::rename(entrez_id=Entrez_id,mutations=Genome_Change,variant_classification=Variant_Classification)|>
        dplyr::select(entrez_id,improve_sample_id,mutations,variant_classification)|>
         dplyr::mutate(source='DepMap',study='CCLE')|>
        dplyr::distinct()

        write_csv(full,file=fname)
        return(fi)
    }
    else if(value=='transcriptomics'){ #if gene expression
      exp_file <- readr::read_csv(fi,skip=2)

      res = tidyr::pivot_longer(data=exp_file,cols=c(2:ncol(exp_file)),
                                names_to='gene_symbol',values_to='transcriptomics',
                                values_transform=list(expression=as.numeric))|>
        dplyr::left_join(genes)|>
        dplyr::distinct()
      colnames(res)[1]<-'other_id'
      vars=c('transcriptomics')

    }
    else if(value=='miRNA'){ #if mirna expression
      exp_file <- readr::read_csv(fi)

      res = tidyr::pivot_longer(data=exp_file,cols=c(2:ncol(exp_file)),
                                names_to='gene_symbol',values_to='miRNA',values_transform=list(mirnas=as.numeric))

      gmod<-genes
      gmod$gene_symbol<-tolower(gmod$gene_symbol)

      res$gene_symbol<-tolower(res$gene_symbol)

      res<-res|>
        dplyr::left_join(gmod)|>
        dplyr::distinct()
        #dplyr::mutate(gene_symbol=tolower(gene_symbol))

      missed<-subset(res,is.na(entrez_id))
      notmissed<-subset(res,!is.na(entrez_id))
      print(paste("matched",length(unique(notmissed$gene_symbol)),'miRNAs and missed',length(unique(missed$gene_symbol))))

      ##mirfix
      mirmap<-mirnaFixing(unique(missed$gene_symbol))
      fixed<-missed|>
        dplyr::rename(old='gene_symbol')|>
        dplyr::select(-entrez_id)|>
        left_join(mirmap)|>
        left_join(genes)
      notmissed<-fixed|>
        subset(!is.na(entrez_id))|>
        dplyr::select(-old)|>
        rbind(notmissed)

      missed<-subset(fixed,is.na(entrez_id))
      print(paste("after second pass, matched",length(unique(notmissed$gene_symbol)),'miRNAs and missed',length(unique(missed$gene_symbol))))

      colnames(res)[1]<-'other_id'
      vars=c('miRNA')

    }

  ##do the last join with samples
  full<-res|>
    dplyr::left_join(samples)|>
    dplyr::select(c('entrez_id','improve_sample_id',vars))|>
    dplyr::distinct()|>
    dplyr::mutate(source='DepMap',study='CCLE')

  write_csv(full,file=gzfile(fname))
  return(fi)

})

