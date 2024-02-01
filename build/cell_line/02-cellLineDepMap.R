
library(readr)
library(tidyr)
library(dplyr)
library(rio)

Sys.setenv(VROOM_CONNECTION_SIZE=100000000)
###PATH TO depmap 23q2 on FigShare
#basename="https://figshare.com/articles/dataset/DepMap_23Q2_Public/22765112"

#basename='https://ftp.mcs.anl.gov/pub/candle/public/improve/Data/Omics/Curated_CCLE_Multiomics_files/'
filenames=list(transcriptomics='https://figshare.com/ndownloader/files/40449128',
                            copy_number='https://figshare.com/ndownloader/files/40448840',
                            mutations='https://figshare.com/ndownloader/files/40449638')


# the dictionary started with
# CCLE data
# CCLE data
variant_schema =list(`3'UTR`=c("3'UTR",'THREE_PRIME_UTR','3prime_UTR_variant','3prime_UTR_ess_splice'),
                     `5'Flank`=c("FIVE_PRIME_FLANK","5'Flank",'upstream'),
                     `5'UTR`=c("5'UTR",'5prime_UTR_variant','5prime_UTR_variant','5prime_UTR_ess_splice'),
                     Undetermined=c('COULD_NOT_DETERMINE'),
                     De_novo_Start_InFrame=c('DE_NOVO_START_IN_FRAME','De_novo_Start_InFrame'),
                     De_novo_Start_OutOfFrame=c('DE_NOVO_START_OUT_FRAME','De_novo_Start_OutOfFrame'),
                     Frame_Shift_Del=c('FRAME_SHIFT_DEL','Frame_Shift_Del','frameshift'),
                     Frame_Shift_Ins=c('FRAME_SHIFT_INS','Frame_Shift_Ins'),
                     IGR=c('IGR','nc_variant'),
                     In_Frame_Del=c('IN_FRAME_DEL','In_Frame_Del','inframe'),
                     In_Frame_Ins=c('IN_FRAME_INS','In_Frame_Ins'),
                     Intron=c('INTRON','Intron','intronic'),
                     Missense_Mutation=c('Missense_Mutation','MISSENSE','missense'),
                     Nonsense_Mutation=c('Nonsense_Mutation','NONSENSE','nonsense'),
                     Nonstop_Mutation=c('Nonstop_Mutation','NONSTOP'),
                     RNA=c('RNA'),
                     Start_Codon_SNP=c('START_CODON_SNP','Start_Codon_SNP'),
                     Start_Codon_Del=c('Start_Codon_Del','START_CODON_DEL','start_lost'),
                     Start_Codon_Ins=c('Start_Codon_Ins','START_CODON_INS'),
                     Stop_Codon_Del=c('Stop_Codon_Del','stop_lost'),
                     Stop_Codon_Ins=c('Stop_Codon_Ins'),
                     Silent=c('Silent','SILENT','silent'),
                     Splice_Site=c('Splice_Site','SPLICE_SITE','splice_region'),
                     Translation_Start_Site=c('Translation_Start_Site','start_lost'))

vtab<-do.call('rbind',sapply(names(variant_schema),function(x) cbind(rep(x,length(variant_schema[[x]])),unlist(variant_schema[[x]]))))
colnames(vtab)<-c('variant_classification','VariantInfo')

getProteomics<-function(){
  #pull directly from gygi lab
  proteomics <- 'https://gygi.hms.harvard.edu/data/ccle/Table_S2_Protein_Quant_Normalized.xlsx'
  options(timeout=300)

#  res<-download.file(proteomics,'prot.xlsx')
  pdat<-rio::import(proteomics,which=2)#readxl::read_xlsx('prot.xlsx',sheet = 'Normalized Protein Expression')
  pdat[,7:ncol(pdat)]<-apply(pdat[,7:ncol(pdat)],2,as.numeric)
  pdat<-pdat|>
    dplyr::select(!starts_with('Ten'))
  plong<-pdat|>
    tidyr::pivot_longer(cols=7:ncol(pdat),names_to='cellLine',values_to='proteomics',values_drop_na=TRUE)

  pfilt<-plong|>
    dplyr::select(gene_symbol='Gene_Symbol',cellLine,proteomics)|>
    dplyr::distinct()|>
    tidyr::separate(cellLine,into=c('other_id','res'),sep='_Ten')

    smap<-samples|>
        subset(other_id%in%pfilt$other_id)|>
        distinct()

    gmap<-genes|>
        subset(gene_symbol%in%pfilt$gene_symbol)|>
        dplyr::select(gene_symbol,entrez_id)|>
        distinct()

  res<-pfilt|>
    dplyr::left_join(smap)|>
    dplyr::left_join(gmap)|>
      dplyr::select(improve_sample_id,entrez_id,proteomics)|>
      subset(!is.na(entrez_id))|>
      dplyr::distinct()
    res$study='DepMap'
    res$source='Broad'
  write_csv(res,file=gzfile('/tmp/proteomics.csv.gz'))
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


do_all<-function(values=names(filenames)){
###run through each file and rewrite
  newres<-lapply(values,function(value){

    fi=filenames[[value]]
    fname=paste0('/tmp/depmap_',value,'.csv.gz')
    print(paste('now reading',fi,'to store as',fname))
    ##now every data type is parsed slightly differently, so we need to change our formatting
    ##and mapping to get it into a unified 3 column schema
    if(value=='copy_number'){
      exp_file <- readr::read_csv(fi)

      res = exp_file|>
        tidyr::pivot_longer(cols=c(2:ncol(exp_file)),
                            names_to='gene_entrez',values_to='copy_number',
                            values_transform=list(copy_number=as.numeric))|>
        dplyr::distinct()

      res<-res|>
        tidyr::separate_wider_delim(gene_entrez,' ',names=c('gene','entrez_id'))|>
        dplyr::select(-gene)
      res$entrez_id<-stringr::str_replace(res$entrez_id,'\\)','')
      res$entrez_id<-stringr::str_replace(res$entrez_id,'\\(','')

      res<-res|> ##deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
        dplyr::mutate(copy_call=ifelse(copy_number<0.5210507,'deep del',
                                       ifelse(copy_number<0.7311832,'het loss',
                                              ifelse(copy_number<1.214125,'diploid',
                                                     ifelse(copy_number<1.422233,'gain','amp')))))|>
      #  dplyr::left_join(genes)|>
          dplyr::distinct()|>
          subset(!is.na(copy_number))


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
        exp_file <- readr::read_csv(fi)|>
          dplyr::select(EntrezGeneID,HgncName,other_id='ModelID',VariantInfo,mutation='DNAChange')|>
          distinct()

        res<-exp_file|>
          mutate(entrez_id=as.numeric(EntrezGeneID))|>
            left_join(as.data.frame(vtab))|>
            dplyr::select(-c(EntrezGeneID,VariantInfo))|>
          subset(!is.na(entrez_id)) ##removes thos with unknonw entrez


        smap<-samples|>
            dplyr::select(improve_sample_id,other_id)|>distinct()


        #res$variant_classification=unlist(lapply(res$VariantInfo,function(x) names(variant_schema)[grep(x,variant_schema)]))
        full<-res|>  ###since we're already in ENTREZ we skip the mapping below
          dplyr::left_join(smap)|>
          #dplyr::rename(entrez_id=Entrez_id,mutations=Genome_Change,variant_classification=Variant_Classification)|>
          dplyr::select(entrez_id,improve_sample_id,mutation,variant_classification)|>
           dplyr::mutate(source='DepMap',study='CCLE')|>
          dplyr::distinct()

          write_csv(full,file=fname)
          return(fi)
      }
      else if(value=='transcriptomics'){ #if gene expression
        exp_file <- readr::read_csv(fi)

        res = tidyr::pivot_longer(data=exp_file,cols=c(2:ncol(exp_file)),
                                  names_to='gene_entrez',values_to='transcriptomics',
                                  values_transform=list(expression=as.numeric))|>
          tidyr::separate_wider_delim(gene_entrez,' ',names=c('gene','entrez_par'))|>
          mutate(entrez_id=stringr::str_replace_all(entrez_par,'\\)|\\(',''))|>
          dplyr::select(-c(entrez_par,gene))|>
          distinct()

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
      dplyr::left_join(samples)

    missed<-full|>subset(is.na(improve_sample_id))|>
      dplyr::select(improve_sample_id,other_id)|>
      distinct()
    print(paste('missing',nrow(missed),'identifiers'))
    print(missed)

    full<-full|>dplyr::select(c('entrez_id','improve_sample_id',vars))|>
      subset(!is.na(improve_sample_id))|>
      dplyr::distinct()|>
      dplyr::mutate(source='DepMap',study='Broad')

    write_csv(full,file=gzfile(fname))
    return(fi)

  })
}


main<-function(){
	args = commandArgs(trailingOnly=TRUE)
	if(length(args)!=2){
	  print('Usage: Rscript 02-cellLineDepMap.R [genefile] [samplefile]')
	  exit()
	  }
	gfile = args[1]
	sfile = args[2]
	###first reqad in all gene information so we can map appropriately
	allgenes = read_csv(gfile)
	genes <<- allgenes|>
	  dplyr::select(gene_symbol,entrez_id)|>
	    dplyr::distinct()

	    ##here are the improve sample id indices
	 samples <<- read_csv(sfile,
                   quote='"')|>
		     dplyr::select(other_id,improve_sample_id)|>
		       unique()
	do_all(names(filenames))
	getProteomics()

}

main()
