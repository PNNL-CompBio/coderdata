##getSanger omics data

library(readr)
library(tidyr)
library(dplyr)

Sys.setenv(VROOM_CONNECTION_SIZE=100000000)


#basename='https://ftp.mcs.anl.gov/pub/candle/public/improve/Data/Omics/Curated_CCLE_Multiomics_files/'
filenames=list(transcriptomics='https://cog.sanger.ac.uk/cmp/download/rnaseq_all_20220624.zip',
               copy_number='https://cog.sanger.ac.uk/cmp/download/WES_pureCN_CNV_genes_latest.csv.gz',
               proteomics='https://cog.sanger.ac.uk/cmp/download/Proteomics_20221214.zip',
               mutation='https://cog.sanger.ac.uk/cmp/download/mutations_all_20230202.zip')
# the dictionary started with
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
colnames(vtab)<-c('variant_classification','effect')

getAll<-function(dt=names(filenames)){

  options(timeout=10000)
  ###run through each file and rewrite
  newres<-lapply(dt,function(value){

    fi=filenames[[value]]
     fname=paste0('/tmp/sanger_',value,'.csv.gz')
    print(paste('now reading',fi,'to store as',fname))
    ##now every data type is parsed slightly differently, so we need to change our formatting
    ##and mapping to get it into a unified 3 column schema
    if(value=='copy_number'){
      #read in file
      exp_file <- readr::read_csv(fi) ##already in long form <3 <3 <3

      smap<-samples|>
          subset(other_id_source=='Sanger')|>
          subset(other_id%in%exp_file$model_id)|>
          dplyr::select(improve_sample_id,other_id)|>
          distinct()

        gmap<-genes|>
            subset(gene_symbol%in%exp_file$symbol)|>
            distinct()

        print('wide to long')

      res<-exp_file|>
        dplyr::select(other_id='model_id',gene_symbol='symbol',gatk_mean_log2_copy_ratio,source,data_type,cn_category)|>
        mutate(copy_number=2^gatk_mean_log2_copy_ratio,.keep='all')|>
        distinct()|>
        left_join(gmap)|>
        dplyr::select(other_id,source,copy_number,entrez_id,sanger_copy_call='cn_category')|>
        left_join(smap)|>
          distinct()
       rm(exp_file)

        print('copy call')
      ##calibrate the copy call
      res<-res|> ##deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
        dplyr::mutate(improve_copy_call=ifelse(copy_number<0.5210507,'deep del',
                                       ifelse(copy_number<0.7311832,'het loss',
                                              ifelse(copy_number<1.214125,'diploid',
                                                     ifelse(copy_number<1.422233,'gain','amp')))))|>
        dplyr::distinct()|>
        mutate(study='Sanger')

      full<-res|>
        tidyr::pivot_longer(cols=c(improve_copy_call,sanger_copy_call),
                            names_to='copy_call_source',
                            values_to='copy_call')
      rm(res)
#      full<-lres

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


    }else if(value=='mutation'){ ####IF DATA REPRESENTS MUTATIONS#####
      res=download.file(fi,'/tmp/tmp.zip')
      filist<-unzip('/tmp/tmp.zip',exdir='/tmp')
      fi= "/tmp/mutations_all_20230202.csv"

      exp_file <- readr::read_csv(fi)|>
        dplyr::select(symbol='gene_symbol',other_id='model_id',effect,mutation='cdna_mutation',source)|>
        distinct()

      smap<-samples|>
        dplyr::select(improve_sample_id,other_id)|>distinct()

      res<-exp_file|>
        left_join(smap)|>
          mutate(study='Sanger')|>
          left_join(as.data.frame(vtab))|>
          dplyr::select(-effect)|>
          subset(!is.na(improve_sample_id))|>
          distinct()

      rm(exp_file)

#      res$variant_classification=unlist(lapply(res$effect,function(x) names(variant_schema)[grep(x,variant_schema)]))
#      res<-res|>dplyr::select(-effect)
      write_csv(res,file=fname)
      rm(res)
      return(fi)
    }
    else if(value=='transcriptomics'){ #if gene expression
      res=download.file(fi,'/tmp/tmp.zip')
      filist<-unzip('/tmp/tmp.zip',exdir='/tmp')
      fi= "/tmp/rnaseq_tpm_20220624.csv"

      exp_file <- readr::read_csv(fi)
      ##the rows have metadata
      samps<-t(exp_file[1:3,])
      colnames(samps)<-samps[1,]
      samps<-samps[-c(1:2),]|>as.data.frame()|>
        tibble::rownames_to_column('other_id')|>
        left_join(samples)|>
        dplyr::rename(source='data_source',study='dataset_name')

      missing<-subset(samps,is.na(improve_sample_id))|>
        dplyr::select(-c(other_id,improve_sample_id))|>
        dplyr::rename(other_id='model_name')

      ##the data starts at row 5
      dat<-exp_file[-c(1:4),-1]
      colnames(dat)[1]<-'gene_symbol'

      rm(exp_file)

      ddat<-apply(dat,2,unlist)|>
          as.data.frame()|>
#      ddat<-ddat|>
        left_join(genes)|>
        dplyr::select(-gene_symbol)

      rm(dat)
      ddat$entrez_id<-as.numeric(ddat$entrez_id)
      res = tidyr::pivot_longer(data=as.data.frame(ddat),cols=c(1:(ncol(ddat)-1)),
                                names_to='other_id',values_to='transcriptomics',
                                values_transform=list(expression=as.numeric))|>
          distinct()

      rm(ddat)
      smap<-samps|>
        dplyr::select(improve_sample_id,other_id,study,source)|>distinct()

      full<-res|>
        left_join(smap)
      rm(res)
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
      full<-res

    }else if(value=='proteomics'){
      res=download.file(fi,'/tmp/tmp.zip')
      filist<-unzip('/tmp/tmp.zip',exdir='/tmp')

      fi='/tmp/Protein_matrix_averaged_zscore_20221214.tsv'
      exp_file <- readr::read_tsv(fi,skip=1)[-1,-1]
      colnames(exp_file)[1]<-'other_id'

      smap<-samples|>
        dplyr::select(improve_sample_id,other_id)|>distinct()

      res<-exp_file|>
        tidyr::pivot_longer(cols=(-c(other_id)),names_to='gene_symbol',values_to='proteomics')|>
        subset(!is.na(proteomics))|>
        left_join(genes)|>
        dplyr::select(-gene_symbol)|>
        left_join(smap)|>
          mutate(source='Sanger',study='Sanger')
      rm(exp_file)
      full<-res
      rm(res)

    }

    ##do the last join with samples
    #full<-res#|>
    #  left_join(samples)
    missed<-full|>subset(is.na(improve_sample_id))|>
      dplyr::select(improve_sample_id,other_id)|>
      distinct()
    print(paste('missing',nrow(missed),'identifiers'))
    print(missed)

    full<-full|>
      subset(!is.na(improve_sample_id))|>
      dplyr::select(-other_id)

    write_csv(full,file=gzfile(fname))
      rm(full)
      return(fi)

  })

}


main<-function(){
	args = commandArgs(trailingOnly=TRUE)
	if(length(args)!=2){
	  print('Usage: Rscript 02a-pullSanger.R [genefile] [samplefile]')
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
		     dplyr::select(other_id,improve_sample_id,other_id_source)|>
		       unique()
	getAll()

}

main()
