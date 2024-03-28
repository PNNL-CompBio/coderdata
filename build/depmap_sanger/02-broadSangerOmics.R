# coderdata cell line omics data capture
# pulls down static versions of depmap and sanger cell line omics measurements
# concatenates each into sinle file
library(readr)
library(tidyr)
library(dplyr)
library(rio)

Sys.setenv(VROOM_CONNECTION_SIZE=100000000)

##### DEPMAP FILES

depmap_filenames=list(   copy_number='https://figshare.com/ndownloader/files/40448840',
               transcriptomics='https://figshare.com/ndownloader/files/40449128',
                              mutations='https://figshare.com/ndownloader/files/40449638')
##### SANGER FILES
sanger_filenames=list(transcriptomics='https://cog.sanger.ac.uk/cmp/download/rnaseq_all_20220624.zip',
               copy_number='https://cog.sanger.ac.uk/cmp/download/WES_pureCN_CNV_genes_latest.csv.gz',
               mutations='https://cog.sanger.ac.uk/cmp/download/mutations_all_20230202.zip')


###### VARIANT SCHEMA HARMONIZATION
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

depmap_vtab<-do.call('rbind',sapply(names(variant_schema),function(x) cbind(rep(x,length(variant_schema[[x]])),unlist(variant_schema[[x]]))))
colnames(depmap_vtab)<-c('variant_classification','VariantInfo')

sanger_vtab<-do.call('rbind',sapply(names(variant_schema),function(x) cbind(rep(x,length(variant_schema[[x]])),unlist(variant_schema[[x]]))))
colnames(sanger_vtab)<-c('variant_classification','effect')


mirnaFixing<-function(mirlist){

  ##first let's get the prefix off
  trimmed<-lapply(mirlist,function(x)
    stringr::str_replace(x,'hsa-','')|>
      stringr::str_replace_all('-','')|>
      toupper())

  newmap<-data.frame(old=mirlist,gene_symbol=unlist(trimmed))
  return(newmap)
}


#### pull down ad process each of the sanger files
sanger_files<-function(fi,value){

  options(timeout=10000)

#  newres<-lapply(dt,function(value){

 #   fi=sanger_filenames[[value]]
    # fname=paste0('/tmp/sanger_',value,'.csv.gz')
    #print(paste('now reading',fi,'to store as',fname))
    ##now every data type is parsed slightly differently, so we need to change our formatting
    ##and mapping to get it into a unified 3 column schema
    if(value=='copy_number'){
      #read in file
      exp_file <- readr::read_csv(fi) ##already in long form <3 <3 <3
      file.remove(fi)
      smap<-sanger_samples|>
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
        dplyr::select(other_id,copy_number,entrez_id,Sanger='cn_category')|>
        left_join(smap)|>
          distinct()
       rm(exp_file)

        print('copy call')
      ##calibrate the copy call
      res<-res|> ##deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
        dplyr::mutate(IMPROVE=ifelse(copy_number<0.5210507,'deep del',
                                       ifelse(copy_number<0.7311832,'het loss',
                                              ifelse(copy_number<1.214125,'diploid',
                                                     ifelse(copy_number<1.422233,'gain','amp')))))|>
        dplyr::distinct()|>
        mutate(study='Sanger')

      full<-res|>
        tidyr::pivot_longer(cols=c(IMPROVE,Sanger),
                            names_to='source',
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


    }else if(value=='mutations'){ ####IF DATA REPRESENTS MUTATIONS#####
      res=download.file(fi,'/tmp/tmp.zip')
      filist<-unzip('/tmp/tmp.zip',exdir='/tmp')
      fi= "/tmp/mutations_all_20230202.csv"
      if(file.exists("/tmp/tmp.zip"))
          file.remove('/tmp/tmp.zip')

      exp_file <- readr::read_csv(fi)|>
        dplyr::select(gene_symbol,other_id='model_id',effect,mutation='cdna_mutation',source)|>
          distinct()
      if(file.exists(fi))
          file.remove(fi)

      smap<-sanger_samples|>
        dplyr::select(improve_sample_id,other_id)|>distinct()

      res<-exp_file|>
          left_join(genes)|>
        left_join(smap)|>
          mutate(study='Sanger')|>
          dplyr::select(-c(other_id,gene_symbol))|>
          left_join(as.data.frame(sanger_vtab))|>
          dplyr::select(-effect)|>
          subset(!is.na(improve_sample_id))|>
          distinct()

      rm(exp_file)

#      res$variant_classification=unlist(lapply(res$effect,function(x) names(variant_schema)[grep(x,variant_schema)]))
#      res<-res|>dplyr::select(-effect)
      #write_csv(res,file=fname)
      #rm(res)
                                        #return(fi)
      print(head(res))
      return(res)
    }else if(value=='transcriptomics'){ #if gene expression
      res=download.file(fi,'/tmp/tmp.zip')
      filist<-unzip('/tmp/tmp.zip',exdir='/tmp')
      fi= "/tmp/rnaseq_tpm_20220624.csv"
      if(file.exists("/tmp/tmp.zip"))
         file.remove('/tmp/tmp.zip')
         exp_file <- readr::read_csv(fi)

      ##the rows have metadata
      samps<-t(exp_file[1:3,])
      colnames(samps)<-samps[1,]
      samps<-samps[-c(1:2),]|>as.data.frame()|>
        tibble::rownames_to_column('other_id')|>
        left_join(sanger_samples)|>
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
      file.remove(fi)
    }else if(value=='miRNA'){ #if mirna expression
      exp_file <- readr::read_csv(fi)

      res = tidyr::pivot_longer(data=exp_file,cols=c(2:ncol(exp_file)),
                                names_to='gene_symbol',values_to='miRNA',values_transform=list(mirnas=as.numeric))
      file.remove(fi)
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
      file.remove(fi)
      smap<-sanger_samples|>
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

    #write_csv(full,file=gzfile(fname))
#      rm(full)
      print(paste('Sanger',value,':'))
      print(head(full))
      return(full)
  #  })
  #    names(newres)<-names(sanger_filenames)

}


#### pull dodwn and process each depmap file per the 2023 q4 schema
depmap_files<-function(fi,value){

   ##runs through entire list of omcis types and files and returns  a list
#  newres<-lapply(values,function(value){

#    fi=depmap_filenames[[value]]
#    fname=paste0('/tmp/depmap_',value,'.csv.gz')
#    print(paste('now reading',fi,'to store as',fname))
    ##now every data type is parsed slightly differently, so we need to change our formatting
    ##and mapping to get it into a unified 3 column schema
    if(value=='copy_number'){
      exp_file <- readr::read_csv(fi)

      print('Long to wide')
      res = exp_file|>
        tidyr::pivot_longer(cols=c(2:ncol(exp_file)),
                            names_to='gene_entrez',values_to='copy_number',
                            values_transform=list(copy_number=as.numeric))|>
        dplyr::distinct()
      rm(exp_file)

      colnames(res)[1]<-'other_id'

      print('String manipulations')
      res<-res|>
          tidyr::separate_wider_delim(gene_entrez,' ',names=c('gene_symbol','entrez_id'))

      print('join with gene')
      res<-res|>
          dplyr::select(-entrez_id)|>
          left_join(genes)|>
          dplyr::select(other_id,entrez_id,copy_number)|>
          distinct()
      ##these are messing things up
    #  res$entrez_id<-stringr::str_replace(res$entrez_id,'\\)','')
    #  res$entrez_id<-stringr::str_replace(res$entrez_id,'\\(','')

      print('Adding copy calls')
      res<-res|> ##deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
        dplyr::mutate(copy_call=ifelse(copy_number<0.5210507,'deep del',
                                       ifelse(copy_number<0.7311832,'het loss',
                                              ifelse(copy_number<1.214125,'diploid',
                                                     ifelse(copy_number<1.422233,'gain','amp')))))|>
          dplyr::distinct()|>
          subset(!is.na(copy_number))


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
            left_join(as.data.frame(depmap_vtab))|>
            dplyr::select(-c(EntrezGeneID,VariantInfo))|>
            distinct()|>
          subset(!is.na(entrez_id)) ##removes thos with unknonw entrez

        rm(exp_file)

        smap<-depmap_samples|>
            dplyr::select(improve_sample_id,other_id)|>
            distinct()|>
            subset(other_id%in%res$other_id)

        full<-res|>  ###since we're already in ENTREZ we skip the mapping below
          dplyr::left_join(smap)|>
          #dplyr::rename(entrez_id=Entrez_id,mutations=Genome_Change,variant_classification=Variant_Classification)|>
          dplyr::select(entrez_id,improve_sample_id,mutation,variant_classification)|>
           dplyr::mutate(source='Broad',study='DepMap')|>
          dplyr::distinct()

          #write_csv(full,file=fname)
                                        #return(fi)
        print(head(full))
        return(full)
      }else if(value=='transcriptomics'){ #if gene expression
        exp_file <- readr::read_csv(fi)
        print("wide to long")
        res = tidyr::pivot_longer(data=exp_file,cols=c(2:ncol(exp_file)),
                                  names_to='gene_entrez',values_to='transcriptomics',
                                  values_transform=list(expression=as.numeric))
        colnames(res)[1]<-'other_id'

        print('fixing gene names')
        res<-res|>
          tidyr::separate_wider_delim(gene_entrez,' ',names=c('gene_symbol','entrez_par'))

      print('join with gene')
      res<-res|>
          dplyr::select(-entrez_par)|>
          left_join(genes)|>
          dplyr::select(other_id,entrez_id,transcriptomics)|>
          distinct()

    #mutate(entrez_id=stringr::str_replace_all(entrez_par,'\\)|\\(',''))|>
    #      dplyr::select(-c(entrez_par,gene))|>
    #      distinct()
        rm(exp_file)
        vars=c('transcriptomics')

      }else if(value=='miRNA'){ #if mirna expression
        exp_file <- readr::read_csv(fi)

        res = tidyr::pivot_longer(data=exp_file,cols=c(2:ncol(exp_file)),
                                  names_to='gene_symbol',values_to='miRNA',values_transform=list(mirnas=as.numeric))

        rm(exp_file)
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
      smap<-depmap_samples|>
          subset(other_id%in%res$other_id)|>
          distinct()

      print('joining with samples')
      full<-res|>
          dplyr::left_join(smap)
      rm(res)

      missed<-full|>subset(is.na(improve_sample_id))|>
          dplyr::select(improve_sample_id,other_id)|>
          distinct()
      print(paste('missing',nrow(missed),'identifiers'))
      print(missed)

      full<-full|>dplyr::select(c('entrez_id','improve_sample_id',vars))|>
          subset(!is.na(improve_sample_id))|>
          dplyr::distinct()|>
          dplyr::mutate(source='Broad',study='DepMap')

      #write_csv(full,file=gzfile(fname))
      #rm(full)
      print(paste('Depmap',value,':'))
      print(head(full))
      return(full)


}


main<-function(){
    args = commandArgs(trailingOnly=TRUE)
    if(length(args)!=2){
        stop('Usage: Rscript 02-depmap-sanger-omics.R [genefile] [samplefile]')
    }
    gfile = args[1]
    sfile = args[2]
###first reqad in all gene information so we can map appropriately
    allgenes = read_csv(gfile)

    genes <<- allgenes|>
        dplyr::select(gene_symbol,entrez_id)|>
        dplyr::distinct()

    ##here are the improve sample id indices
    depmap_samples <<- read_csv(sfile,
                               quote='"')|>
        dplyr::select(other_id,improve_sample_id)|>
        unique()

    sanger_samples <<- read_csv(sfile,
                                quote='"')|>
        dplyr::select(other_id,improve_sample_id,other_id_source)|>
        unique()

    alltypes<-c('mutations','transcriptomics','copy_number')

    lapply(alltypes,function(dt){
        print(dt)
        temps<-sanger_files(sanger_filenames[[dt]],dt)
        tempd<-depmap_files(depmap_filenames[[dt]],dt)
        readr::write_csv(rbind(tempd,temps),file=paste0('/tmp/depmap_sanger_',dt,'.csv'))
        rm(tempd)
        rm(temps)
    })
    system(paste0('/opt/venv/bin/python 02a-depmap_sanger_proteomics.py --gene ',gfile,' --sample ',sfile))

}
main()
