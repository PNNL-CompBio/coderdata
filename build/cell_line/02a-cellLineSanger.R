##getSanger omics data

library(readr)
library(tidyr)
library(dplyr)
###first reqad in all gene information so we can map appropriately
allgenes = read_csv("../genes.csv")
genes = allgenes|>
  dplyr::select(gene_symbol,entrez_id)|>
  dplyr::distinct()



##here are the improve sample id indices
samples = read_csv('samples.csv',
                   quote='"')|>
  dplyr::select(other_id,improve_sample_id)|>
  unique()


#basename='https://ftp.mcs.anl.gov/pub/candle/public/improve/Data/Omics/Curated_CCLE_Multiomics_files/'
filenames=list(transcriptomics='https://cog.sanger.ac.uk/cmp/download/rnaseq_all_20220624.zip',
               copy_number='https://cog.sanger.ac.uk/cmp/download/WES_pureCN_CNV_segments_20221214.zip',
               #methylation='CCLE_AID_RRBS_TSS_1kb_20180614.csv',
               proteomics='https://cog.sanger.ac.uk/cmp/download/Proteomics_20221214.zip',
               #miRNA='CCLE_AID_miRNA_20180525.csv',
               mutations='https://cog.sanger.ac.uk/cmp/download/mutations_all_20230202.zip')
# the dictionary started with
# CCLE data
variant_schema =list(`3'UTR`=c("3'UTR",'THREE_PRIME_UTR'),
                     `5'Flank`=c("FIVE_PRIME_FLANK","5'Flank"),`5'UTR`=c("5'UTR"),
                     Undetermined=c('COULD_NOT_DETERMINE'),
                     De_novo_Start_InFrame=c('DE_NOVO_START_IN_FRAME','De_novo_Start_InFrame'),
                     De_novo_Start_OutOfFrame=c('DE_NOVO_START_OUT_FRAME','De_novo_Start_OutOfFrame'),
                     Frame_Shift_Del=c('FRAME_SHIFT_DEL','Frame_Shift_Del'),
                     Frame_Shift_Ins=c('FRAME_SHIFT_INS','Frame_Shift_Ins'),
                     IGR=c('IGR'),In_Frame_Del=c('IN_FRAME_DEL','In_Frame_Del'),
                     In_Frame_Ins=c('IN_FRAME_INS','In_Frame_Ins'),Intron=c('INTRON','Intron'),
                     Missense_Mutation=c('Missense_Mutation','MISSENSE'),
                     Nonsense_Mutation=c('Nonsense_Mutation','NONSENSE'),
                     Nonstop_Mutation=c('Nonstop_Mutation','NONSTOP'),
                     RNA=c('RNA'),
                     Start_Codon_SNP=c('START_CODON_SNP','Start_Codon_SNP'),
                     Start_Codon_Del=c('Start_Codon_Del','STARD_CODON_DEL'),
                     Start_Codon_Ins=c('Start_Codon_Ins','START_CODON_INS'),
                     Stop_Codon_Del=c('Stop_Codon_Del'),
                     Stop_Codon_Ins=c('Stop_Codon_Ins'),
                     Silent=c('Silent','SILENT'),
                     Splice_Site=c('Splice_Site','SPLICE_SITE'),
                     Translation_Start_Site=c('Translation_Start_Site'))

getAll<-function(dt=names(filenames)){
  
  ###run through each file and rewrite
  newres<-lapply(df,function(value){
    
    fi=filenames[[value]]
    fname=paste0(value,'.csv.gz')
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
      exp_file <- readr::read_csv(fi)|>
        dplyr::select(EntrezGeneID,HgncName,other_id='ModelID',VariantInfo,mutations='DNAChange')|>
        distinct()
      
      res<-exp_file|>
        mutate(entrez_id=as.numeric(EntrezGeneID))|>
        dplyr::select(-EntrezGeneID)|>
        subset(!is.na(entrez_id)) ##removes thos with unknonw entrez
      
      res$variant_classification=unlist(lapply(res$VariantInfo,function(x) names(variant_schema)[grep(x,variant_schema)]))
      full<-res|>  ###since we're already in ENTREZ we skip the mapping below
        dplyr::left_join(samples)|>
        #dplyr::rename(entrez_id=Entrez_id,mutations=Genome_Change,variant_classification=Variant_Classification)|>
        dplyr::select(entrez_id,improve_sample_id,mutations,variant_classification)|>
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
      dplyr::left_join(samples)|>
      dplyr::select(c('entrez_id','improve_sample_id',vars))|>
      dplyr::distinct()|>
      dplyr::mutate(source='DepMap',study='CCLE')
    
    write_csv(full,file=gzfile(fname))
    return(fi)
    
  })
  
}