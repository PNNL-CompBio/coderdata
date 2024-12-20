# Load required libraries
library(data.table)
# library(biomaRt)# biomart issues still exist
library(synapser)
library(dplyr)

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if a token was provided
if (length(args) == 0) {
  stop("No token or sample file provided. Usage: Rscript my_script.R <PAT> [samples] [genes]", call. = FALSE)
}

# Set your personal access token
PAT <- args[1]
patients <- args[2]
genefile <- args[3]

# Log in to Synapse
synLogin(authToken = PAT)

# Define the Ensembl mart # biomart issues still exist
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # biomart issues still exist; fix later...

# Path to the directory to save .sf files
#path <- "./tmp"
#dir.create(path, showWarnings = FALSE)

# Read the sample mapping CSV and genes.csv
samples_df <- fread(patients)|>
    dplyr::select(improve_sample_id,common_name,model_type)|>
                                        distinct()#"mpnst/synapse_NF-MPNST_samples.csv")

pdx_samps<-subset(samples_df,model_type=='patient derived xenograft')
tumor_samps<-subset(samples_df,model_type=='tumor')
mt_samps<-subset(samples_df,model_type=='organoid')

##now get the manifest from synapse
manifest<-synapser::synTableQuery("select * from syn53503360")$asDataFrame()|>
                                                             as.data.frame()|>
                                                             dplyr::rename(common_name='Sample')


##for now we only have tumor and PDX data
##they each get their own sample identifier
pdx_data<-manifest|>dplyr::select(common_name,starts_with("PDX"))|>
    left_join(pdx_samps)|>
    dplyr::select(improve_sample_id,common_name,model_type,RNASeq='PDX_RNASeq',Mutations='PDX_Somatic_Mutations',CopyNumber='PDX_CNV',Proteomics='PDX_Proteomics')|>
        subset(!is.na(improve_sample_id))

tumor_data<- manifest|>dplyr::select(common_name,starts_with("Tumor"))|>
    left_join(tumor_samps)|>
    dplyr::select(improve_sample_id,common_name,model_type,RNASeq='Tumor_RNASeq',Mutations='Tumor_Somatic_Mutations',CopyNumber='Tumor_CNV')|>
    mutate(Proteomics='')|>
    subset(!is.na(improve_sample_id))
           ##we dont have tumor proteomics from these samples
#print(tumor_data)

mt_data<- manifest|>dplyr::select(common_name,starts_with("PDX"))|>
    left_join(mt_samps)|>
    dplyr::select(improve_sample_id,common_name,model_type, RNASeq='PDX_RNASeq',Mutations='PDX_Somatic_Mutations',CopyNumber='PDX_CNV',Proteomics='PDX_Proteomics')|>##we dont have mt data yet, so collecting PDX instead
    subset(!is.na(improve_sample_id))
#print(tumor_data)


combined<-rbind(pdx_data,tumor_data,mt_data)|>distinct()

# gene mapping table
genes_df <- fread(genefile)


##added proteomics first
proteomics<-do.call('rbind',lapply(setdiff(mt_data$Proteomics,c('',NA,"NA")),function(x){
                                        # if(x!=""){
    #print(x)
    sample<-subset(mt_data,Proteomics==x)
    #print(sample)
    res<-fread(synGet(x)$path)|>
        #tidyr::separate(Name,into=c('other_id','vers'),sep='\\.')|>
                                        #dplyr::select(-vers)|>
        dplyr::rename(gene_symbol='Gene')|>
        left_join(genes_df)|>
        dplyr::select(entrez_id,proteomics='logRatio')|>
        distinct()|>
        subset(!is.na(entrez_id))|>
        subset(proteomics!=0)

    res$improve_sample_id=rep(sample$improve_sample_id[1],nrow(res))
    res$source=rep('NF Data Portal',nrow(res))
    res$study=rep('MPNST PDX MT',nrow(res))
    return(distinct(res))
                                        # }
}))

fwrite(proteomics,'/tmp/mpnst_proteomics.csv.gz')


#### FIRST WE GET RNASeq Data

rnaseq<-do.call('rbind',lapply(setdiff(mt_data$RNASeq,c(NA,"NA")),function(x){
                                        # if(x!=""){
    #print(x)
    sample<-subset(mt_data,RNASeq==x)
    #print(sample)
    res<-fread(synGet(x)$path)|>
        tidyr::separate(Name,into=c('other_id','vers'),sep='\\.')|>
        dplyr::select(-vers)|>
        left_join(genes_df)|>
        dplyr::select(entrez_id,transcriptomics='TPM')|>
        subset(!is.na(entrez_id))|>
        subset(transcriptomics!=0)

    res$improve_sample_id=rep(sample$improve_sample_id[1],nrow(res))
    res$source=rep('NF Data Portal',nrow(res))
    res$study=rep('MPNST PDX MT',nrow(res))
    return(distinct(res))
                                        # }
}))

fwrite(rnaseq,'/tmp/mpnst_transcriptomics.csv.gz')



#####NEXT WE DO WES DATA
print("Getting WES")
wes<-do.call(rbind,lapply(setdiff(mt_data$`Mutations`,c(NA,"NA")),function(x){

    x2=x#gsub('"','',gsub("[",'',gsub("]",'',x,fixed=T),fixed=T),fixed=T)
    print(x)
    sample<-subset(mt_data,Mutations==x)
    print(sample$improve_sample_id)
    res<-NULL
    try(res<-fread(synGet(x2)$path)|>
            dplyr::select(entrez_id='Entrez_Gene_Id',mutation='HGVSc',variant_classification='Variant_Classification')|>
            subset(entrez_id%in%genes_df$entrez_id)|>
            distinct())
    if(is.null(res))
        return(NULL)

    res$improve_sample_id=rep(sample$improve_sample_id[1],nrow(res))
    res$source=rep('NF Data Portal',nrow(res))
    res$study=rep('MPNST PDX MT',nrow(res))

    return(distinct(res))
                                        # }
}))

fwrite(wes,'/tmp/mpnst_mutations.csv.gz')


print(paste("getting CNV"))
##next let's do CNVs!
cnv<-do.call(rbind,lapply(setdiff(mt_data$CopyNumber,c(NA,"NA")),function(x){

    x2=x#gsub('"','',gsub("[",'',gsub("]",'',x,fixed=T),fixed=T),fixed=T)
    print(x)
    sample<-subset(mt_data,CopyNumber==x)
    print(sample$improve_sample_id)
    res<-fread(synGet(x2)$path)

    long_df<- res|>
      tidyr::separate_rows(gene,sep=',')|>
      dplyr::rename(gene_symbol='gene')|>
      dplyr::left_join(genes_df)|>
      subset(!is.na(entrez_id))|>
      dplyr::select(entrez_id,log2)|>
      dplyr::distinct()|>
        dplyr::mutate(copy_number=2^log2)|>
        dplyr::select(-log2)

  res<-long_df|> ##deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
      dplyr::mutate(copy_call=ifelse(copy_number<0.5210507,'deep del',
                                     ifelse(copy_number<0.7311832,'het loss',
                                            ifelse(copy_number<1.214125,'diploid',
                                                   ifelse(copy_number<1.422233,'gain','amp')))))|>
    mutate(study='MPNST PDX MT',source='NF Data Portal',improve_sample_id=sample$improve_sample_id[1])|>
    dplyr::distinct()

    # long_df <- res[, strsplit(as.character(gene), ","), by = .(chromosome, start, end, depth, log2)]
    # filtered_df <- long_df |>
    #     subset(is.finite(log2))|>
    #     filter(V1 %in% genes_df$gene) # get only protein coding genes and remove empty gene symbols
    # filtered_df <- filtered_df[, .(gene_symbol = V1,
    #                        improve_sample_id = sample$improve_sample_id[1],
    #                        copy_number = 2^log2,
    #                        source = "NF Data Portal",
    #                        study = "MPNST PDX MT")]
    # res<-filtered_df|> ##deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
    #     dplyr::mutate(copy_call=ifelse(copy_number<0.5210507,'deep del',
    #                                    ifelse(copy_number<0.7311832,'het loss',
    #                                           ifelse(copy_number<1.214125,'diploid',
    #                                           ifelse(copy_number<1.422233,'gain','amp')))))|>
    #     left_join(genes_df)|>
    #     dplyr::select(entrez_id,improve_sample_id,copy_number,copy_call,study,source)|>
    #     subset(!is.na(entrez_id))|>
    #     distinct()
    # res|>group_by(copy_call)|>summarize(n_distinct(entrez_id))
    return(res)
                                        # }
}))

fwrite(cnv,'/tmp/mpnst_copy_number.csv.gz')

##TODO: get proteomics!!!
