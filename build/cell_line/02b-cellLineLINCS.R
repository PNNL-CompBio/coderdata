
library(readr)
library(tidyr)
library(dplyr)
library(R.utils)
library(cmapR)
#library(webchem)

#### Step 1: Get IMPROVE gene, drug, and sample IDs for mapping ####
# IMPROVE on FigShare: https://figshare.com/articles/dataset/IMPROVE_Cell_Line_Data_Files/22822286

allgenes = readr::read_csv("https://figshare.com/ndownloader/files/40576109")
genes = allgenes|>
  dplyr::select(gene_symbol,entrez_id)|>
  dplyr::distinct()

download.file("https://figshare.com/ndownloader/files/42357210",
              "drugs_by_structure.tsv.gz")
alldrugs = readr::read_tsv("drugs_by_structure.tsv.gz")
drugs = alldrugs|>
  dplyr::select(chem_name,pubchem_id,improve_drug_id)|>
  dplyr::distinct()

samples = readr::read_csv('https://figshare.com/ndownloader/files/40576103',
                   quote='"')|>
  dplyr::select(other_id,improve_sample_id)|>
  unique()

Sys.setenv(VROOM_CONNECTION_SIZE=100000000)
#### Step 2: get data from source ####
getL1000 <- function() {
  options(timeout = 300)
  
  # identify URLs
  #basename="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101406"
  L1000 <- 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE101nnn/GSE101406/suppl/GSE101406%5FBroad%5FLINCS%5FL1000%5FLevel4%5FZSPCINF%5Fmlr12k%5Fn1667x12328.gctx.gz'
  L1000.genes <- 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE101nnn/GSE101406/suppl/GSE101406%5FBroad%5FLINCS%5FL1000%5Fgene%5Finfo.txt.gz'
  L1000.inst <- 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE101nnn/GSE101406/suppl/GSE101406%5FBroad%5FLINCS%5FL1000%5Finst%5Finfo.txt.gz'
  
  # download source sample & perturbation info
  L1000.genes.info <- readr::read_delim(L1000.genes, "\t") # pr_gene_symbol, pr_gene_id
  L1000.inst.info <- readr::read_delim(L1000.inst, "\t") # pert_iname; pert_type = "trt_cp" if drug, "ctl_vehicle" if DMSO; cell_id
  L1000.inst.info <- L1000.inst.info[L1000.inst.info$pert_type == "trt_cp", ] # only keep drug treatments (remove DMSO controls)
  
  # download L1000 data
  res<-download.file(L1000,'L1000.gctx.gz', mode="wb")
  R.utils::gunzip('L1000.gctx.gz', 'L1000.gctx') # unzip file before reading it
  L1000.df<-cmapR::parse_gctx('L1000.gctx') 

  # put data into long format
  L1000.long <- cmapR::melt_gct(L1000.df)
  colnames(L1000.long) <- c("inst_id", "pr_gene_id", "data_value")
  L1000.long$pr_gene_id <- as.numeric(L1000.long$pr_gene_id)

  # join with L1000.gene.info based on pr_gene_id;
  # L1000.inst.info based on inst_id;
  # keep only cell_id, pert_type, pert_iname, gene_symbol, data_value
  L1000.full <- L1000.long |> 
    dplyr::left_join(L1000.gene.info) |>
    dplyr::left_join(L1000.inst.info) |>
    dplyr::select(cell_id,pert_type,pert_iname,pr_gene_symbol,data_value)|>
    dplyr::distinct()
  L1000.long <- NULL # save space

  # add relevant columns
  L1000.full <- na.omit(L1000.full)
  colnames(L1000.full)[1] <- "other_id" # match samples column name
  colnames(L1000.full)[4] <- "gene_symbol" # match genes column name
  colnames(L1000.full)[3] <- "chem_name"
  L1000.full$data_type <- "transcriptomics"
  L1000.full$source <- "CMap"
  L1000.full$study <- "LINCS"
  L1000.full$perturbation_type <- "drug" # all entries have pert_type="trt_cp"
  
  # join with IMPROVE IDs:
  # samples by "other_id"
  # genes by "gene_symbol"
  # drugs by "chem_name"
  res<-L1000.full|>
    dplyr::left_join(samples)|>
    dplyr::left_join(genes)|>
    dplyr::left_join(drugs)|>
    dplyr::select(entrez_id,improve_sample_id,data_value,data_type,
                  improve_drug_id,perturbation_type,source,study)|>
    dplyr::distinct()
  colnames(res)["improve_drug_id"] <- "perturbation"
  write_csv(res,file=gzfile('perturbations.csv.gz'))
}

getCMap <- function() {
  getL1000()
  #getCRISPR()
  #getScreen()
  #getP100()
  #getGCP()
}

#### Step 3: rewrite each file ####
filenames=list(perturbations='insert-figshare-link')
newres<-lapply(names(filenames),function(value){

  fi=filenames[[value]]
  fname=paste0(value,'.csv.gz')
  print(paste('now reading',fi,'to store as',fname))
  ## now every data type is parsed slightly differently,
  ## so we need to change our formatting and mapping
  ## to get it into a unified 3 column schema

    if(value=='perturbations'){ # if perturbations in transcriptomics
      exp_file <- readr::read_csv(fi)

      res = tidyr::pivot_longer(data=exp_file,cols=c(2:ncol(exp_file)),
                                names_to='gene_entrez',values_to='data_value',
                                values_transform=list(expression=as.numeric))|>
        tidyr::separate_wider_delim(gene_entrez,' ',names=c('gene','entrez_par'))|>
        mutate(entrez_id=stringr::str_replace_all(entrez_par,'\\)|\\(',''))|>
        dplyr::select(-c(entrez_par,gene))|>
        distinct()

      colnames(res)[1]<-'other_id'
      vars=c('data_value','data_type','perturbation','perturbation_type')
    }
    
  ##do the last join with samples
  full<-res|>
    dplyr::left_join(samples)|>
    dplyr::select(c('entrez_id','improve_sample_id',vars))|>
    dplyr::distinct()|>
    dplyr::mutate(source='Broad',study='LINCS')

  write_csv(full,file=gzfile(fname))
  return(fi)

})

getCMap()
