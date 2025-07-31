###Here is a script that runs through all the data files one by one.
library('reticulate')
use_python("/opt/venv/bin/python3", required = TRUE)
library('tidyr')
library(dplyr)
library(data.table)
#this is a helper file that loads the data
source_python("pubchem_retrieval.py")

library('PharmacoGx')


all.dsets<-PharmacoGx::availablePSets()



#' getCellLineData - gets cell line dose response data
getDepMapDrugData<-function(cell.lines=c('CTRPv2','FIMM','gCSI','PRISM','GDSC','CCLE'),efile=''){

    if(efile!=''){
        existing_ids=readr::read_tsv(efile)
    }else{
        existing_ids=NULL
        }
    output_file_path <- '/tmp/depmap_sanger_drugs.tsv'
    ignore_file_path <- '/tmp/ignore_chems.txt'
    for(cel in cell.lines){

        files<-subset(all.dsets,`Dataset Name`==cel)%>%
          dplyr::select(`PSet Name`)%>%
            unlist()

        for(f in files){
            print(f)
            if(f=='GDSC_2020(v2-8.2)')
                cel='GDSCv2'
            if(f=='GDSC_2020(v1-8.2)')
                cel='GDSCv1'
            #tmpfile<-paste0(cel,'doseResponse')

            dset<<-downloadPSet(f,saveDir='.',timeout=10000)


            mapping <- sensitivityInfo(dset)##get the dataset dose response data

            ##map the rownames
            if(!'exp_id'%in%names(mapping))
                mapping<-mapping%>%
                    tibble::rownames_to_column('exp_id')

            ##fix up the NSC ids
            if('NSC'%in%names(mapping))
                mapping<-mapping|>
                    dplyr::select(-treatmentid)|>
                    dplyr::mutate(treatmentid=paste0('NSC-',NSC))

            ##move drug to treatment id
            if("drugid"%in%names(mapping))
                mapping<-dplyr::rename(mapping,treatmentid='drugid')

            ##query to build the drug ids
            chem_list <- unique(mapping$treatmentid)
            if (!is.null(existing_ids)) {
                chem_list <- setdiff(chem_list, existing_ids$chem_name)
            }
            print(paste('Found',length(chem_list),'chemicals for dataset',cel))


            #build prev_drug_filepaths for this iteration (include original and accumulated) ---
            prev_list <- c()
            if (!is.na(efile) && nzchar(efile)) prev_list <- c(prev_list, efile)
            if (file.exists(output_file_path)) prev_list <- c(prev_list, output_file_path)
            prev_drug_filepaths <- if (length(prev_list) > 0) paste(prev_list, collapse = ",") else NULL

            # write this PSet's results to a temp file
            temp_out <- paste0(output_file_path, ".part.tsv")
            update_dataframe_and_write_tsv(
                unique_names           = chem_list,
                output_filename        = temp_out,
                ignore_chems           = ignore_file_path,
                batch_size             = 50,
                isname                 = TRUE,
                prev_drug_filepaths    = prev_drug_filepaths,
                restrict_to_raw_names  = chem_list
            )

            # --- merge temp_out into cumulative output_file_path ---
            if (file.exists(temp_out)) {
                if (file.exists(output_file_path)) {
                    agg <- fread(output_file_path, sep="\t", header=TRUE)
                    new_part <- fread(temp_out, sep="\t", header=TRUE)
                    combined <- unique(rbindlist(list(agg, new_part), use.names=TRUE, fill=TRUE))
                } else {
                    combined <- fread(temp_out, sep="\t", header=TRUE)
                }
                fwrite(combined, output_file_path, sep="\t")
                file.remove(temp_out)
            } else {
                warning("Expected temporary output not found: ", temp_out)
            }

            ##clean up file when done
            file.remove(paste0(f,'.rds'))

        }
    }
}



main<-function(){
	args = commandArgs(trailingOnly=TRUE)
	if(length(args)<2){
	  print('Usage: Rscript 03-createDrugFile.R [datasets] [existing file]')
# 	  exit()
	  }
#	sfile = args[1]
        dsets<-unlist(strsplit(args[1],split=','))
        if(length(args)==2)
            efile=args[2]
        else
            efile=''
       dl1<-getDepMapDrugData(dsets,efile)


}

main()
