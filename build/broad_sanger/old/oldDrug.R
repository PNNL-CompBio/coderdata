###Here is a script that runs through all the data files one by one.

#this is a helper file that loads the data
source("mapDrugsToPubchem.R")

#if(!require('PharmacoGx')){
#  BiocManager::install("PharmacoGx",force=TRUE)
library('PharmacoGx')
#}

all.dsets<-PharmacoGx::availablePSets()



#' getCellLineData - gets cell line dose response data
getDepMapDrugData<-function(cell.lines=c('CTRPv2','FIMM','gCSI','PRISM','GDSC','NCI60','CCLE')){


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
            drug.map<-buildDrugTable(unique(mapping$treatmentid),'/tmp/drugs.tsv.gz')%>%
                dplyr::select(common_drug_name='chem_name',improve_drug_id)%>%
                distinct()

            ##clean up file when done
            file.remove(paste0(f,'.rds'))

        }
    }
}





main<-function(){
	args = commandArgs(trailingOnly=TRUE)
	if(length(args)!=1){
	  print('Usage: Rscript 03-createDrugFile.R [datasets]')
	  exit()
	  }
#	sfile = args[1]
        dsets<-unlist(strsplit(args[1],split=','))



       dl1<-getDepMapDrugData(dsets)


}

main()
