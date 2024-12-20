##remap drugs
library(readr)
library(dplyr)

main<-function(drugfile='/tmp/drugs.tsv.gz',expfile='/tmp/experiments.tsv.gz'){

  drugs<-readr::read_tsv(drugfile)
  exp<-readr::read_tsv(expfile)

  ##first copy files to same path with 'orig'
  drugids<-unique(drugs$improve_drug_id)
  pubchem<-drugids[grep('PC',drugids)]

  print(paste("Found",length(drugids),'drug ids from original dataset, of which',length(pubchem),'of those are in pubChem'))

  smiles<-unique(drugs$canSMILES)
  print(paste("Found",length(smiles),'unique SMILES'))

  newdrugids<-data.frame(canSMILES=smiles,newid=paste("SMI",seq(1:length(smiles)),sep='_'))|>
    subset(!is.na(canSMILES))

  ##drugs with new id
  newdrugs<-drugs|>dplyr::right_join(newdrugids)|>
    dplyr::rename(pubchem_id='improve_drug_id',improve_drug_id='newid')|>
    dplyr::distinct()

  newexp<-exp|>
    dplyr::rename(pubchem_id='Drug')|>
    dplyr::right_join(dplyr::select(newdrugs,c(improve_drug_id,pubchem_id)))|>
    dplyr::select(-pubchem_id)|>
      dplyr::distinct()

  colnames(newexp)<-tolower(colnames(newexp))

  readr::write_tsv(newdrugs,drugfile)
  readr::write_csv(newexp,expfile)
}

args<-commandArgs()
print(args)
if(length(args)==7){
  df=args[6]
  ef=args[7]
  main(df,ef)
}else{
  main()
}

