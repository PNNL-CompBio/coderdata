##remap drugs
library(readr)
drugs<-readr::read_tsv("cell_line/drugs.tsv.gz")
exp<-readr::read_tsv('cell_line/experiments.tsv.gz')


drugids<-unique(drugs$improve_drug_id)
pubchem<-drugids[grep('PC',drugids)]

print(paste("Found",length(drugids),'drug ids from original dataset, of which',length(pubchem),'of those are in pubChem'))


smiles<-unique(drugs$canSMILES)
print(paste("Found",length(smiles),'unique SMILES'))

newdrugids<-data.frame(canSMILES=smiles,newid=paste("SMI",seq(1:length(smiles)),sep='_'))

##drugs with new id
newdrugs<-drugs|>left_join(newdrugids)
