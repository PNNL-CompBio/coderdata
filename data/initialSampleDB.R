##this file creates a sample database table from Priyanka's depmap file

require(remotes)
if(!require(TidyComb)){
  BiocManager::install("synergyfinder")
  remotes::install_github("DrugComb/TidyComb")
}


tab<-read.table('DepMap_Argonne_Mapping.csv',sep=',',header=T)


##query for cellosaurus

##add column for cell line type

##add column for disease