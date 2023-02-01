#takes the pubchem list and mapps all the drug

#the challenge is that there are millions of chemicals out there, we need to restrict somehow
##

require(remotes)
if(!require(TidyComb)){
  BiocManager::install("synergyfinder")
  remotes::install_github("DrugComb/TidyComb")
}

library(curl)
library(dplyr)

ids<-c()

TidyComb::GetPubNames(ids)
TidyComb::GetPubchemPro(ids)