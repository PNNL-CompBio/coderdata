##cdata summary

library(ggplot2)
library(dplyr)
library(readr)
##first we get sample files

###get all samples
cptac<-readr::read_csv('cptac/samples.csv')
cell_line<-readr::read_csv('cell_line/samples.csv')
hcmi<-c()



##now we join them into a single table, with cancer type
fulldat<-rbind(cptac,cell_line)

##how many sample by cancer type and model type?

##looking for exacat matches
stats<-fulldat|>subset(cancer_type%in%cptac$cancer_type)|>
  group_by(cancer_type,model_type)|>summarize(numSamps=n_distinct(improve_sample_id))

fig1<-ggplot(stats,aes(x=cancer_type,y=numSamps,fill=model_type))+
  geom_bar(stat='identity',position='dodge')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggitle('exact tumor type matches')

##now get molecular data

gex<-readr::read_csv('cptac/transcriptomics.csv.gz')|>
  dplyr::rename(expression='transcriptomics')|>
  rbind(readr::read_csv('cell_line/expression.csv.gz')|>
        dplyr::select(entrez_id,improve_sample_id,expression))

gsamps<-unique(gex$improve_sample_id)
fulldat$hasGex=fulldat$improve_sample_id%in%gsamps

prot<-readr::read_csv('cptac/proteomics.csv.gz')|>
  rbind(readr::read_csv('cell_line/proteomics.csv.gz'))
psamps<-unique(prot$improve_sample_id)
fulldat$hasProt=fulldat$improve_sample_id%in%psamps

cp<-readr::read_csv('cptac/CNV.csv.gz')
cc<-readr::read_csv('cell_line/copy_number.csv.gz')

csamps<-union(cp$improve_sample_id,cc$improve_sample_id)
fulldat$hasCNV=fulldat$improve_sample_id%in%csamps

gstats<- fulldat|>
  subset(cancer_type%in%cptac$cancer_type)|>
  group_by(cancer_type,model_type,hasGex,hasCNV,hasProt)|>
  summarize(numSamps=n_distinct(improve_sample_id))

fig2<-ggplot(fulldat)

