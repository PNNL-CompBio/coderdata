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

##how many cancers can we cover in both diseases
dual_dat <-subset(fulldat,cancer_type%in%setdiff(cptac$cancer_type,'Glioblastoma'))
##looking for exacat matches
stats<-dual_dat|>
  group_by(cancer_type,model_type)|>summarize(numSamps=n_distinct(improve_sample_id))

fig1<-ggplot(stats,aes(x=cancer_type,y=numSamps,fill=model_type))+
  geom_bar(stat='identity',position='dodge')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle('exact tumor type matches')

ggsave('cancerDataAvailable.pdf',fig1,height=5)

##now get molecular data

gex<-readr::read_csv('cptac/transcriptomics.csv.gz')|>
 # dplyr::rename(expression='transcriptomics')|>
  rbind(readr::read_csv('cell_line/transcriptomics.csv.gz')|>
        dplyr::select(entrez_id,improve_sample_id,transcriptomics))|>
  subset(improve_sample_id%in%dual_dat$improve_sample_id)


mx<-function(x){  mean(x[is.finite(x)],na.rm=T)}

gmat<-gex|>
  subset(!is.na(entrez_id))|>
  tidyr::pivot_wider(names_from='improve_sample_id',values_from='transcriptomics',
                     values_fn=list(transcriptomics=mx),values_fill=0.0)|>
  tibble::column_to_rownames('entrez_id')

library(ggfortify)
#rmovenas
#gmat[which(!is.finite(gmat),arr.ind=T)]<-NA

rav<-which(apply(gmat,1,function(x) any(is.na(x))))
nna<-which(colnames(gmat)!='NA')
gmat2<-gmat[-rav,nna]
pres<-prcomp(t(as.matrix(gmat2)),.scale=TRUE)

md<-fulldat|>dplyr::select(improve_sample_id,model_type,cancer_type)|>
  distinct()

p<-autoplot(pres,data=subset(md,improve_sample_id%in%gex$improve_sample_id),color='cancer_type')
p2<-autoplot(pres,data=subset(md,improve_sample_id%in%gex$improve_sample_id),color='model_type')

ggsave('samples_by_cancer.pdf',p)
ggsave('samples_by_model.pdf',p2)


gsamps<-unique(gex$improve_sample_id)
fulldat$Gex=fulldat$improve_sample_id%in%gsamps

prot<-readr::read_csv('cptac/proteomics.csv.gz')|>
  rbind(readr::read_csv('cell_line/proteomics.csv.gz'))|>
  subset(improve_sample_id%in%dual_dat$improve_sample_id)

psamps<-unique(prot$improve_sample_id)
fulldat$Proteomics=fulldat$improve_sample_id%in%psamps

##now we can do pca of the proteomics
pmat<-prot|>
  subset(!is.na(entrez_id))|>
  tidyr::pivot_wider(names_from='improve_sample_id',values_from='proteomics',
                     values_fn=list(proteomics=mx),values_fill=list(proteomics=0.0))|>
  tibble::column_to_rownames('entrez_id')

rav<-which(apply(pmat,1,function(x) any(is.na(x))))
nna<-which(colnames(pmat)!='NA')
pmat2<-pmat[-rav,nna]

pres<-prcomp(t(as.matrix(pmat2)),.scale=TRUE)
p<-autoplot(pres,data=subset(md,improve_sample_id%in%prot$improve_sample_id),color='cancer_type')
p2<-autoplot(pres,data=subset(md,improve_sample_id%in%prot$improve_sample_id),color='model_type')

ggsave('prot_samples_by_cancer.pdf',p)
ggsave('prot_samples_by_model.pdf',p2)


cp<-readr::read_csv('cptac/CNV.csv.gz')
cc<-readr::read_csv('cell_line/copy_number.csv.gz')

##pca of the copy number

csamps<-union(cp$improve_sample_id,cc$improve_sample_id)
fulldat$CNV=fulldat$improve_sample_id%in%csamps

gstats<- fulldat|>
  subset(cancer_type%in%cptac$cancer_type)|>
  group_by(cancer_type,model_type,Proteomics,CNV,Gex)|>
  summarize(numSamps=n_distinct(improve_sample_id))|>
  tidyr::pivot_longer(c(Gex,Proteomics,CNV),names_to='expr',values_to='val')|>
  subset(val)|>
  dplyr::select(-val)

fig2<-ggplot(gstats,aes(x=cancer_type,fill=model_type,y=numSamps))+geom_bar(position='dodge',stat='identity')+
  facet_grid(expr~.)+
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('allDataTypes.pdf',fig2)



###now can we compare omics measurements?ß

