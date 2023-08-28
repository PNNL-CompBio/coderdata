##cdata summary

library(ggplot2)
library(dplyr)
library(readr)
##first we get sample files

###get all samples and their metadata
cptac<-readr::read_csv('cptac/samples.csv')
cell_line<-readr::read_csv('cell_line/samples.csv')
hcmi<-c()

##now we join them into a single table, with cancer type
fulldat<-rbind(cptac,cell_line)

##how many cancers can we cover in tumor types
dual_dat <-subset(fulldat,cancer_type%in%setdiff(cptac$cancer_type,'Glioblastoma'))


##looking for exact matches
stats<-dual_dat|>
  group_by(cancer_type,model_type)|>
  summarize(numSamps=n_distinct(improve_sample_id))

fig1<-ggplot(stats,aes(x=cancer_type,y=numSamps,fill=model_type))+
  geom_bar(stat='identity',position='dodge')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle('exact tumor type matches')
print(fig1)
ggsave('fig0_cancerDataAvailable.pdf',fig1,height=5)

##now get molecular data
##takes a data frame from the data and plots the reduced dimensions
doPlotFromDF<-function(gex,value){
  mx<-function(x){  mean(x[is.finite(x)],na.rm=T)}

  ##this is the worst part! why does it take so long???
  print('converting to matrix...')
  gmat<-gex|>
    subset(!is.na(entrez_id))|>
    dplyr::rename(eval=value)|>
    dplyr::select(entrez_id,improve_sample_id,eval)|>
    tidyr::pivot_wider(names_from='improve_sample_id',values_from='eval',
                       values_fn=list(eval=mx),values_fill=0.0)|>
    tibble::column_to_rownames('entrez_id')

  library(umap)
  library(Rtsne)
  library(ggplot2)
#library(ggfortify)
#rmovenas
#gmat[which(!is.finite(gmat),arr.ind=T)]<-NA

  md<-fulldat|>dplyr::select(improve_sample_id,model_type,cancer_type)|>
    subset(improve_sample_id%in%dual_dat$improve_sample_id)|>
    distinct()

  rav<-which(apply(gmat,1,function(x) any(is.na(x))))
  nav<-which(apply(gmat,1,function(x) any(is.nan(x))))
  missing <- union(rav,nav)
  if(length(missing)>0)
    gmat1<-gmat[-missing,]
  else
    gmat1<-gmat
  nna<-which(colnames(gmat1)!='NA')
  gmat2<-gmat1[,nna]

  ##zero variance values could be a challenge
  zvar<-which(apply(gmat2,1,var)==0)
 # pres<-prcomp(t(as.matrix(gmat2)),.scale=TRUE)

  print('calculating umap')
  ures<-umap(t(as.matrix(gmat2)))
  udf <- data.frame(x = ures$layout[,1],
                   y = ures$layout[,2],
                 subset(md,improve_sample_id%in%rownames(ures$layout)))

  p1<-ggplot(udf,aes(x=x,y=y,col=cancer_type))+geom_point()
  p2<-ggplot(udf,aes(x=x,y=y,col=model_type))+geom_point()
  ggsave(paste0(value,'umapCancerTypeTypeDimRed.pdf'),p1,width=14)
  ggsave(paste0(value,'umapCellTypeDimRed.pdf'),p2,width=10)


  print('calculating tsne')
  tsne<-Rtsne(t(as.matrix(gmat2)))
  tdf <- data.frame(x = tsne$Y[,1],
                    y = tsne$Y[,2],
                    subset(md,improve_sample_id%in%rownames(ures$layout)))

  p3<-ggplot(tdf,aes(x=x,y=y,col=cancer_type))+geom_point()
  ggsave(paste0(value,'tsneCancerTypeTypeDimRed.pdf'),p3,width=14)

  p4<-ggplot(tdf,aes(x=x,y=y,col=model_type))+geom_point()
  ggsave(paste0(value,'tsneCellTypeDimRed.pdf'),p4,width=10)


 # fp<-cowplot::plot_grid(p2,p4,nrow=1)

  #  ggsave(paste0(value,'CellTypeDimRed.pdf'),fp,width=10)

  #fp<-cowplot::plot_grid(p1,p3,nrow=1)
  #ggsave(paste0(value,'CancerTypeTypeDimRed.pdf'),fp,width=14)
}

#PLOT gene expression, then proteomics, then mirnas, then copy number
#gsamps<-unique(gex$improve_sample_id)
#fulldat$Gex=fulldat$improve_sample_id%in%gsamps
gex<-readr::read_csv('cptac/transcriptomics.csv.gz')|>
  # dplyr::rename(expression='transcriptomics')|>
  rbind(readr::read_csv('cell_line/transcriptomics.csv.gz')|>
          dplyr::select(entrez_id,improve_sample_id,transcriptomics,source,study))|>
  subset(improve_sample_id%in%dual_dat$improve_sample_id)

doPlotFromDF(gex,'transcriptomics')

prot<-readr::read_csv('cptac/proteomics.csv.gz')|>
  dplyr::select(improve_sample_id,entrez_id,proteomics)|>
  rbind(readr::read_csv('cell_line/proteomics.csv.gz'))|>
  subset(improve_sample_id%in%dual_dat$improve_sample_id)

doPlotFromDF(prot,'proteomics')

cp<-readr::read_csv('cptac/CNV.csv.gz')|>
  dplyr::rename(copy_number='CNV')
cc<-readr::read_csv('cell_line/copy_number.csv.gz')
cnv<-rbind(cp,cc)

##first let's look at the CNVs a bit
cdat<-cnv|>
  subset(!is.na(copy_call))|>
  left_join(md)

cdat<-cdat|>tidyr::replace_na(list(model_type='cell line'))|>
  tidyr::replace_na(list(cancer_type='other'))

p1<-ggplot(cdat,aes(x=cancer_type,fill=copy_call))+
  geom_bar(position='dodge')+
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle('Copy number spread')+
  facet_grid(~model_type)

ggsave('copyNumberSpread.pdf',p1)

print('removing diploid measurements')
nd<-subset(cnv,copy_call!='diploid')
doPlotFromDF(nd,'copy_number')

##TODO: explore somatic variants

fullStats<-function(gex,prot,cnv,mut){
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

}


###now can we compare omics measurements?ß

