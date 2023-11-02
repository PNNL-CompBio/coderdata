##cdata summary

library(ggplot2)
library(dplyr)
library(readr)
##first we get sample files

mergeSamples<-function(){

  ##########################
  ## CANCER NAME MAPPING FILE
  ##this file tries to map as many cancer types as possible to groups that span model systems
  ###########################
  cmaps<-readr::read_csv('../cellLineTypes.csv')

  ###########################
  ## CPTAC SAMPLE DATA
  ## Primary task here is to fix the capitalizations
  ###########################
  cptac<-readr::read_csv('../build/cptac/samples.csv')|>
    mutate(cancer_type=stringr::str_replace_all(cancer_type,'Head and Neck','Head and neck'))|>
    mutate(cancer_type=stringr::str_replace_all(cancer_type,'Colon','Colorectal'))|>
    mutate(cancer_type=stringr::str_replace_all(cancer_type,'Uterine Corpus Endometrial Carcinoma','Uterine corpus endometrial carcinoma'))|>
    dplyr::mutate(`CPTAC Cancer type`=cancer_type)|>
    left_join(cmaps)|>
    mutate(sampleSource='CPTAC')|>
    dplyr::select(improve_sample_id,`CPTAC Cancer type`,model_type,species,sampleSource)|>
    distinct()

  ###########################
  ## CCLE Cell Line data
  ## We have many more cancer types here, so we try to map what we have to CPTAC names and then adjust the rest
  ##
  ###########################
  cell_line<-readr::read_csv('../build/cell_line/samples.csv')|>
    dplyr::mutate(`Cell line cancer type`=cancer_type)|>
    mutate(sampleSource='CCLE')

  allec<-grep('Endometrial',cell_line$`Cell line cancer type`)
  cell_line$`Cell line cancer type`[allec]<-'Uterine corpus endometrial carcinoma'

  cell_line<-cell_line|>
    left_join(cmaps)

  ##first we collect the names of the cancers that are NOT in CPTAC
  other_cans<-which(is.na(cell_line$`CPTAC Cancer type`))
  cell_line$`CPTAC Cancer type`[other_cans]<-cell_line$`Cell line cancer type`[other_cans]
  cell_line<-cell_line|>
    dplyr::select(improve_sample_id,`CPTAC Cancer type`,model_type,species,sampleSource)|>
     distinct()

  #then we rename the NA values to 'Other' if we want
  #other_cans<-which(is.na(cell_line$`CPTAC Cancer type`))
  #cell_line$`CPTAC Cancer type`[other_cans]<-'Other'
  # or just remove them
  cell_line<-cell_line|>
    subset(!is.na(`CPTAC Cancer type`))

  ###########################
  ## HCMI SAMPLE DATA
  ## Here we lean heavily on the sample mapping file
  ###########################
  hcmi<-readr::read_csv('../build/hcmi/samples.csv')|>
    dplyr::rename(id_source='other_id_source')|>
    mutate(species='human')|>
    subset(model_type%in%c('3D Organoid','Solid Tissue','Adherent Cell Line'))|>
    dplyr::mutate(`HCMI Cancer type`=cancer_type,`HCMI Common name`=common_name)|>
    dplyr::mutate(model_type=stringr::str_replace_all(model_type,'Solid Tissue','Tumor'))|>
    dplyr::mutate(model_type=stringr::str_replace_all(model_type,'Adherent Cell Line','cell line'))|>
    left_join(cmaps)|>
    mutate(sampleSource='HCMI')|>
    dplyr::select(improve_sample_id,`CPTAC Cancer type`,model_type,species,sampleSource)|>
#    dplyr::rename(cancer_type='CPTAC Cancer type')|>
    distinct()

  ##now rename samples
  ##next up: beatAMLdata
  ###########################
  ## BeatAML SAMPLE DATA
  ## TBD
  ###########################
  baml<-readr::read_csv("../build/beatAML/samples.csv")|>
    mutate(cancer_type='Acute myeloid leukemia')|>
    mutate(species='Human')|>
    mutate(model_type='Tumor')|>
    mutate(sampleSource='BeatAML')|>
    dplyr::select(improve_sample_id,species,cancer_type,sampleSource,model_type)|>
    distinct()
  ###########################
  ## TCGA SAMPLE DATA
  ## TBD
  ###########################


  ##now we join thomdelsem into a single table, with cancer type
  fulldat<<-rbind(cptac,cell_line,hcmi)|>
    dplyr::rename(cancer_type=`CPTAC Cancer type`)|>
    subset()
  #  subset(!is.na(cancer_type))


#  fulldat<-fulldat|>
#    dplyr::rename(orig_cancer_type='cancer_type')|>
#    dplyr::rename(cancer_type='CPTAC Cancer type')

  fulldat<-fulldat|>
    dplyr::select(improve_sample_id,species,cancer_type,sampleSource,model_type)|>
    distinct()|>
    rbind(baml)

  models<-fulldat|>
    group_by(cancer_type)|>
    summarize(num_models=n_distinct(model_type))|>
    subset(num_models>1)

  alldat<-fulldat|>
    subset(cancer_type%in%models$cancer_type)

  other_can <- fulldat|>
    subset(!is.na(cancer_type))|>
    subset(!cancer_type%in%models$cancer_type)|>
    mutate(cancer_type='Other')

#    subset(!cancer_type%in%models
#  fulldat[which(fulldat$cancer_type%in%other_can),'cancer_type']<-'Other cancer'

  return(rbind(alldat,other_can))
}
# ###get all samples and their metadata
# cptac<-readr::read_csv('../build/cptac/samples.csv')|>
#   mutate(cancer_type=stringr::str_replace_all(cancer_type,'Head and Neck','Head and neck'))
#
# #cell_line<-readr::read_csv('../build/cell_line/samples_curated.csv')
# cell_line<-readr::read_csv("../build/cell_line/samples.csv")
# allec<-grep('Endometrial',cell_line$cancer_type)
# cell_line$cancer_type[allec]<-'Uterine Corpus Endometrial Carcinoma'
#
#
# hcmi<-c()
#
# ##now we join them into a single table, with cancer type
# fulldat<<-rbind(cptac,cell_line)|>
#   subset(!is.na(cancer_type))
#
# other_can <- setdiff(fulldat$cancer_type,cptac$cancer_type)
# fulldat[which(fulldat$cancer_type%in%other_can),'cancer_type']<-'Other cancer'
#

fulldat<-mergeSamples()

##looking for exact matches
stats<-fulldat|>
  subset(!is.na(cancer_type))|>
  group_by(cancer_type,model_type,sampleSource)|>
  summarize(numSamps=n_distinct(improve_sample_id))|>
  subset(model_type!='Not Reported')|>
  subset(numSamps>1)


fig1<-ggplot(stats,aes(x=cancer_type,y=numSamps,fill=model_type))+
  geom_bar(stat='identity',position='dodge')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_log10()+scale_fill_viridis_d()+
  ggtitle('Samples by tumor type')

print(fig1)
ggsave('fig0_cancerDataAvailable.pdf',fig1,height=8,width=10)

fig2<-ggplot(stats,aes(x=cancer_type,y=numSamps,fill=model_type))+
  geom_bar(stat='identity',position='dodge')+
  facet_grid(~sampleSource)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle('exact tumor type matches')

runCombat<-function(mat, fulldat){
  library(sva)
  #mat<-readr::read_csv(fname)
  batch<-fulldat|>
    subset(improve_sample_id%in%colnames(mat))|>
    select(improve_sample_id,model_type)|>distinct()
  bl<-batch$model_type
  names(bl)<-batch$improve_sample_id

  ttype<-fulldat|>select(improve_sample_id,cancer_type)|>
    subset(improve_sample_id%in%colnames(mat))|>
    distinct()
  tlist<-as.factor(ttype$cancer_type)
  names(tlist)<-ttype$improve_sample_id
  tlist<-tlist[colnames(mat)]
  bl<-bl[colnames(mat)]
  res<-ComBat(mat,batch=bl)#,mod=tlist)

  res
}


##standard making a matrix from long table
matFromDF<-function(gex,value){
  mx<-function(x){  mean(x[is.finite(x)],na.rm=T)}

  if(value%in%c('mutation','copy_number'))
    fill=0.0
  else
    fill=NA
  ##this is the worst part! why does it take so long???
  print('converting to matrix...')
  gmat<-gex|>
    dplyr::select('entrez_id','improve_sample_id',value)|>
    distinct()|>
    subset(!is.na(entrez_id))|>
    dplyr::rename(eval=value)|>
   # dplyr::select(entrez_id,improve_sample_id,eval)|>
    tidyr::pivot_wider(names_from='improve_sample_id',values_from='eval',
                       values_fn=list(eval=mx),values_fill=fill)|>
    tibble::column_to_rownames('entrez_id')


  md<-fulldat|>
    dplyr::select(improve_sample_id,model_type,cancer_type,sampleSource)|>
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

  pats = intersect(colnames(gmat2),md$improve_sample_id)
  gmat2<-gmat2[,pats]

  ##now let's save the matrix file
  #readr::write_csv(tibble::rownames_to_column(as.data.frame(gmat2),'gene'),file=paste0(value,'_inMatrixForm.csv'))

  return(gmat2)

}

##now get molecular data
##takes a data frame from the data and plots the reduced dimensions
doPlotFromMat<-function(gmat2,value,suff=''){
  library(umap)
  library(Rtsne)
  library(ggplot2)


  md<-fulldat|>
    dplyr::select(improve_sample_id,model_type,cancer_type,sampleSource)|>
    distinct()|>
    subset(improve_sample_id%in%colnames(gmat2))

  print('calculating umap')

  ures<-umap(t(as.matrix(gmat2)))

  over<-intersect(md$improve_sample_id,colnames(gmat2))
  udf <- data.frame(x = ures$layout[over,1],
                   y = ures$layout[over,2],
                 subset(md,improve_sample_id%in%over)|>dplyr::distinct())

  p1<-ggplot(udf,aes(x=x,y=y,col=cancer_type,shape=sampleSource))+geom_point()+
    ggtitle(paste('UMAP of',value,'by cancer type'))
  p2<-ggplot(udf,aes(x=x,y=y,col=model_type,shape=sampleSource))+geom_point()+
    ggtitle(paste('UMAP of',value,'by model type'))

  ggsave(paste0(value,suff,'_umap_CancerTypeTypeDimRed.pdf'),p1,width=14)
  ggsave(paste0(value,suff,'_umap_CellTypeDimRed.pdf'),p2,width=10)

  rudf<-subset(udf,cancer_type!='Other')
  p1<-ggplot(rudf,aes(x=x,y=y,col=cancer_type,shape=sampleSource))+geom_point()+
    ggtitle(paste('UMAP of',value,'by cancer type'))

  p2<-ggplot(rudf,aes(x=x,y=y,col=model_type,shape=sampleSource))+geom_point()+
    ggtitle(paste('UMAP of',value,'by model type'))

  ggsave(paste0(value,suff,'_umap_CancerTypeTypeDimRed_no_other.pdf'),p1,width=14)
  ggsave(paste0(value,suff,'_umap_CellTypeDimRed_no_other.pdf'),p2,width=10)


  print('calculating tsne')
  tsne<-Rtsne(t(as.matrix(gmat2)),check_duplicates=FALSE)
  tdf <- data.frame(x = tsne$Y[,1],
                    y = tsne$Y[,2],
                    subset(md,improve_sample_id%in%rownames(ures$layout)))

  p3<-ggplot(tdf,aes(x=x,y=y,col=cancer_type))+geom_point()+
    ggtitle(paste('TSNE of',value,'by cancer type'))

  ggsave(paste0(value,suff,'_TSNE_CancerTypeTypeDimRed.pdf'),p3,width=14)

  p4<-ggplot(tdf,aes(x=x,y=y,col=model_type,shape=sampleSource))+geom_point()+
    ggtitle(paste('TSNE of',value,'by model type'))
  ggsave(paste0(value,suff,'_TSNE_CellTypeDimRed.pdf'),p4,width=10)

  rtdf<-subset(tdf,cancer_type!='Other')
  p1<-ggplot(rtdf,aes(x=x,y=y,col=cancer_type,shape=sampleSource))+geom_point()+
    ggtitle(paste('TSNE of',value,'by cancer type'))
  p2<-ggplot(rtdf,aes(x=x,y=y,col=model_type,shape=sampleSource))+geom_point()+
    ggtitle(paste('TSNE of',value,'by model type'))
  ggsave(paste0(value,suff,'_TSNE_CancerTypeTypeDimRed_no_other.pdf'),p1,width=14)
  ggsave(paste0(value,suff,'_TSNE_CellTypeDimRed_no_other.pdf'),p2,width=10)

}

plotTranscripts<-function(){
  #PLOT gene expression, then proteomics, then mirnas, then copy number
  #gsamps<-unique(gex$improve_sample_id)
  #fulldat$Gex=fulldat$improve_sample_id%in%gsamps
  gex<-readr::read_csv('../build/cptac/transcriptomics.csv.gz')|>
    # dplyr::rename(expression='transcriptomics')|>
    rbind(readr::read_csv('../build/cell_line/transcriptomics.csv.gz')|>
            dplyr::select(entrez_id,improve_sample_id,transcriptomics,source,study))|>
    rbind(readr::read_csv('../build/hcmi/transcriptomics.csv'))|>
    rbind(readr::read_csv("../build/beatAML/transcriptomics.csv"))|>
    subset(improve_sample_id%in%fulldat$improve_sample_id)

  ##filter for genes expressed in all samples
  nsamples<-length(unique(gex$improve_sample_id))
  gcounts<-gex|>group_by(entrez_id)|>
    summarize(gcounts=n_distinct(improve_sample_id))|>
    subset(gcounts==nsamples)
  gex<-subset(gex,entrez_id%in%gcounts$entrez_id)
  print(paste('evaluating',nrow(gcounts),'genes and ',length(unique(gex$improve_sample_id)),'samples for expression'))
  mat<-matFromDF(gex,'transcriptomics')
  doPlotFromMat(mat,'transcriptomics')
  newmat<-runCombat(mat,fulldat)
  doPlotFromMat(newmat,'transcriptomics','_COMBAT_')

  gex
}


plotProteomics<-function(){
  prot<-readr::read_csv('../build/cptac/proteomics.csv.gz')|>
    rbind(readr::read_csv('../build/beatAML/proteomics.csv'))|>
    dplyr::select(improve_sample_id,entrez_id,proteomics)|>
    rbind(readr::read_csv('../build/cell_line/proteomics.csv.gz'))|>
    subset(improve_sample_id%in%fulldat$improve_sample_id)



  ##filter for genes expressed in all samples
  nsamples<-length(unique(prot$improve_sample_id))
  pcounts<-prot|>group_by(entrez_id)|>
    summarize(pcounts=n_distinct(improve_sample_id))|>
    subset(pcounts==nsamples)
  pex<-subset(prot,entrez_id%in%pcounts$entrez_id)
  print(paste('evaluating',nrow(pcounts),'proteins for expression'))
  mat<-matFromDF(pex,'proteomics')
  doPlotFromMat(mat,'proteomics')
  newmat<-runCombat(mat,fulldat)
  doPlotFromMat(newmat,'proteomics','_combat')

  prot
}

plotCopyNumber<-function(){
  cp<-readr::read_csv('../build/cptac/CNV.csv.gz')|>
    dplyr::rename(copy_number='CNV')
  cc<-readr::read_csv('../build/cell_line/copy_number.csv.gz')
  #ca<-readr::read_csv("../beatAML/")
  ch<-readr::read_csv('../build/hcmi/copy_number.csv')
  cnv<-rbind(cp,cc,ch)

  md<-fulldat|>
    dplyr::select(improve_sample_id,model_type,cancer_type,sampleSource)|>
    #subset(improve_sample_id%in%fulldat$improve_sample_id)|>
    distinct()

  ##first let's look at the CNVs a bit
  cdat<-cnv|>
    subset(!is.na(copy_call))|>
    left_join(md)|>
    subset(!copy_call%in%c('diploid','normal'))

 # cdat<-cdat|>tidyr::replace_na(list(model_type='cell line'))|>
#    tidyr::replace_na(list(cancer_type='other'))

  p1<-ggplot(cdat,aes(x=cancer_type,fill=copy_call))+
    geom_bar(position='dodge')+
    scale_y_log10()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle('Copy number spread')+
    facet_grid(~model_type)

  ggsave('copyNumberSpread.pdf',p1)

  print('removing diploid measurements')
  nd<-subset(cnv,copy_call!='diploid')

  mat<-matFromDF(nd,'copy_number')
  doPlotFromMat(mat,'copy_number')
  newmat<-runCombat(mat,fulldat)
  doPlotFromMat(newmat,'copy_number','_combat')

    nd
}
##TODO: explore somatic variants
plotMutations<-function(){
  library(ggplot2)
  mut<-readr::read_csv('../build/cptac/somatic_mutation.csv.gz')
   # dplyr::rename(entrez_id='entrez_gene')
  mut2<-readr::read_csv('../build/cell_line/mutations.csv.gz')|>
    dplyr::rename(mutation='mutations')

  mut3<-readr::read_csv('../build/hcmi/mutations.csv')|>
    dplyr::rename(mutation='mutations')|>
    dplyr::rename(variant_classification='variant_class')

  mut4<-readr::read_csv("../build/beatAML/mutations.csv") ##need variant classification

  allmut<-rbind(mut,mut2,mut3)|>
    left_join(fulldat)

  mp<-ggplot(subset(allmut,!is.na(model_type)),aes(x=variant_classification,fill=cancer_type))+geom_bar()+
    scale_y_log10()+
    facet_grid(~model_type)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave('mutationDataSpread.pdf',mp)

  newmut<-allmut
  newmut$mutation<-rep(1.0,nrow(newmut))
  ##now we can compute distances based on each type of mutation
  varclass<-unique(newmut$variant_classification)
  res<-lapply(varclass,function(vart){
    df<-subset(newmut,variant_classification==vart)
    mat<-matFromDF(df,'mutation')
    doPlotFromMat(mat,'mutation',vart)
  })

  newmut

}


##now collect full stats
fullStats<-function(gex,prot,cnv,mut){


  gstats<- fulldat|>
    mutate(Proteomics=improve_sample_id%in%prot$improve_sample_id)|>
    mutate(Transcriptomics=improve_sample_id%in%gex$improve_sample_id)|>
    mutate(CNV=improve_sample_id%in%cnv$improve_sample_id)|>
    mutate(Mutation=improve_sample_id%in%mut$improve_sample_id)

  summ <- gstats|>
    group_by(cancer_type,model_type,Proteomics,Transcriptomics,CNV,Mutation)|>
    tidyr::pivot_longer(c(Proteomics,Transcriptomics,CNV,Mutation),names_to='dataType',values_to='val')|>
    subset(val)|>
    group_by(cancer_type,model_type,dataType)|>summarize(numSamples=n_distinct(improve_sample_id))

  fig2<-ggplot(summ,aes(x=cancer_type,fill=model_type,y=numSamples))+geom_bar(position='dodge',stat='identity')+
    facet_grid(dataType~.)+
    scale_y_log10()+
    scale_fill_viridis_d()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave('allDataTypes.pdf',fig2)

}

main<-function(){
  prot<-plotProteomics()

  gex<-plotTranscripts()
  cnv<-plotCopyNumber()
  ##need to figure out mutations!!!

  fullStats(gex,prot,cnv,NA)
}





