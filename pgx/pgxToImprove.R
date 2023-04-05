###Here is a script that runs through all the data files one by one. 

#this is a helper file that loads the data
source("improveDataFiles.R")


#if(!require('PharmacoGx')){
#  BiocManager::install("PharmacoGx",force=TRUE)
  library('PharmacoGx')
#}

all.dsets<-PharmacoGx::availablePSets()
##first define a generic dose response function

#' getDoseRespData
#' Generic function to get dose and response data from PGX object
#' out of dataset object, and store with dataset name
getDoseRespData<-function(dset,studyName){
  
  mapping <- sensitivityInfo(dset)
  if(!'exp_id'%in%names(mapping))
    mapping<-mapping%>%
      tibble::rownames_to_column('exp_id')
  
  if("drugid"%in%names(mapping))
    mapping<-dplyr::rename(mapping,treatmentid='drugid')
  if('cellid'%in%names(mapping))
    mapping<-dplyr::rename(mapping,sampleid='cellid')
  
  ##now get the drug ids
  drug.map<-buildDrugTable(unique(mapping$treatmentid))%>%
    dplyr::select(common_drug_name='chem_name',improve_chem_id)%>%
    distinct()
  
  ##first get the sample id
  samp.map<-mapping%>%
    dplyr::select(sampleid,exp_id,treatmentid)%>%distinct()
  
  comm.map<-samp.map%>%
    dplyr::rename(other_names='sampleid')%>%
    left_join(improve_samples)%>%
    subset(!is.na(improve_sample_id))
  
  print(paste('By common name, found',length(unique(comm.map$improve_sample_id)),
              'matches out of',length(unique(mapping$sampleid)),'for study',studyName))
  
  ldmap<-drug.map%>%
    mutate(chem_name=tolower(common_drug_name))
  ##then join the sample id
  full.map<-comm.map%>%
    dplyr::select(exp_id,improve_sample_id,treatmentid)%>%
    mutate(chem_name=tolower(treatmentid))%>%
    distinct()%>%
    left_join(ldmap)
  
  #all_ids<-unique(samps$other_id)
  #  dplyr::rename(other_id='sampleid')%>%##this maps it to the file
  
  ##now map the sample information to sample file. here, grep is hte most reliable
  
  
  ##get the raw data
  ##now 
  alldat <- sensitivityRaw(dset)
  
  doseDat<-alldat[,,1]%>%
    as.data.frame()%>%
    tibble::rownames_to_column('exp_id')%>%
    tidyr::pivot_longer(cols=starts_with('dose'),names_to='doseNum',values_to='Dose')
  
  
  respDat<-alldat[,,2]%>%
    as.data.frame()%>%
    tibble::rownames_to_column('exp_id')%>%
    tidyr::pivot_longer(cols=starts_with('dose'),names_to='doseNum',values_to='Response')
  
  
  doseRep<-doseDat%>%
    dplyr::full_join(respDat,by=c('doseNum','exp_id'))%>%
    left_join(full.map)%>%
    dplyr::select(DRUG=improve_chem_id,CELL=improve_sample_id,DOSE=Dose,GROWTH=Response)%>%
    #dplyr::mutate(DOSE=-log10(Dose/1000))###curve fitting code requires -log10(M), these are mM
    #rename(GROWTH=RESPONSE)%>%
    mutate(SOURCE='pharmacoGX')%>%
    mutate(STUDY=studyName)
  
  print(head(doseRep))
  
  ##now map  drugs, samples, genes to ids in database files
  
  write.table(doseRep,file=paste0(tolower(studyName),'DoseResponse'),sep='\t',row.names=F,quote=F)
  
}

###################################################################
##now we iterate through the drug files one by one
####




#' getCellLineData - gets cell line dose response data
getCellLineDoseData<-function(cell.lines=c('CTRPv2','FIMM','gCSI','PRISM','GDSC','NCI60','CCLE')){
  ###first get cell lines
  all.dose.rep<-do.call(rbind,lapply(cell.lines,function(cel){
  
   print(cel)
   files<-subset(all.dsets,`Dataset Name`==cel)%>%
      dplyr::select(`PSet Name`)%>%
      unlist()
  
    res<-do.call(rbind,lapply(files,function(f){
      print(f)
      if(f=='GDSC_2020(v2-8.2)')
        cel='GDSCv2'
      if(f=='GDSC_2020(v1-8.2)')
        cel='GDSCv1'
      tmpfile<-paste0(cel,'doseResponse')

      if(file.exists(tmpfile))
        dres<-read.table(tmpfile,sep='\t',header=T)
      else{
        dset<<-downloadPSet(f)
    
        url=subset(all.dsets,`PSet Name`==f)$Download
      #print(url)

        dres<-getDoseRespData(dset,cel)
      #p  rint(dres)
      }  
      return(dres)
      
    }))
      #print(res)
      return(res)
  }))
  
  #print(head(all.dose.rep))
  all.dose.rep%>%
    subset(!is.na(DOSE))
  #write.table(all.dose.rep,file='allDoseRepPreCalc.tsv',sep='\t')
}

if(FALSE){

  cl1<-c('CTRPv2','FIMM','gCSI','PRISM')
  cl2<-c('GDSC','NCI60','CCLE')
  dl1<-getCellLineDoseData(cl1)
  dl2<-getCellLineDoseData(cl2)

  write.table(dl1,'doseRep1.tsv',sep='\t',row.names=F,quote=F)
  write.table(dl2,'doseRep2.tsv',sep='\t',row.names=F,quote=F)
}

#' get cell line mutation data
getCellLineMutData<-function(){

  
  cmap<-improve_samples%>%
    dplyr::select(other_id,improve_sample_id)%>%
    distinct()
  
    ##GDSCv1
  dset<-downloadPSet('GDSC_2020(v1-8.2)')
  edat<-molecularProfiles(dset)$mutation
  sampmap<-colData(edat)%>%
    as.data.frame()%>%
    tibble::rownames_to_column('samples')%>%
    dplyr::select('samples','sampleid')%>%
    distinct()
  ##now we want to get the rna expression
  geneExp<-assays(edat)$exprs%>%
    as.data.frame()%>%
    tibble::rownames_to_column('gene_symbol')%>%
     tidyr::pivot_longer(cols=seq(2,1+ncol(assays(edat)$exprs)),
                        names_to='samples',
                        values_to='mutation')%>%
    left_join(sampmap)%>%
    subset(!is.na(mutation))
  
  ##now we map samples
  gdscV1Exp<-geneExp%>%
    dplyr::rename(other_id='sampleid')%>%
    left_join(cmap)%>%
    subset(!is.na(improve_sample_id))%>%
    select(gene_symbol,mutation,improve_sample_id)%>%
    distinct()%>%
#    dplyr::rename(other_id='gene')%>%
    left_join(improve_genes)%>%
    dplyr::select(entrez_id,improve_sample_id,mutation)%>%
    subset(!is.na(entrez_id))%>%
    mutate(source='PharmacoGX',study='GDSCv1')
  
  ##gdscv2
  dset<-downloadPSet('GDSC_2020(v2-8.2)')
  edat<-molecularProfiles(dset)$mutation
  sampmap<-colData(edat)%>%
    as.data.frame()%>%
    tibble::rownames_to_column('samples')%>%
    dplyr::select('samples','sampleid')%>%
    distinct()
  ##now we want to get the rna expression
  geneExp<-assays(edat)$exprs%>%
    as.data.frame()%>%
    tibble::rownames_to_column('gene_symbol')%>%
    tidyr::pivot_longer(cols=seq(2,1+ncol(assays(edat)$exprs)),
                        names_to='samples',
                        values_to='mutation')%>%
    left_join(sampmap)%>%
    subset(!is.na(mutation))
  
  ##now we map samples
  gdscV2Exp<-geneExp%>%
    dplyr::rename(other_id='sampleid')%>%
    left_join(cmap)%>%
    subset(!is.na(improve_sample_id))%>%
    select(gene_symbol,mutation,improve_sample_id)%>%
    distinct()%>%
    #    dplyr::rename(other_id='gene')%>%
    left_join(improve_genes)%>%
    dplyr::select(entrez_id,improve_sample_id,mutation)%>%
    subset(!is.na(entrez_id))%>%
    mutate(source='PharmacoGX',study='GDSCv2')
  
  
  ##gcsi
  dset<-downloadPSet('gCSI_2019')
  edat<-molecularProfiles(dset)$mutation
  sampmap<-colData(edat)%>%
    as.data.frame()%>%
    tibble::rownames_to_column('samples')%>%
    dplyr::select('samples','sampleid')%>%
    distinct()
  ##now we want to get the rna expression
  geneExp<-assays(edat)$exprs%>%
    as.data.frame()%>%
    tibble::rownames_to_column('gene_symbol')%>%
    tidyr::pivot_longer(cols=seq(2,1+ncol(assays(edat)$exprs)),
                        names_to='samples',
                        values_to='mutation')%>%
    left_join(sampmap)%>%
    subset(!is.na(mutation))
  
  ##now we map samples
  gcsiGeneExp<-geneExp%>%
    dplyr::rename(other_id='sampleid')%>%
    left_join(cmap)%>%
    subset(!is.na(improve_sample_id))%>%
    select(gene_symbol,mutation,improve_sample_id)%>%
    distinct()%>%
  #  dplyr::rename(other_id='gene')%>%
    left_join(improve_genes)%>%
    dplyr::select(entrez_id,improve_sample_id,mutation)%>%
    subset(!is.na(entrez_id))%>%
    mutate(source='PharmacoGX',study='gCSI') 
  
  ##CCLE
  dset<-downloadPSet('CCLE_2015')
  edat<-molecularProfiles(dset)$mutation
  sampmap<-colData(edat)%>%
    as.data.frame()%>%
    tibble::rownames_to_column('samples')%>%
    dplyr::select('samples','sampleid')%>%
    distinct()
  ##now we want to get the rna expression
  geneExp<-assays(edat)$exprs%>%
    as.data.frame()%>%
    tibble::rownames_to_column('gene_symbol')%>%
    tidyr::pivot_longer(cols=seq(2,1+ncol(assays(edat)$exprs)),
                        names_to='samples',
                        values_to='mutation')%>%
    left_join(sampmap)%>%
    subset(!is.na(mutation))
  
  ##now we map samples
  ccleExp<-geneExp%>%
    dplyr::rename(other_id='sampleid')%>%
    left_join(cmap)%>%
    subset(!is.na(improve_sample_id))%>%
    select(gene_symbol,mutation,improve_sample_id)%>%
    distinct()%>%
  #  dplyr::rename(other_id='gene')%>%
    left_join(improve_genes)%>%
    dplyr::select(entrez_id,improve_sample_id,mutation)%>%
    subset(!is.na(entrez_id))%>%
    mutate(source='PharmacoGX',study='CCLE') 
  
  fullmuts<-rbind(ccleExp,gcsiGeneExp,gdscV1Exp,gdscV2Exp)%>%
    unique()
  write.table(fullmuts,sep=',',row.names=F,col.names=T,file='data/mutations.csv')
}

#' get cell line copy number data
getCellLineCNVdata<-function(){
  
}

#' get cell line expression
getCellLineExpData<-function(cell.lines=c('gCSI','GDSC','NCI60','CCLE')){
  

  gex_val<-list(gCSI='Kallisto_0.46.1.rnaseq.counts',GDSCv2='Kallisto_0.46.1.rnaseq.counts',
                NCI60='rnaseq.iso',CCLE='Kallisto_0.46.1.rnaseq.counts')
  gene_val<-list(gCSI='ensgene',GDSCv2='ensgene',NCI60='gene_symbol',CCLE='ensgene')
  
  
  
  cmap<-improve_samples%>%
    dplyr::select(other_id,improve_sample_id)%>%
    distinct()
  
  
  ##GDSCv1
  dset<-downloadPSet('GDSC_2020(v1-8.2)')
  edat<-molecularProfiles(dset)$Kallisto_0.46.1.rnaseq.counts
  sampmap<-colData(edat)%>%
    as.data.frame()%>%
    tibble::rownames_to_column('samples')%>%
    dplyr::select('samples','sampleid')%>%
    distinct()
  ##now we want to get the rna expression
  geneExp<-assays(edat)$exprs%>%
    as.data.frame()%>%
    tibble::rownames_to_column('ensgene')%>%
    tidyr::separate(ensgene,into=c('gene','vers'),'\\.')%>%
    tidyr::pivot_longer(cols=seq(3,2+ncol(assays(edat)$exprs)),
                        names_to='samples',
                        values_to='counts')%>%
    left_join(sampmap)
  
  ##now we map samples
  gdscV1Exp<-geneExp%>%
    dplyr::rename(other_id='sampleid')%>%
    left_join(cmap)%>%
    subset(!is.na(improve_sample_id))%>%
    select(gene,counts,improve_sample_id)%>%
    distinct()%>%
#    dplyr::rename(other_id='gene')%>%
    left_join(improve_genes)%>%
    dplyr::select(entrez_id,improve_sample_id,counts)%>%
    subset(!is.na(entrez_id))%>%
    mutate(source='PharmacoGX',study='GDSCv1')
  
  
  ##GDSCv2
  dset<-downloadPSet('GDSC_2020(v2-8.2)')
  edat<-molecularProfiles(dset)$Kallisto_0.46.1.rnaseq.counts
  sampmap<-colData(edat)%>%
    as.data.frame()%>%
    tibble::rownames_to_column('samples')%>%
    dplyr::select('samples','sampleid')%>%
    distinct()
  ##now we want to get the rna expression
  geneExp<-assays(edat)$exprs%>%
    as.data.frame()%>%
    tibble::rownames_to_column('ensgene')%>%
    tidyr::separate(ensgene,into=c('gene','vers'),'\\.')%>%
    tidyr::pivot_longer(cols=seq(3,2+ncol(assays(edat)$exprs)),
                        names_to='samples',
                        values_to='counts')%>%
    left_join(sampmap)
  
  ##now we map samples
  gdscV2Exp<-geneExp%>%
    dplyr::rename(other_id='sampleid')%>%
    left_join(cmap)%>%
    subset(!is.na(improve_sample_id))%>%
    select(gene,counts,improve_sample_id)%>%
    distinct()%>%
    dplyr::rename(other_id='gene')%>%
    left_join(improve_genes)%>%
    dplyr::select(entrez_id,improve_sample_id,counts)%>%
    subset(!is.na(entrez_id))%>%
    mutate(source='PharmacoGX',study='GDSCv2')
  
  dset<-downloadPSet('CCLE_2015')
  
  ##CCLE only
  edat<-molecularProfiles(dset)$Kallisto_0.46.1.rnaseq.counts
  sampmap<-colData(edat)%>%
    as.data.frame()%>%
    tibble::rownames_to_column('samples')%>%
    dplyr::select('samples','sampleid')%>%
    distinct()
  ##now we want to get the rna expression
  geneExp<-assays(edat)$exprs%>%
    as.data.frame()%>%
    tibble::rownames_to_column('ensgene')%>%
    tidyr::separate(ensgene,into=c('gene','vers'),'\\.')%>%
    tidyr::pivot_longer(cols=seq(3,2+ncol(assays(edat)$exprs)),
                        names_to='samples',
                        values_to='counts')%>%
    left_join(sampmap)
  
  ##now we map samples
  ccleExp<-geneExp%>%
    dplyr::rename(other_id='sampleid')%>%
    left_join(cmap)%>%
    subset(!is.na(improve_sample_id))%>%
    select(gene,counts,improve_sample_id)%>%
    distinct()%>%
    dplyr::rename(other_id='gene')%>%
    left_join(improve_genes)%>%
    dplyr::select(entrez_id,improve_sample_id,counts)%>%
    subset(!is.na(entrez_id))%>%
    mutate(source='PharmacoGX',study='CCLE')
  
  #now we map gene names
 # write.table(geneExp2,file='ccleGeneExp.csv',sep=',',row.names=F,quote=F)
  
  
  ##NCI60 expression data
  dset<-PharmacoGx::downloadPSet('NCI60_2021')
  
  ##now we want to get the rna expression
  geneExp<-assays(molecularProfiles(dset)$rnaseq.iso)$exprs%>%
    as.data.frame()%>%
    tibble::rownames_to_column('gene')%>%
    tidyr::pivot_longer(cols=seq(2,1+ncol(assays(molecularProfiles(dset)$rnaseq.iso)$exprs)),
                        names_to='samples',
                        values_to='counts')
  
  ##now we map samples
  nciExp<-geneExp%>%
    dplyr::rename(other_names='samples')%>%
    left_join(improve_samples)%>%
    subset(!is.na(improve_sample_id))%>%
    select(gene,counts,improve_sample_id)%>%
    distinct()%>%
    dplyr::rename(gene_symbol='gene')%>%
    left_join(improve_genes)%>%
    dplyr::select(entrez_id,improve_sample_id,counts)%>%
    subset(!is.na(entrez_id))%>%
    mutate(source='PharmacoGX',study='NCI60')
  
  ##rbind data and then save
  
  
  ##gCSI expression data
  dset<-downloadPSet('gCSI_2019')
  edat<-molecularProfiles(dset)$Kallisto_0.46.1.rnaseq.counts
  sampmap<-colData(edat)%>%
    as.data.frame()%>%
    tibble::rownames_to_column('samples')%>%
    dplyr::select('samples','sampleid')%>%
    distinct()
  ##now we want to get the rna expression
  geneExp<-assays(edat)$exprs%>%
    as.data.frame()%>%
    tibble::rownames_to_column('ensgene')%>%
    tidyr::separate(ensgene,into=c('gene','vers'),'\\.')%>%
    tidyr::pivot_longer(cols=seq(3,2+ncol(assays(edat)$exprs)),
                        names_to='samples',
                        values_to='counts')%>%
    left_join(sampmap)

  
  ##now we map samples
  gCSIExp<-geneExp%>%
    dplyr::rename(other_id='sampleid')%>%
    left_join(cmap)%>%
    subset(!is.na(improve_sample_id))%>%
    select(gene,counts,improve_sample_id)%>%
    distinct()%>%
    dplyr::rename(other_id='gene')%>%
    left_join(improve_genes)%>%
    dplyr::select(entrez_id,improve_sample_id,counts)%>%
    subset(!is.na(entrez_id))%>%
    mutate(source='PharmacoGX',study='gCSI')
  
  cellLineGex<-rbind(gdscV2Exp,ccleExp,nciExp,gCSIExp)%>%
    unique()
  
  write.table(cellLineGex,file='expression.csv',sep=',',quote=F,row.names=F)
  return(cellLineGex)
}

 

#third get other data

