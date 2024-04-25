
library(RColorBrewer)
library(dplyr)
library(ggplot2)

mergeSamples<-function(){
  
  ##########################
  ## CANCER NAME MAPPING FILE
  ##this file tries to map as many cancer types as possible to groups that span model systems
  ###########################
  cmaps<-readr::read_csv('cellLineTypes.csv')
  
  ###########################
  ## CPTAC SAMPLE DATA
  ## Primary task here is to fix the capitalizations
  ###########################
  cptac<-readr::read_csv('cptac_samples.csv')|>
    mutate(cancer_type=stringr::str_replace_all(cancer_type,'Head and Neck','Head and neck'))|>
    mutate(cancer_type=stringr::str_replace_all(cancer_type,'Colon','Colorectal'))|>
    mutate(cancer_type=stringr::str_replace_all(cancer_type,'Uterine Corpus Endometrial Carcinoma','Uterine corpus endometrial carcinoma'))|>
    dplyr::mutate(`CPTAC Cancer type`=cancer_type)|>
    left_join(cmaps)|>
    mutate(sampleSource='CPTAC')|>
    dplyr::select(improve_sample_id,`CPTAC Cancer type`,model_type,species,sampleSource)|>
    distinct()
  
  ###########################
  ## Broad Sanger data
  ## We have many more cancer types here, so we try to map what we have to CPTAC names and then adjust the rest
  ##
  ###########################
  broad_sanger<-readr::read_csv('broad_sanger_samples.csv')|>
    dplyr::mutate(`Cell line cancer type`=cancer_type)|>
    mutate(sampleSource='CCLE')
  
  allec<-grep('Endometrial',broad_sanger$`Cell line cancer type`)
  broad_sanger$`Cell line cancer type`[allec]<-'Uterine corpus endometrial carcinoma'
  
  broad_sanger<-broad_sanger|>
    left_join(cmaps)
  
  ##first we collect the names of the cancers that are NOT in CPTAC
  other_cans<-which(is.na(broad_sanger$`CPTAC Cancer type`))
  broad_sanger$`CPTAC Cancer type`[other_cans]<-broad_sanger$`Cell line cancer type`[other_cans]
  broad_sanger<-broad_sanger|>
    dplyr::select(improve_sample_id,`CPTAC Cancer type`,model_type,species,sampleSource)|>
    distinct()
  
  #then we rename the NA values to 'Other' if we want
  #other_cans<-which(is.na(broad_sanger$`CPTAC Cancer type`))
  #broad_sanger$`CPTAC Cancer type`[other_cans]<-'Other'
  # or just remove them
  broad_sanger<-broad_sanger|>
    subset(!is.na(`CPTAC Cancer type`))
  
  ###########################
  ## HCMI SAMPLE DATA
  ## Here we lean heavily on the sample mapping file
  ###########################
  hcmi<-readr::read_csv('hcmi_samples.csv')|>
    dplyr::rename(id_source='other_id_source')|>
    mutate(species='human')|>
    subset(model_type%in%c('organoid','tumor','cell line','Patient derived xenograft'))|>
    dplyr::mutate(`HCMI Cancer type`=cancer_type,`HCMI Common name`=common_name)|>
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
  baml<-readr::read_csv("beataml_samples.csv")|>
    mutate(cancer_type='Acute myeloid leukemia')|>
    mutate(species='Human')|>
    mutate(model_type='tumor')|>
    mutate(sampleSource='BeatAML')|>
    dplyr::select(improve_sample_id,species,cancer_type,sampleSource,model_type)|>
    distinct()
  ###########################
  ## MPNST SAMPLE DATA
  ## TBD
  ###########################
  mpnst<-readr::read_csv("mpnst_samples.csv")|>
    mutate(cancer_type='Neurofibromatosis')|>
    mutate(species='Human')|>
    mutate(sampleSource='MPNST')|>
    dplyr::select(improve_sample_id,species,cancer_type,sampleSource,model_type)|>
    distinct()
  
  ##now we join thomdelsem into a single table, with cancer type
  fulldat<<-rbind(cptac,broad_sanger,hcmi)|>
    dplyr::rename(cancer_type=`CPTAC Cancer type`)|>
    subset()
  
  fulldat<-fulldat|>
    dplyr::select(improve_sample_id,species,cancer_type,sampleSource,model_type)|>
    distinct()|>
    rbind(baml)
  
  fulldat<-fulldat|>
    dplyr::select(improve_sample_id,species,cancer_type,sampleSource,model_type)|>
    distinct()|>
    rbind(mpnst)
  
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
  
  return(rbind(alldat,other_can))
}

fulldat<-mergeSamples()

##looking for exact matches
stats<-fulldat|>
  subset(!is.na(cancer_type))|>
  group_by(cancer_type,model_type,sampleSource)|>
  summarize(numSamps=n_distinct(improve_sample_id))|>
  subset(model_type!='Not Reported')|>
  subset(numSamps>1)

color_palette <- brewer.pal(n = 4, name = "Set2")

# Assign colors to the model types
names(color_palette) <- c("tumor", "cell line", "organoid",'Patient derived xenograft')

fig0<-ggplot(stats,aes(x=cancer_type,y=numSamps,fill=model_type))+
  geom_bar(stat='identity',position='dodge')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
  plot.background = element_rect(fill = background_color, color = background_color),
  legend.background = element_rect(fill = background_color, color = background_color))+
  scale_y_log10()+scale_fill_manual(values=color_palette)+
  ggtitle('Samples by tumor type')

print(fig0)
ggsave('Fig0_Overview.png',fig0,height=8,width=10)

# Subset data for each type
data_type1 <- subset(stats, sampleSource == 'HCMI')
data_type2 <- subset(stats, sampleSource == 'BeatAML')
data_type3 <- subset(stats, sampleSource == 'CPTAC')
data_type4 <- subset(stats, sampleSource == 'CCLE')
data_type5 <- subset(stats, sampleSource == 'MPNST')

# Create separate plots for each type, with colorblind-friendly colors
background_color <- "#E0F2F1"
fig1 <- ggplot(data_type1, aes(x=cancer_type, y=numSamps, fill=model_type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values=color_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
  plot.background = element_rect(fill = background_color, color = background_color),
  legend.background = element_rect(fill = background_color, color = background_color)) +
  ggtitle('Cancer and Tissue Types - HCMI')

fig2 <- ggplot(data_type2, aes(x=cancer_type, y=numSamps, fill=model_type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values=color_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
  plot.background = element_rect(fill = background_color, color = background_color),
  legend.background = element_rect(fill = background_color, color = background_color)) +
  ggtitle('Cancer and Tissue Types - BeatAML')

fig3 <- ggplot(data_type3, aes(x=cancer_type, y=numSamps, fill=model_type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values=color_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
  plot.background = element_rect(fill = background_color, color = background_color),
  legend.background = element_rect(fill = background_color, color = background_color))+
  ggtitle('Cancer and Tissue Types - CPTAC')

fig4 <- ggplot(data_type4, aes(x=cancer_type, y=numSamps, fill=model_type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values=color_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
  plot.background = element_rect(fill = background_color, color = background_color),
  legend.background = element_rect(fill = background_color, color = background_color)) +
  scale_y_log10() +
  ggtitle('Cancer and Tissue Types - Broad Sanger')

fig5 <- ggplot(data_type5, aes(x=cancer_type, y=numSamps, fill=model_type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values=color_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
  plot.background = element_rect(fill = background_color, color = background_color),
  legend.background = element_rect(fill = background_color, color = background_color)) +
  ggtitle('Cancer and Tissue Types - MPNST')

ggsave('Fig1_HCMI.png', fig1, height=8, width=10)
ggsave('Fig2_BeatAML.png', fig2, height=8, width=10)
ggsave('Fig3_CPTAC.png', fig3, height=8, width=10)
ggsave('Fig4_Broad_Sanger.png', fig4, height=8, width=10)
ggsave('Fig5_MPNST.png', fig5, height=8, width=10)


