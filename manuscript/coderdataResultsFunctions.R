##figure ploting functions incorporates and standardizes the calls made across figures

library(ggplot2)
library(dplyr)
library(arrow)
library(ggridges)
library(synapser)
library(RColorBrewer)
##COLORS: standardize here

modelcolors <- RColorBrewer::brewer.pal(n=6,name='RdYlBu')
names(modelcolors) <- c('deepttc','graphdrp','lgbm','pathdsp','uno')


exvivo = c('mpnst','beataml','sarcpdo','pancpdo','bladderpdo','liverpdo')
cellline = c('nci60','ctrpv2','fimm','gcsi','gdscv1','gdscv2','prism','ccle')

ccols = RColorBrewer::brewer.pal(n = length(cellline),name = 'RdBu')
names(ccols) <- cellline

ecols = RColorBrewer::brewer.pal(n = length(exvivo),name = 'PRGn')
names(ecols) <- exvivo

datasetcolors <- c(ccols,ecols)

synapser::synLogin()

getProteomicsData <- function(){

  allscoreslist <- list(lgbm = 'syn68156330')
  fullres <- do.call(rbind,lapply(names(allscoreslist),function(mod)
    readr::read_csv(synapser::synGet(allscoreslist[[mod]])$path) |> mutate(model = mod)))

  fullres <- fullres |>
    mutate(withinDataset = ifelse(src == trg,TRUE,FALSE))

  #lets remove same-dataset data
 # cdres <- subset(fullres,!withinDataset)

  ##lets remove ex vivo training
  fullres <- subset(fullres,!src %in% c('mpnst','beataml'))

  return(fullres)

}

getModelPerformanceData <- function(){

  allscoreslist <- list(deepttc = 'syn65880080',graphdrp = 'syn65928973',lgbm = 'syn65880116',pathdsp = 'syn65880133',uno = 'syn65676159')

  ##with pancpdo
  allscoreslist <- list(deepttc = 'syn66323471',graphdrp = 'syn66323492',lgbm = 'syn66323510',pathdsp = 'syn66326173',uno = 'syn66323527')

  allscoreslist <- list(deepttc = 'syn68155944', graphdrp = 'syn68155957', lgbm = 'syn68155973', pathdsp = 'syn66326173', uno = 'syn68155991')
  fullres <- do.call(rbind,lapply(names(allscoreslist),function(mod)
    readr::read_csv(synapser::synGet(allscoreslist[[mod]])$path) |> mutate(model = mod)))

  fullres <- fullres |>
    mutate(withinDataset = ifelse(src == trg,TRUE,FALSE))

  #lets remove same-dataset data
  cdres <- subset(fullres,!withinDataset)

  ##lets remove ex vivo training
  cdres <- subset(cdres,!src %in% c('mpnst','beataml'))

  return(cdres)
}


# this currently only retrieves one dataset at a time and returns an appache
# "arrow" tabular dataset object that can be interacted / queried via dplyr
getModelPredictionData <- function(dset='lgbm') {

  preds <- list(lgbm = 'syn68176033')

  dataset <- arrow::open_dataset(
    sources = synapser::synGet(preds[[dset]])$path,
    format = "parquet"
    )

  return(dataset)
}

###these files are very big so i'm not sure how to deal with them.
# getModelPredictionData <- function(dset='lgbm') {
#
#   preds <- list(deepttc = 'syn68149793', graphdrp = 'syn68146828', lgbm = 'syn68149807', pathdsp = 'syn66772452', uno = 'syn68149809')
#
#   fullres <- do.call(rbind,lapply(dset,function(mod)
#     readr::read_csv(synapser::synGet(preds[[mod]])$path) |> mutate(model = mod)))
#
#   return(preds)
# }

#this function plots a single metric by all the possible values
#
ridgelineMetricPlots <- function(metric,dataset=cdres, prefix='all'){

    sr <- dataset |>
    subset(met == metric)


  ##facet by source - compare performance across a single source

  ##re-rank src samples by mean metrics
  mvals <- sr |> group_by(src) |>
    summarize(mvals = mean(value)) |>
    arrange(mvals)

  if (metric == 'r2') {
    sr$value <- sapply(sr$value,function(x) ifelse(x < (-1),-1,x))
  }

  sr$src = factor(sr$src,levels = mvals$src)

  #compare models by source dataset
  p1 <- sr |>
    ggplot(aes(x = value,y = trg,fill = model)) +
    ggridges::geom_density_ridges(alpha = 0.5) +
    facet_grid(src~.) +
    ggtitle(paste0(metric,' by source dataset'))+
    scale_fill_manual(values=modelcolors)

  ##now we rerank by target dataset and evaluate by target
  mvals <- sr |> group_by(trg) |>
    summarize(mvals = mean(value)) |>
    arrange(mvals)
  sr$trg = factor(sr$trg,levels = mvals$trg)

  #plot source by target data
  p3 <- sr |>
    ggplot(aes(x = value,y = src,fill = model)) +
    ggridges::geom_density_ridges(alpha = 0.5) +
    facet_grid(trg~.) +
    ggtitle(paste0(metric,' by target dataset'))+
    scale_fill_manual(values=modelcolors)


  return(list(src=p1,trg=p3))
}


##here we have to interrogate the results to visualize how specific drugs are behaving
performanceByDrugOrSample <- function(){

}

## calculate source dataset statistics
## how do features of the source dataset impact  performance?
calcSourceStatistics<-function(metric, dataset=cdres){
  #number of combos
  combos = c(ccle = 10911,ctrpv2 = 303520,fimm = 2457 ,gcsi = 12320,
             gdscv1 = 105808,gdscv2 = 45323, nci60 = 2317205,prism = 633169)

  numsamples = c(ccle = 503, ctrpv2 = 847, fimm=52, gcsi = 571, gdscv1 = 984,
                 gdscv2 = 806, nci60 = 83, prism = 478)
  numdrugs = c(ccle = 24, ctrpv2 = 461, fimm=52, gcsi = 43, gdscv1 = 296, gdscv2 = 169,
               nci60 = 54707, prism = 1418)

  stats = data.frame(Samples = numsamples, Drugs =numdrugs, Combos = combos)
  stats$src = rownames(stats)

  #todo: we can also evaluate number of samples or drugs

  #e can get performance summaries
  gres <- dataset  |>
    #subset(model!='uno')|>
    subset(met==metric) |>
    group_by(met,src,trg,model) |>
    summarize(meanVal=mean(value,na.rm=TRUE)) |>
    left_join(stats) |>
    arrange(meanVal)

 # mom <- gres|>
#    group_by(src,Combos)|>
#    summarize(mv = mean(meanVal))|>
#    arrange(mv)

#  gres$src = factor(gres$src,levels = unique(mom$src))

  p1 <- ggplot(gres, aes(x=Samples, y = meanVal,col=model))+
    geom_point()+scale_x_log10()+scale_color_manual(values=modelcolors)+geom_smooth(method=lm, alpha=0.2)+theme_bw()

  p2 <- ggplot(gres, aes(x=Drugs, y = meanVal,col=model))+
    geom_point()+scale_x_log10()+scale_color_manual(values=modelcolors)+geom_smooth(method=lm, alpha=0.2)+theme_bw()

  p3 <- ggplot(gres, aes(x=Combos, y = meanVal,col=model))+
    geom_point()+scale_x_log10()+scale_color_manual(values=modelcolors)+geom_smooth(method=lm, alpha=0.2)+theme_bw()

  corvals <- gres |>
    ungroup() |>
    group_by(model) |>
    summarize(Sample=cor(Samples,meanVal),Drugs=cor(Drugs,meanVal),Combinations=cor(Combos,meanVal,use='pairwise.complete.obs'))|>
    tidyr::pivot_longer(cols=c(2,3,4),names_to='statistic',values_to='correlation')

  p4 <- ggplot(corvals, aes(x=statistic,y=correlation,fill=model)) + geom_bar(position='dodge',stat='identity') +
    scale_fill_manual(values=modelcolors) + theme_bw()

  return(cowplot::plot_grid(p1,p2,p3,p4,nrow=4))
}

##do we still need this function?

doModelPlot <- function(metric, dataset=cdres){

  sr <- dataset |>
    subset(met == metric)
  ##re-rank src samples by mean metric
  mvals <- sr |>
    group_by(trg) |>
    summarize(mvals = mean(value)) |>
    arrange(mvals)

  if(metric == 'r2') {
    sr$value <- sapply(sr$value,function (x) ifelse(x<(-1),-1,x))
  }

  sr$trg = factor(sr$trg,levels = mvals$trg)

  sr |>
    subset(trg %in% exvivo) |>
    ggplot(aes(x = value,alpha = 0.8)) +
    geom_histogram() + facet_grid(model~trg) +
    ggtitle(paste0(metric,' evaluated on ex vivo data'))


  ggsave(paste0(metric,'exVivoPerformance.png'))


  sr |> subset(!trg %in% exvivo) |>
    ggplot(aes(x = value,alpha = 0.8)) +
    geom_histogram() + facet_grid(model~trg) +
    ggtitle(paste0(metric,' evaluated on cell line data'))


  ggsave(paste0(metric,'CellLinePerformance.png'))

}
