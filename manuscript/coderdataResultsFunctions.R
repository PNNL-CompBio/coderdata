##figure ploting functions incorporates and standardizes the calls made across figures

library(ggplot2)
library(dplyr)
library(ggridges)
library(synapser)

##COLORS: standardize here
modelcolors <- c()
datasetcolors <- c()

exvivo = c('mpnst','beataml','sarcpdo','pancpdo','bladderpdo')

synapser::synLogin()

getModelPerformanceData <- function(){

  allscoreslist <- list(deepttc = 'syn65880080',graphdrp = 'syn65928973',lgbm = 'syn65880116',pathdsp = 'syn65880133',uno = 'syn65676159')

  ##with pancpdo
  allscoreslist <- list(deepttc = 'syn66323471',graphdrp = 'syn66323492',lgbm = 'syn66323510',pathdsp = 'syn66326173',uno = 'syn66323527')

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



###these files are very big so i'm not sure how to deal with them.
getModelPredictionData<-function(dset='lgbm'){

  preds <- list(deepttc = 'syn68149793', graphdrp = 'syn68146828', lgbm = 'syn68149807', pathdsp = 'syn66772452', uno = 'syn68149809')


  fullres <- do.call(rbind,lapply(dset,function(mod)
    readr::read_csv(synapser::synGet(preds[[mod]])$path) |> mutate(model = mod)))

  return(preds)
}

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
    ggtitle(paste0(metric,' by source dataset'))

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
    ggtitle(paste0(metric,' by target dataset'))

  return(list(src=p1,trg=p3))
}


##here we have to interrogate the results to visualize how specific drugs are behaving
performanceByDrugOrSample<-function(){

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
