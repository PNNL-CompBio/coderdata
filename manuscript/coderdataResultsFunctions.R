##figure ploting functions incorporates and standardizes the calls made across figures

library(ggplot2)
library(dplyr)
library(arrow)
library(ggridges)
library(synapser)
library(RColorBrewer)
##COLORS: standardize here

models <- c('deepttc','graphdrp','lgbm','pathdsp','uno')
modelcolors <- RColorBrewer::brewer.pal(n=6,name='RdYlBu')
names(modelcolors) <- models


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
#
#
# getProteomicsData <- function(){
#
#   allscoreslist <- list(lgbm = 'syn68156330')
#   fullres <- do.call(rbind,lapply(names(allscoreslist),function(mod)
#     readr::read_csv(synapser::synGet(allscoreslist[[mod]])$path) |> mutate(model = mod)))
#
#   fullres <- fullres |>
#     mutate(withinDataset = ifelse(src == trg,TRUE,FALSE))
#
#   #lets remove same-dataset data. - I'm going to keep it actually...
#  # cdres <- subset(fullres,!withinDataset)
#
#   ##lets remove ex vivo training
#   fullres <- subset(fullres,!src %in% c('mpnst','beataml'))
#
#   return(fullres)
#
# }
#
 getModelPerformanceData <- function(include_within = TRUE){
   allscoreslist <- list(deepttc='syn69960485', graphdrp='syn69953705',
                         lgbm='syn69931939', xgboostdrp='syn69931970', uno='syn69953101')

   fullres <- purrr::imap_dfr(allscoreslist, \(id, mod)
     readr::read_csv(synapser::synGet(id)$path, show_col_types = FALSE) |> mutate(model = mod))

   fullres <- fullres |>
     mutate(withinDataset = src == trg)

   cdres <- if (include_within) fullres else subset(fullres, !withinDataset)

   cdres <- subset(cdres, !src %in% c('mpnst','beataml'))
   cdres
 }
#
#
#
# Reticulate vs Arrow - this is where it matters.
if (use_reticulate) {
  getModelPredictionData <- function(dset = "lgbm") {
    preds <- list(
      deepttc="syn69960497", graphdrp="syn69953716", lgbm="syn69931953",
      xgboostdrp="syn69931981", uno="syn69953114"
    )
    p <- synapser::synGet(preds[[dset]])$path
    if (dir.exists(p)) p <- file.path(p, "*.parquet")

    pd <- import("pandas", convert = FALSE)
    df_py <- pd$read_parquet(p)
    as_tibble(reticulate::py_to_r(df_py))
  }
} else {
  getModelPredictionData <- function(dset='lgbm') {
    preds <- list(
      deepttc = "syn69960497",
      graphdrp = "syn69953716",
      lgbm = "syn69931953",
      xgboostdrp = "syn69931981",
      uno = "syn69953114"
    )

    dataset <- arrow::open_dataset(
      sources = synapser::synGet(preds[[dset]])$path,
      format = "parquet"
      )

    return(dataset)
  }
}
# getModelPerformanceData <- function(){
#
#   allscoreslist <- list(deepttc = 'syn65880080',graphdrp = 'syn65928973',lgbm = 'syn65880116',pathdsp = 'syn65880133',uno = 'syn65676159')
#
#   ##with pancpdo
#   allscoreslist <- list(deepttc = 'syn66323471',graphdrp = 'syn66323492',lgbm = 'syn66323510',pathdsp = 'syn66326173',uno = 'syn66323527')
#
#   allscoreslist <- list(deepttc = 'syn68155944', graphdrp = 'syn68155957', lgbm = 'syn68155973', pathdsp = 'syn66326173', uno = 'syn68155991')
#   fullres <- do.call(rbind,lapply(names(allscoreslist),function(mod)
#     readr::read_csv(synapser::synGet(allscoreslist[[mod]])$path) |> mutate(model = mod)))
#
#   fullres <- fullres |>
#     mutate(withinDataset = ifelse(src == trg,TRUE,FALSE))
#
#   #lets remove same-dataset data
#   cdres <- subset(fullres,!withinDataset)
#
#   ##lets remove ex vivo training
#   cdres <- subset(cdres,!src %in% c('mpnst','beataml'))
#
#   return(cdres)
# }

#
# # this currently only retrieves one dataset at a time and returns an appache
# # "arrow" tabular dataset object that can be interacted / queried via dplyr
# getModelPredictionData <- function(dset='lgbm') {
#
#   preds <- list(
#     deepttc = "syn68176968",
#     graphdrp = "syn68176977",
#     lgbm = "syn68176033",
#     pathdsp = "syn68176970",
#     uno = "syn68176971"
#   )
#
#   dataset <- arrow::open_dataset(
#     sources = synapser::synGet(preds[[dset]])$path,
#     format = "parquet"
#     )
#
#   return(dataset)
# }


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


paired_heatmap_triangular_mean_by_target_4x4 <- function(
    cdres,
    metric = "scc",
    model_order = c("deepttc","graphdrp","lgbm","xgboostdrp","uno"),
    limits = c(-0.15, 0.15),
    show_labels = TRUE,
    diag_mode = c("na","zero")
) {
  diag_mode <- match.arg(diag_mode)
  have_models <- intersect(model_order, unique(cdres$model))
  stopifnot(length(have_models) >= 2)

  dom_of_trg <- function(trg) ifelse(trg %in% exvivo, "ex vivo", "cell line")

  # Per-(model, src, trg) mean + domain
  cell_means <- cdres %>%
    dplyr::filter(met == metric, src != trg, model %in% have_models) %>%
    dplyr::group_by(model, src, trg) %>%
    dplyr::summarise(mu = mean(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(trg_domain = dom_of_trg(trg))

  # ---- Per-target average (for annotation) ----
  # avg over models and sources for that target
  avg_by_trg <- cell_means %>%
    dplyr::group_by(trg, trg_domain) %>%
    dplyr::summarise(avg_mu = mean(mu, na.rm = TRUE), .groups = "drop")

  # Helper to compute pairwise Δ mean per target (fill only domain-matching triangle)
  compute_for_target <- function(df_trg) {
    trgnm <- df_trg$trg[1]; dom <- df_trg$trg_domain[1]
    wide <- tidyr::pivot_wider(df_trg, names_from = model, values_from = mu)
    prs  <- t(combn(have_models, 2))

    rows <- apply(prs, 1, function(pr) {
      a <- pr[1]; b <- pr[2]
      sub <- wide %>% dplyr::select(src, !!rlang::sym(a), !!rlang::sym(b)) %>% tidyr::drop_na()
      if (nrow(sub) < 3) return(NULL)
      x <- sub[[a]]; y <- sub[[b]]
      d  <- y - x
      md <- mean(d, na.rm = TRUE)
      tt <- tryCatch(stats::t.test(y, x, paired = TRUE), error = function(e) NULL)
      tibble::tibble(trg = trgnm, trg_domain = dom,
                     a = a, b = b, MD = md,
                     p = if (!is.null(tt)) tt$p.value else NA_real_,
                     n = sum(stats::complete.cases(d)))
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (!length(rows)) return(NULL)
    out <- dplyr::bind_rows(rows)

    # Tile grid for this target
    grid <- expand.grid(row_model = have_models, col_model = have_models, stringsAsFactors = FALSE)
    grid$row_model <- factor(grid$row_model, levels = have_models)
    grid$col_model <- factor(grid$col_model, levels = have_models)
    grid$i <- as.integer(grid$row_model); grid$j <- as.integer(grid$col_model)

    # Δ mean = (Column − Row), using either (a=row,b=col) or flipped
    get_md <- function(row, col) {
      m1 <- out %>% dplyr::filter(a == row, b == col)
      if (nrow(m1) == 1) return(c(MD = -m1$MD, p = m1$p))   # flip sign
      m2 <- out %>% dplyr::filter(a == col, b == row)
      if (nrow(m2) == 1) return(c(MD = m2$MD, p = m2$p))    # flip sign
      c(MD = NA_real_, p = NA_real_)
    }

    purrr::pmap_dfr(
      list(grid$i, grid$j, as.character(grid$row_model), as.character(grid$col_model)),
      function(i, j, rm, cm) {
        if (i == j) {
          if (diag_mode == "zero") {
            tibble::tibble(trg = trgnm, trg_domain = dom,
                           row_model = rm, col_model = cm,
                           MD = 0, p = NA_real_, which_tri = "diag")
          } else {
            tibble::tibble(trg = trgnm, trg_domain = dom,
                           row_model = rm, col_model = cm,
                           MD = NA_real_, p = NA_real_, which_tri = "diag")
          }
        } else if (i < j && dom == "cell line") {
          mdp <- get_md(rm, cm)
          tibble::tibble(trg = trgnm, trg_domain = dom,
                         row_model = rm, col_model = cm,
                         MD = as.numeric(mdp["MD"]), p = as.numeric(mdp["p"]),
                         which_tri = "upper_cellline")
        } else if (i > j && dom == "ex vivo") {
          mdp <- get_md(rm, cm)
          tibble::tibble(trg = trgnm, trg_domain = dom,
                         row_model = rm, col_model = cm,
                         MD = as.numeric(mdp["MD"]), p = as.numeric(mdp["p"]),
                         which_tri = "lower_exvivo")
        } else {
          tibble::tibble(trg = trgnm, trg_domain = dom,
                         row_model = rm, col_model = cm,
                         MD = NA_real_, p = NA_real_, which_tri = "offdomain")
        }
      }
    )
  }

  mats <- cell_means %>%
    dplyr::group_by(trg) %>%
    tidyr::nest() %>%
    dplyr::mutate(res = purrr::map(data, compute_for_target)) %>%
    tidyr::unnest(res) %>%
    dplyr::mutate(
      row_model = factor(row_model, levels = have_models),
      col_model = factor(col_model, levels = have_models)
    )

  # Drop NA tiles → blanks for off-domain cells
  mats_plot <- mats %>% dplyr::filter(is.finite(MD))

  # Facet order: cell-line first, then ex vivo
  if (exists("cellline") && exists("exvivo")) {
    facet_levels <- c(cellline, exvivo)
  } else {
    trg_domains <- mats %>% dplyr::distinct(trg, trg_domain)
    facet_levels <- c(
      trg_domains %>% dplyr::filter(trg_domain == "cell line") %>% dplyr::arrange(trg) %>% dplyr::pull(trg),
      trg_domains %>% dplyr::filter(trg_domain == "ex vivo")   %>% dplyr::arrange(trg) %>% dplyr::pull(trg)
    )
  }
  mats_plot$trg <- factor(mats_plot$trg, levels = facet_levels)
  avg_by_trg$trg <- factor(avg_by_trg$trg, levels = facet_levels)

  # Build corner-annotation data:
  #  - cell line: top-left -> (x=first model, y=first model)
  #  - ex vivo:   bottom-right -> (x=last model,  y=last model)
  # Build corner-annotation data (fix y positions):
  #  - cell line: top-left  -> x = first model, y = LAST model
  #  - ex vivo:   bottom-right -> x = last model,  y = FIRST model
  first_m <- have_models[1]
  last_m  <- tail(have_models, 1)

  ann_cell <- avg_by_trg %>%
    dplyr::filter(trg_domain == "cell line") %>%
    dplyr::transmute(
      trg,
      x = factor(first_m, levels = have_models),   # left
      y = factor(last_m,  levels = have_models),   # top
      lab = sprintf("Avg=%.3f", avg_mu),
      hjust = 0, vjust = 1.1
    )

  ann_ex <- avg_by_trg %>%
    dplyr::filter(trg_domain == "ex vivo") %>%
    dplyr::transmute(
      trg,
      x = factor(last_m,  levels = have_models),   # right
      y = factor(first_m, levels = have_models),   # bottom
      lab = sprintf("Avg=%.3f", avg_mu),
      hjust = 1, vjust = -0.1
    )


  # Labels inside tiles
  star <- function(p) dplyr::case_when(
    is.na(p) ~ "",
    p < 1e-4 ~ "****",
    p < 1e-3 ~ "***",
    p < 1e-2 ~ "**",
    p < 0.05 ~ "*",
    TRUE     ~ ""
  )
  mats_plot$lab <- if (show_labels) sprintf("%.3f%s", mats_plot$MD, star(mats_plot$p)) else ""

  ggplot(mats_plot, aes(x = col_model, y = row_model, fill = MD)) +
    # Tile grid
    geom_tile(color = "white", linewidth = 0.25, width = 0.999, height = 0.999) +
    {if (show_labels) geom_text(aes(label = lab), size = 2.0)} +
    # Corner annotations (do not inherit tile aesthetics)
    geom_text(data = ann_cell, aes(x = x, y = y, label = lab, hjust = hjust, vjust = vjust),
              inherit.aes = FALSE, size = 3, fontface = "bold") +
    geom_text(data = ann_ex,   aes(x = x, y = y, label = lab, hjust = hjust, vjust = vjust),
              inherit.aes = FALSE, size = 3, fontface = "bold") +
    scale_fill_distiller(
      palette = "RdBu", direction = 1,
      limits = limits, oob = scales::squish,
      name = "Mean Difference (row − col)") +
    coord_fixed() +
    theme_bw(base_size = 9) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      strip.text = element_text(face = "bold", size = 9, margin = margin(2,2,2,2)),
      panel.grid = element_blank(),
      panel.spacing = grid::unit(0.15, "lines"),
      plot.margin = margin(4,4,4,4),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text  = element_text(size = 7),
      legend.key.size = grid::unit(0.28, "cm"),
      legend.spacing.y = grid::unit(0.12, "cm")
    ) +
    guides(
      fill = guide_colorbar(
        title.position = "top",
        label.position = "right",
        barheight = grid::unit(3.0, "cm"),
        barwidth  = grid::unit(0.26, "cm")
      )
    ) +
    labs(
      title = paste0("Model pair differences per target"),
    ) +
    facet_wrap(~ trg, ncol = 4, nrow = 4)
}

# Per-dataset Self scores heatmap
model_self_scores_heatmap <- function(
    cdres,
    metric = "scc",
    model_order = c("deepttc","graphdrp","lgbm","xgboostdrp","uno"),
    datasets_cellline = if (exists("cellline")) cellline else NULL,
    datasets_exvivo   = if (exists("exvivo"))   exvivo   else NULL,
    show_labels = TRUE,
    limits = NULL,
    palette = "RdBu"
) {
  have_models <- intersect(model_order, unique(cdres$model))

  self_tbl <- cdres %>%
    filter(met == metric, src == trg, model %in% have_models) %>%
    group_by(model, trg) %>%
    summarise(mu = mean(value, na.rm=TRUE), .groups="drop")

  # Order: cell lines, ex vivo
  col_order <- c(datasets_cellline, datasets_exvivo)
  self_tbl$model <- factor(self_tbl$model, levels = rev(have_models))
  self_tbl$trg   <- factor(self_tbl$trg,   levels = col_order)

  if (is.null(limits)) {
    rng <- range(self_tbl$mu, na.rm=TRUE); pad <- diff(rng)*0.05
    limits <- c(rng[1]-pad, rng[2]+pad)
  }
  self_tbl$lab <- if (show_labels) sprintf("%.3f", self_tbl$mu) else ""

  ggplot(self_tbl, aes(x=trg, y=model, fill=mu)) +
    geom_tile(width=0.999, height=0.999, color="white", linewidth=0.25) +
    {if(show_labels) geom_text(aes(label=lab), size=3)} +
    scale_fill_distiller(palette=palette, direction=1,
                         limits=c(.65,.9), oob=scales::squish,
                         name="Spearman" )+
    coord_fixed() +
    theme_bw(11) +
    theme(axis.title=element_blank(),
          axis.text.x=element_text(angle=45,hjust=1),
          panel.grid=element_blank(),
          plot.title=element_text(hjust=0.5)) +
    labs(
      title="Performance within Datasets, Spearman",
      subtitle=""
    )
}

plot_drug_consistency <- function(method = "spearman",
                                  stat_y = c("IQR","SD"),
                                  label_top = 10,
                                  models = c("deepttc","graphdrp","lgbm","xgboostdrp","uno")) {
  stat_y <- match.arg(stat_y)

  d <- drug_transfer_scores(models = models, method = method)

  agg <- d %>%
    dplyr::group_by(drug) %>%
    dplyr::summarise(
      mu  = mean(perf, na.rm = TRUE),
      sd  = stats::sd(perf, na.rm = TRUE),
      iqr = stats::IQR(perf, na.rm = TRUE),
      n   = dplyr::n(),
      .groups = "drop"
    )
  agg$varstat <- if (stat_y == "IQR") agg$iqr else agg$sd

  # Label a few most interesting (highest mu and lowest varstat)
  labset <- dplyr::bind_rows(
    agg %>% dplyr::arrange(desc(mu)) %>% head(label_top),
    agg %>% dplyr::arrange(varstat)   %>% head(label_top)
  ) %>% dplyr::distinct(drug)

  p <- ggplot(agg, aes(x = mu, y = varstat)) +
    geom_point(alpha = 0.8) +
    ggrepel::geom_text_repel(data = dplyr::semi_join(agg, labset, by = "drug"),
                             aes(label = drug), size = 3, max.overlaps = 100) +
    theme_bw(11) +
    labs(title = "Drug transfer: mean vs variability across datasets/models",
         subtitle = paste0("Method: ", method, "; variability = ", stat_y),
         x = "Mean cross-dataset performance",
         y = paste0(stat_y, " of cross-dataset performance"))
}


# Build per-drug cross-dataset performance from predictions
# method = "spearman" (scc) or "pearson"
drug_transfer_scores <- function(models = c("deepttc","graphdrp","lgbm","xgboostdrp","uno"),
                                 exclude_src = c("beataml","mpnst"),
                                 targets_subset = c(cellline, exvivo),
                                 method = "spearman",
                                 min_n = 8) {
  stopifnot(method %in% c("spearman","pearson"))
  purrr::map_dfr(models, function(m) {
    df <- getModelPredictionData(dset = m) |>
      dplyr::mutate(model = m) |>
      dplyr::filter(!source %in% exclude_src,
                    target %in% targets_subset,
                    source != target)
    df <- .standardize_pred_cols(df)

    # per (src, trg, model, drug) correlation
    df %>%
      dplyr::group_by(source, target, model, drug) %>%
      dplyr::summarise(
        n = dplyr::n(),
        perf = suppressWarnings(stats::cor(y, yhat, method = method, use = "pairwise.complete.obs")),
        .groups = "drop"
      ) %>%
      dplyr::filter(is.finite(perf), n >= min_n)
  })
}
