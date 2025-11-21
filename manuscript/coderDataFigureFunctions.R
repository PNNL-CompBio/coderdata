# this is originally from "coderdataResultsFunctions.R
# Ridgeline plots

library(gridExtra)
library(patchwork)
library(ggplot2)
library(ggridges)
library(RColorBrewer)
library(plotly)

ridgelineMetricPlots <- function(metric, dataset = cdres, print_plots = TRUE, hide_y_grid = TRUE) {

  sr <- subset(dataset, met == metric)

  # Clamp only for r2
  if (identical(metric, "r2")) {
    sr$value <- ifelse(sr$value < -1, -1, sr$value)
  }

  # ---- order facets (unchanged look) ----
  mvals_src <- sr |>
    dplyr::group_by(src) |>
    dplyr::summarize(m = mean(value, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(m)
  sr$src <- factor(sr$src, levels = mvals_src$src)

  # common theme tweak to hide horizontal grid if desired
  grid_tweak <- if (hide_y_grid) ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank()
  ) else ggplot2::theme()

  # by source dataset (remove ridge outlines with color = NA)
  p1 <- sr |>
    ggplot2::ggplot(ggplot2::aes(x = value, y = trg, fill = model)) +
    ggridges::geom_density_ridges(alpha = 0.5, color = NA) +
    ggplot2::facet_grid(src ~ .) +
    ggplot2::ggtitle(paste0(metric, " by source dataset")) +
    ggplot2::scale_fill_manual(values = modelcolors) +
    grid_tweak

  # reorder targets
  mvals_trg <- sr |>
    dplyr::group_by(trg) |>
    dplyr::summarize(m = mean(value, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(m)
  sr$trg <- factor(sr$trg, levels = mvals_trg$trg)

  # by target dataset (remove ridge outlines with color = NA)
  p3 <- sr |>
    ggplot2::ggplot(ggplot2::aes(x = value, y = src, fill = model)) +
    ggridges::geom_density_ridges(alpha = 0.5, color = NA) +
    ggplot2::facet_grid(trg ~ .) +
    ggplot2::ggtitle(paste0(metric, " by target dataset")) +
    ggplot2::scale_fill_manual(values = modelcolors) +
    grid_tweak

  if (isTRUE(print_plots)) {
    print(p1)
    print(p3)
  }
  invisible(list(src = p1, trg = p3))
}


plot_transfer_heatmap <- function(mat_df, title, limits = c(-.2,.6)) {
  ggplot(mat_df, aes(x = trg, y = src, fill = mu)) +
    geom_tile(color = "white") +
    scale_fill_distiller(palette = "RdBu", limits = limits, direction = 1, na.value = "grey95") +
    geom_point(aes(size = n), shape = 21, stroke = .2, alpha = .9, fill = NA) +
    coord_fixed() + theme_bw(base_size = 11) +
    labs(title = title, x = "Target dataset", y = "Source dataset", fill = "Perf", size = "#pairs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


paired_grid_all_pretty <- function(cdres, metric = "scc",
                                   model_order = c("deepttc","graphdrp","lgbm","xgboostdrp","uno"),
                                   breaks = c(-0.25, 0, 0.25, 0.5, 0.75),
                                   ncol = 5) {
  # --- helpers ---
  safe_palette_for_targets <- function(trg_vec) {
    lvls <- sort(unique(trg_vec))
    stopifnot(length(lvls) > 0)
    if (exists("datasetcolors")) {
      base <- datasetcolors[names(datasetcolors) %in% lvls]
      missing <- setdiff(lvls, names(base))
      if (length(missing) > 0) {
        add <- setNames(scales::hue_pal()(length(missing)), missing)
        base <- c(base, add)
      }
      base[lvls]
    } else {
      setNames(scales::hue_pal()(length(lvls)), lvls)
    }
  }
  pval_fmt <- function(p) ifelse(is.na(p), "p=NA",
                                 ifelse(p < 1e-3, "p<0.001", sprintf("p=%.3g", p)))

  have_models <- intersect(model_order, unique(cdres$model))
  if (length(have_models) < 2) stop("Need at least two models present in cdres.")

  # 1) cell means per (model, src, trg)
  cell_means <- cdres %>%
    dplyr::filter(met == metric, src != trg, model %in% have_models) %>%
    dplyr::group_by(model, src, trg) %>%
    dplyr::summarise(mu = mean(value, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = model, values_from = mu)

  # 2) long table for every unordered model pair
  pairs <- t(combn(have_models, 2))
  long_list <- apply(pairs, 1, function(pr) {
    a <- pr[1]; b <- pr[2]
    if (!all(c(a, b) %in% names(cell_means))) return(NULL)
    cell_means %>%
      dplyr::select(src, trg, !!rlang::sym(a), !!rlang::sym(b)) %>%
      tidyr::drop_na(!!rlang::sym(a), !!rlang::sym(b)) %>%
      dplyr::transmute(src, trg, x = .data[[a]], y = .data[[b]],
                       pair = paste0(a, " vs ", b), a = a, b = b)
  })
  long_list <- long_list[!vapply(long_list, is.null, logical(1))]
  if (!length(long_list)) stop("No overlapping (src→trg) cells across model pairs.")
  long <- dplyr::bind_rows(long_list)

  # 3) stats per pair: winner by HL median (y − x), keep p for evidence
  stats <- long %>%
    dplyr::group_by(pair) %>%
    dplyr::summarise(
      a = dplyr::first(a),
      b = dplyr::first(b),
      wt = list(tryCatch(
        stats::wilcox.test(y, x, paired = TRUE, exact = FALSE,
                           conf.int = TRUE, conf.level = 0.95),
        error = function(e) NULL
      )),
      p  = ifelse(is.null(wt[[1]]), NA_real_, wt[[1]]$p.value),
      HL = ifelse(is.null(wt[[1]]), NA_real_, unname(wt[[1]]$estimate)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      winner = dplyr::case_when(
        is.na(HL)            ~ "NA",
        HL >  0              ~ b,   # y-axis model wins
        HL <  0              ~ a,   # x-axis model wins
        TRUE                 ~ "tie"
      ),
      label = paste0(pair, "\n", winner, ", ", pval_fmt(p))
    )
  lab_map <- setNames(stats$label, stats$pair)

  # 4) palette only for targets that appear
  pal <- safe_palette_for_targets(long$trg)

  # 5) common square limits + nice breaks
  rng <- range(c(long$x, long$y), na.rm = TRUE)
  if (!is.finite(diff(rng)) || diff(rng) == 0) rng <- rng + c(-0.5, 0.5)
  pad <- diff(rng) * 0.05
  lims <- c(rng[1] - pad, rng[2] + pad)
  brx  <- breaks[breaks >= lims[1] & breaks <= lims[2]]
  if (!length(brx)) brx <- scales::pretty_breaks(5)(lims)

  ggplot(long, aes(x = x, y = y, color = trg)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.5, color = "grey55") +
    geom_point(alpha = 0.9, size = 0.7) +
    scale_color_manual(values = pal, name = "Target dataset", drop = FALSE) +
    coord_equal(xlim = lims, ylim = lims, expand = FALSE) +
    scale_x_continuous(breaks = brx) +
    scale_y_continuous(breaks = brx) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey88"),
      axis.text  = element_text(size = 8),
      axis.title = element_text(size = 10),
      legend.position   = "right",
      legend.title      = element_text(size = 10),
      legend.text       = element_text(size = 9),
      legend.key.size   = unit(0.35, "cm"),
      legend.spacing.y  = unit(0.2, "cm"),
      strip.text = element_text(size = 7, face = "bold",  # smaller facet titles
                                margin = margin(2, 2, 2, 2)),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    labs(
      title = paste0("Paired Comparisons across model pairs"),
      x = "Model 1 (average spearman correlation)", y = "Model 2 (average spearman correlation)"
    ) +
    facet_wrap(~ pair, labeller = labeller(pair = lab_map), ncol = ncol)
}


# Symmetric triangular heatmap using MEAN differences (ROW − COLUMN):
# - Upper triangle  (row < col):  cell line   Δmean = mean( row − column )
# - Lower triangle  (row > col):  ex vivo     Δmean = mean( row − column )
# - Diagonal        (row = col):  within-model Δmean = mean(ex vivo) − mean(cell line)
paired_heatmap_triangular_mean <- function(cdres,
                                           metric = "scc",
                                           model_order = c("deepttc","graphdrp","lgbm","xgboostdrp","uno"),
                                           limits = c(-0.15, 0.15),
                                           show_labels = TRUE) {

  have_models <- intersect(model_order, unique(cdres$model))
  stopifnot(length(have_models) >= 2)

  dom_of_trg <- function(trg) ifelse(trg %in% exvivo, "ex vivo", "cell line")

  # Per-(model, src, trg) average performance + domain tag
  cell_means <- cdres %>%
    dplyr::filter(met == metric, src != trg, model %in% have_models) %>%
    dplyr::group_by(model, src, trg) %>%
    dplyr::summarise(mu = mean(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(domain = dom_of_trg(trg))

  # --- Pairwise (between models) MEAN Δ and paired t-test, per domain -----------
  # IMPORTANT CHANGE: d <- x - y (row − column), and t.test(x, y, paired=TRUE)
  compute_pairwise_domain_mean <- function(df_domain) {
    wide <- tidyr::pivot_wider(df_domain, names_from = model, values_from = mu)
    prs  <- t(combn(have_models, 2))
    rows <- apply(prs, 1, function(pr) {
      a <- pr[1]; b <- pr[2]  # 'a' will be the row model; 'b' the column model if (row=a, col=b)
      sub <- wide %>% dplyr::select(src, trg, !!rlang::sym(a), !!rlang::sym(b)) %>% tidyr::drop_na()
      if (nrow(sub) < 3) return(NULL)
      x <- sub[[a]]; y <- sub[[b]]
      d  <- x - y                              # <-- row − column
      md <- mean(d, na.rm = TRUE)
      tt <- tryCatch(stats::t.test(x, y, paired = TRUE), error = function(e) NULL)
      tibble::tibble(a = a, b = b,
                     MEAN_D = md,
                     p      = if (!is.null(tt)) tt$p.value else NA_real_,
                     n      = sum(stats::complete.cases(d)))
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (!length(rows)) return(NULL)
    dplyr::bind_rows(rows)
  }

  by_dom <- cell_means %>%
    dplyr::group_by(domain) %>%
    tidyr::nest() %>%
    dplyr::mutate(res = purrr::map(data, compute_pairwise_domain_mean)) %>%
    tidyr::unnest(res) %>%
    dplyr::select(domain, a, b, MEAN_D, p, n)

  # --- Diagonal (within-model) MEAN gap: mean(ex vivo) − mean(cell line) --------
  diag_tbl <- purrr::map_dfr(have_models, function(m) {
    xe <- cell_means %>% dplyr::filter(model == m, domain == "ex vivo")   %>% dplyr::pull(mu)
    yc <- cell_means %>% dplyr::filter(model == m, domain == "cell line") %>% dplyr::pull(mu)
    if (length(xe) < 3 || length(yc) < 3) {
      tibble::tibble(model = m, MEAN_diag = NA_real_, p_diag = NA_real_)
    } else {
      md <- mean(xe, na.rm = TRUE) - mean(yc, na.rm = TRUE)
      tt <- tryCatch(stats::t.test(xe, yc, paired = FALSE), error = function(e) NULL) # Welch
      tibble::tibble(model = m, MEAN_diag = md,
                     p_diag = if (!is.null(tt)) tt$p.value else NA_real_)
    }
  })

  # Grid of cells
  grid <- expand.grid(row_model = have_models, col_model = have_models, stringsAsFactors = FALSE)
  grid$row_model <- factor(grid$row_model, levels = have_models)
  grid$col_model <- factor(grid$col_model, levels = have_models)
  grid$i <- as.integer(grid$row_model); grid$j <- as.integer(grid$col_model)

  # Lookup now returns Δmean = (ROW − COLUMN) for requested domain
  lookup_mean <- function(dom, row, col) {
    m1 <- by_dom %>% dplyr::filter(domain == dom, a == row, b == col)
    if (nrow(m1) == 1) return(c(MD = m1$MEAN_D, p = m1$p, found = TRUE))      # already row−col
    m2 <- by_dom %>% dplyr::filter(domain == dom, a == col, b == row)
    if (nrow(m2) == 1) return(c(MD = -m2$MEAN_D, p = m2$p, found = TRUE))     # invert if stored reversed
    c(MD = NA_real_, p = NA_real_, found = FALSE)
  }

  vals <- purrr::pmap_dfr(list(grid$i, grid$j, as.character(grid$row_model), as.character(grid$col_model)),
                          function(i, j, rm, cm) {
                            if (i == j) {
                              d <- diag_tbl %>% dplyr::filter(model == rm)
                              tibble::tibble(row_model = rm, col_model = cm,
                                             MD = d$MEAN_diag, p = d$p_diag,
                                             domain_used = "ex vivo − cell line")
                            } else if (i < j) {
                              v <- lookup_mean("cell line", rm, cm)
                              tibble::tibble(row_model = rm, col_model = cm,
                                             MD = as.numeric(v["MD"]), p = as.numeric(v["p"]),
                                             domain_used = "cell line")
                            } else {
                              v <- lookup_mean("ex vivo", rm, cm)
                              tibble::tibble(row_model = rm, col_model = cm,
                                             MD = as.numeric(v["MD"]), p = as.numeric(v["p"]),
                                             domain_used = "ex vivo")
                            }
                          })

  star <- function(p) dplyr::case_when(
    is.na(p) ~ "",
    p < 1e-4 ~ "****",
    p < 1e-3 ~ "***",
    p < 1e-2 ~ "**",
    p < 0.05 ~ "*",
    TRUE     ~ ""
  )
  vals$lab <- if (show_labels) sprintf("%.3f%s", vals$MD, star(vals$p)) else ""

  ggplot(vals, aes(x = col_model, y = row_model, fill = MD)) +
    geom_tile(color = "white", linewidth = 0.3) +
    {if (show_labels) geom_text(aes(label = lab), size = 3)} +
    scale_fill_distiller(palette = "RdBu", direction = 1,
                         limits = limits, oob = scales::squish,
                         name = "Mean Difference\n(row − column)") +
    coord_fixed() +
    theme_bw(11) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = "Average Model Differences",
      subtitle = "upper = cell line, lower = ex vivo; diagonal = ex vivo − cell line")
}


plot_drug_leaderboard <- function(method = "spearman",
                                  top_k = 20,
                                  models = c("deepttc","graphdrp","lgbm","xgboostdrp","uno")) {
  d <- drug_transfer_scores(models = models, method = method)

  d2 <- d %>%
    dplyr::mutate(domain = ifelse(target %in% get0("exvivo", ifnotfound = character()),
                                  "ex vivo", "cell line")) %>%
    dplyr::group_by(drug, domain) %>%
    dplyr::summarise(mu = mean(perf, na.rm = TRUE),
                     se = stats::sd(perf, na.rm = TRUE) / sqrt(pmax(1, dplyr::n())),
                     n  = dplyr::n(),
                     .groups = "drop")

  # choose globally by mean across domains
  keep <- d2 %>%
    dplyr::group_by(drug) %>%
    dplyr::summarise(gmu = mean(mu, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(gmu))
  top <- head(keep$drug, top_k)
  bot <- tail(keep$drug, top_k)
  d3  <- d2 %>% dplyr::filter(drug %in% c(top, bot)) %>%
    dplyr::mutate(drug = factor(drug, levels = unique(c(rev(top), bot))))

  p <- ggplot(d3, aes(x = mu, y = drug, color = domain)) +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = mu - 1.96*se, xmax = mu + 1.96*se), height = 0) +
    scale_color_manual(values = c("cell line" = "#377eb8", "ex vivo" = "#e41a1c")) +
    theme_bw(11) +
    labs(title = "Drug transfer leaderboard (top & bottom)",
         subtitle = paste0("Method: ", method, " (mean ± 95% CI across models & sources)"),
         x = "Cross-dataset performance", y = NULL, color = "Target domain")
  print(p)
  invisible(p)
}


# Main heatmap with aggregated Self scores
model_by_dataset_heatmap <- function(
    cdres,
    metric = "scc",
    model_order = c("deepttc","graphdrp","lgbm","xgboostdrp","uno"),
    datasets_cellline = if (exists("cellline")) cellline else NULL,
    datasets_exvivo   = if (exists("exvivo"))   exvivo   else NULL,
    show_labels = TRUE,
    limits = NULL,
    palette = "RdBu",
    sep_lines = TRUE
) {
  stopifnot(all(c("model","met","src","trg","value") %in% names(cdres)))

  have_models <- intersect(model_order, unique(cdres$model))

  # Cross-prediction means (exclude self)
  base <- cdres %>%
    dplyr::filter(met == metric, src != trg, model %in% have_models) %>%
    dplyr::group_by(model, trg) %>%
    dplyr::summarise(mu = mean(value, na.rm = TRUE), .groups = "drop")

  # Self scores (src==trg)
  self_tbl <- cdres %>%
    dplyr::filter(met == metric, src == trg, model %in% have_models) %>%
    dplyr::group_by(model, trg) %>%
    dplyr::summarise(mu = mean(value, na.rm = TRUE), .groups = "drop")

  # Aggregated Self (cell line/ex vivo)
  self_cell <- self_tbl %>% filter(trg %in% datasets_cellline) %>%
    group_by(model) %>% summarise(trg = "Self (cell line)", mu = mean(mu, na.rm=TRUE), .groups="drop")
  self_ex   <- self_tbl %>% filter(trg %in% datasets_exvivo) %>%
    group_by(model) %>% summarise(trg = "Self (ex vivo)", mu = mean(mu, na.rm=TRUE), .groups="drop")

  # Aggregated averages (cross only)
  avg_cell <- base %>% filter(trg %in% datasets_cellline) %>%
    group_by(model) %>% summarise(trg = "Avg (cell line)", mu = mean(mu, na.rm=TRUE), .groups="drop")
  avg_ex <- base %>% filter(trg %in% datasets_exvivo) %>%
    group_by(model) %>% summarise(trg = "Avg (ex vivo)", mu = mean(mu, na.rm=TRUE), .groups="drop")
  avg_all <- base %>% filter(trg %in% c(datasets_cellline, datasets_exvivo)) %>%
    group_by(model) %>% summarise(trg = "Avg (overall)", mu = mean(mu, na.rm=TRUE), .groups="drop")

  # Column order
  col_order <- c(datasets_cellline, "Self (cell line)", "Avg (cell line)",
                 datasets_exvivo,   "Self (ex vivo)",  "Avg (ex vivo)",
                 "Avg (overall)")

  # Assemble long table
  long_tbl <- bind_rows(
    base %>% filter(trg %in% c(datasets_cellline, datasets_exvivo)),
    self_cell, self_ex, avg_cell, avg_ex, avg_all
  )
  long_tbl$model <- factor(long_tbl$model, levels = rev(have_models))
  long_tbl$trg   <- factor(long_tbl$trg, levels = col_order)

  # Color limits
  if (is.null(limits)) {
    rng <- range(long_tbl$mu, na.rm=TRUE); pad <- diff(rng)*0.05
    limits <- c(rng[1]-pad, rng[2]+pad)
  }
  long_tbl$lab <- if (show_labels) sprintf("%.3f", long_tbl$mu) else ""

  # Vertical separators
  vlines <- tibble::tibble(x = c(
    length(datasets_cellline) + 0.5,
    length(datasets_cellline) + 1 + 0.5,
    length(datasets_cellline) + 2 + length(datasets_exvivo) + 0.5,
    length(datasets_cellline) + 3 + length(datasets_exvivo) + 0.5
  ))

  ggplot(long_tbl, aes(x = trg, y = model, fill = mu)) +
    geom_tile(width=0.999, height=0.999, color="white", linewidth=0.25) +
    {if(show_labels) geom_text(aes(label=lab), size=3)} +
    scale_fill_distiller(palette=palette, direction=1,
                         limits=limits, oob=scales::squish,
                         name=paste0("Mean ", metric)) +
    {if(sep_lines) geom_vline(data=vlines, aes(xintercept=x), inherit.aes=FALSE,
                              color="grey30", linewidth=0.5)} +
    coord_fixed() +
    theme_bw(11) +
    theme(axis.title=element_blank(),
          axis.text.x=element_text(angle=45,hjust=1),
          panel.grid=element_blank(),
          plot.title=element_text(hjust=0.5)) +
    labs(
      title=paste0("Average performance by Model")
    )
}
