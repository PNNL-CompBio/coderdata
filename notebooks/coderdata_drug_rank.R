# coderdata_drug_rank.R
# Rank drugs from CoderData transcriptomics or proteomics signatures.

suppressPackageStartupMessages({
  library(GSVA)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
})


read_signatures <- function(path) {
  if (grepl("\\.gmt$", path, ignore.case = TRUE)) {
    lines <- readLines(path)
    result <- lapply(lines, function(line) {
      fields <- strsplit(line, "\\t", fixed = FALSE)[[1]]
      fields[-c(1, 2)]
    })
    names(result) <- vapply(
      lines,
      function(line) strsplit(line, "\\t", fixed = FALSE)[[1]][1],
      character(1)
    )
    return(result)
  }

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Install jsonlite to read JSON signatures: install.packages('jsonlite')")
  }
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}


normalize_signatures <- function(
    signatures,
    available_features,
    min_signature_features
) {
  available <- unique(as.character(available_features))
  gene_sets <- list()
  mapping <- list()
  reports <- list()

  for (signature_name in names(signatures)) {
    definition <- signatures[[signature_name]]
    two_sided <- is.list(definition) &&
      any(c("up", "down") %in% names(definition))

    if (two_sided) {
      up_requested <- unique(as.character(definition$up %||% character()))
      down_requested <- unique(as.character(definition$down %||% character()))
    } else {
      up_requested <- unique(as.character(unlist(definition)))
      down_requested <- character()
    }

    up <- intersect(up_requested, available)
    down <- intersect(down_requested, available)

    reports[[length(reports) + 1]] <- tibble(
      signature = signature_name,
      component = c("up", "down"),
      requested = c(length(up_requested), length(down_requested)),
      matched = c(length(up), length(down)),
      missing_features = c(
        paste(setdiff(up_requested, available), collapse = ";"),
        paste(setdiff(down_requested, available), collapse = ";")
      )
    )

    if (length(up) < min_signature_features ||
        (two_sided && length(down) < min_signature_features)) {
      warning("Dropping signature '", signature_name,
              "' because too few features matched.")
      next
    }

    up_name <- paste0(signature_name, "__UP")
    gene_sets[[up_name]] <- up

    if (two_sided) {
      down_name <- paste0(signature_name, "__DOWN")
      gene_sets[[down_name]] <- down
    } else {
      down_name <- NA_character_
    }

    mapping[[signature_name]] <- list(
      up = up_name,
      down = down_name
    )
  }

  if (length(mapping) == 0) {
    stop("No signature met the feature-matching threshold.")
  }

  list(
    gene_sets = gene_sets,
    mapping = mapping,
    feature_match = bind_rows(reports)
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

score_signatures <- function(
    matrix,
    normalized,
    method,
    min_signature_features
) {
  if (method == "ssgsea") {
    parameters <- GSVA::ssgseaParam(
      exprData = matrix,
      geneSets = normalized$gene_sets,
      minSize = min_signature_features,
      maxSize = Inf,
      alpha = 0.25,
      normalize = TRUE,
      checkNA = "auto",
      use = "na.rm",
      verbose = FALSE
    )
  } else {
    parameters <- GSVA::gsvaParam(
      exprData = matrix,
      geneSets = normalized$gene_sets,
      minSize = min_signature_features,
      maxSize = Inf,
      kcdf = "Gaussian",
      checkNA = "auto",
      use = "na.rm",
      verbose = FALSE
    )
  }

  raw_scores <- GSVA::gsva(parameters, verbose = FALSE)

  map_dfr(names(normalized$mapping), function(signature_name) {
    keys <- normalized$mapping[[signature_name]]
    score <- raw_scores[keys$up, , drop = TRUE]
    if (!is.na(keys$down)) {
      score <- score - raw_scores[keys$down, , drop = TRUE]
    }
    tibble(
      sample_id = names(score),
      signature = signature_name,
      signature_score = as.numeric(score)
    )
  })
}

#' Rank drugs associated with transcriptomic or proteomic signatures
#'
#' @param omics Long-form CoderData transcriptomics or proteomics table.
#' @param experiments Long-form CoderData experiments table.
#' @param signatures Named list. Values can be feature vectors or
#'   list(up = ..., down = ...).
#' @param modality "transcriptomics" or "proteomics".
#' @param method "ssgsea" or "gsva".
#' @param association "spearman", "pearson", or "quartile".
#' @return List containing ranking, scores, feature_match, and settings.
rank_coderdata_drugs <- function(
    omics_matrix,
    experiment_data,
    signatures,
    modality = c("transcriptomics", "proteomics"),
    method = c("ssgsea", "gsva"),
    association = c("spearman", "pearson", "quartile"),
    sample_col = "improve_sample_id",
    #feature_col = "entrez_id",
    #value_col = NULL,
    drug_col = "improve_drug_id",
    #drug_name_col = NULL,
   # drug_data, 
    response_col = "fit_auc",
    lower_is_sensitive = TRUE,
    min_samples = 20L,
    min_signature_features = 3L,
    max_missing_fraction = NULL,
    adjust_scope = c("signature", "global")
) {
  modality <- match.arg(modality)
  method <- match.arg(method)
  association <- match.arg(association)
  adjust_scope <- match.arg(adjust_scope)

  value_col <- modality
  max_missing_fraction <- max_missing_fraction %||%
    if (modality == "proteomics") 0.75 else 0.50

#  omics_matrix <- prepare_omics(
#    omics,
#    sample_col,
#    feature_col,
#    value_col,
#    max_missing_fraction
#  )

  normalized <- normalize_signatures(
    signatures,
    rownames(omics_matrix),
    min_signature_features
  )

  scores <- score_signatures(
    omics_matrix,
    normalized,
    method,
    min_signature_features
  ) |>
    mutate(improve_sample_id = as.numeric(sample_id))

#  required <- c(sample_col, drug_col, response_col)
#  required <- required[!is.na(required) & nzchar(required)]
#  missing_columns <- setdiff(required, names(experiments))
#  if (length(missing_columns) > 0) {
#    stop("Missing experiment columns: ",
#         paste(missing_columns, collapse = ", "))
#  }

#  experiment_data <- experiments %>%
#    transmute(
#      sample_id = as.character(.data[[sample_col]]),
#      drug_id = as.character(.data[[drug_col]]),
#      drug_name = if (is.null(drug_name_col)) {
#        as.character(.data[[drug_col]])
#      } else {
#        as.character(.data[[drug_name_col]])
#      },
#      response = suppressWarnings(as.numeric(.data[[response_col]]))
#    ) %>%
#    filter(
#      !is.na(sample_id), sample_id != "",
#      !is.na(drug_id), drug_id != "",
#      is.finite(response)
#    ) %>%
#    group_by(sample_id, drug_id, drug_name) %>%
#    summarize(response = median(response, na.rm = TRUE), .groups = "drop")

  joined <- inner_join(scores, experiment_data, by = "improve_sample_id") #|>
 #   inner_join(drug_data, by = 'improve_drug_id')
  if (nrow(joined) == 0) {
    stop("No overlapping samples between omics and experiments.")
  }

  test_group <- function(data) {
    if (nrow(data) < min_samples ||
        n_distinct(data$signature_score) < 3 ||
        n_distinct(data[[response_col]]) < 3) {
      return(tibble())
    }

    if (association %in% c("spearman", "pearson")) {
      test <- suppressWarnings(cor.test(
        data$signature_score,
        data[[response_col]],
        method = association,
        exact = FALSE
      ))
      raw_effect <- unname(test$estimate)
      p_value <- test$p.value
      association_type <- "correlation"
      n_high <- NA_integer_
      n_low <- NA_integer_
    } else {
      q1 <- quantile(data$signature_score, 0.25, na.rm = TRUE)
      q3 <- quantile(data$signature_score, 0.75, na.rm = TRUE)
      low <- data[[response_col]][data$signature_score <= q1]
      high <- data[[response_col]][data$signature_score >= q3]
      if (length(low) < 3 || length(high) < 3) return(tibble())
      test <- suppressWarnings(wilcox.test(high, low, exact = FALSE))
      raw_effect <- mean(high) - mean(low)
      p_value <- test$p.value
      association_type <- "high_minus_low_response"
      n_high <- length(high)
      n_low <- length(low)
    }

    tibble(
      modality = modality,
      n = nrow(data),
      n_high = n_high,
      n_low = n_low,
      association_type = association_type,
      association = raw_effect,
      sensitivity_effect = if (lower_is_sensitive) -raw_effect else raw_effect,
      p_value = p_value,
      mean_response = mean(data[[response_col]])
    )
  }

  ranking <- joined %>%
    group_by(signature, improve_drug_id) %>%
    group_split() %>%
    map_dfr(function(group) {
      ids <- group %>%
        slice(1) %>%
        select(signature, improve_drug_id)
      result <- test_group(group)
      if (nrow(result) == 0) return(tibble())
      bind_cols(ids, result)
    })

  if (nrow(ranking) == 0) {
    stop("No drug-signature pair met the testing requirements.")
  }

  if (adjust_scope == "global") {
    ranking <- ranking %>%
      mutate(fdr = p.adjust(p_value, method = "BH"))
  } else {
    ranking <- ranking %>%
      group_by(signature) %>%
      mutate(fdr = p.adjust(p_value, method = "BH")) %>%
      ungroup()
  }

  ranking <- ranking %>%
    mutate(
      direction = if_else(
        sensitivity_effect > 0,
        "higher signature associated with sensitivity",
        "higher signature associated with resistance"
      )
    ) %>%
    group_by(signature) %>%
    arrange(fdr, desc(sensitivity_effect), .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    select(
      modality, signature, rank, improve_drug_id, n,
      association_type, association, sensitivity_effect,
      p_value, fdr, direction, everything()
    )

  structure(
    list(
      ranking = ranking,
      scores = scores,
      feature_match = normalized$feature_match,
      settings = list(
        modality = modality,
        method = method,
        association = association,
        value_col = value_col,
        response_col = response_col,
        lower_is_sensitive = lower_is_sensitive,
        max_missing_fraction = max_missing_fraction,
        adjust_scope = adjust_scope
      )
    ),
    class = "coderdata_drug_ranking"
  )
}

# Example:
# source("coderdata_drug_rank.R")
# result <- rank_coderdata_drugs(
#   omics = read_coderdata_table("ccle_proteomics.csv.gz"),
#   experiments = read_coderdata_table("prism_experiments.tsv.gz"),
#   signatures = read_signatures("protein_signatures.json"),
#   modality = "proteomics",
#   method = "ssgsea",
#   association = "spearman"
# )
# write_csv(result$ranking, "proteomics_drug_ranking.csv")
# write_csv(result$scores, "proteomics_signature_scores.csv")
# write_csv(result$feature_match, "proteomics_feature_matching.csv")
