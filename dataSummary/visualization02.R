# Required libraries
library(circlize)
# library("tidyverse")
library(dplyr)
library(tools)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ComplexUpset)

# Function to read .csv or .csv.gz files
read_data <- function(file_path) {
  if (grepl("\\.gz$", file_path)) {
    return(read.csv(gzfile(file_path), stringsAsFactors = FALSE))
  } else {
    return(read.csv(file_path, stringsAsFactors = FALSE))
  }
}

read_and_preprocess_data <- function(file_names) {
  allowed_datasets <- c("transcriptomics", "proteomics", "mutations", "miRNA", "methylation", "copy_number")
  datasets <- names(file_names)
  
  if (!all(datasets %in% allowed_datasets)) {
    stop("Provided datasets are not all allowed.")
  }
  
  data_list <- lapply(file_names, function(f, dataset) {
    if (!file.exists(f)) {
      warning(paste("File", f, "does not exist. Creating empty dataset for", dataset))
      return(data.frame(entrez_id = integer(0), count = integer(0)))
    }
    data <- read.csv(f)
    data <- na.omit(data)
    data$X <- NULL
    return(data)
  }, datasets)
  
  range_list <- lapply(data_list, function(data) {
    if (nrow(data) == 0) {
      return(c(NA, NA))
    }
    range(data$entrez_id, na.rm = TRUE)
  })
  
  sector.data <- data.frame(
    sector = datasets,
    xlim1 = sapply(range_list, function(x) x[1]),
    xlim2 = sapply(range_list, function(x) x[2])
  )
  
  count_data <- lapply(data_list, function(data) {
    if (nrow(data) == 0) {
      return(data)
    }
    data %>% group_by(entrez_id) %>% summarise(count=n())
  })
  
  merged_data <- Reduce(function(...) merge(..., by="entrez_id", all=T), count_data)
  merged_data[is.na(merged_data)] <- 0
  merged_data <- as.data.frame(merged_data)
  
  names(merged_data)[-1] <- datasets
  
  list(data = merged_data, sector_data = sector.data)
}

# Circos plot function
generate_circos_plot <- function(processed_data,prefix) {
  background_color <- "#E0F2F1"
  save_name <- paste(prefix,"_circos",".png", sep = "")
  if (is.null(processed_data$data) || nrow(processed_data$data) == 0) {
    # This shouldn't trigger anymore.
    png(save_name, width = 1600, height = 1600, res = 600, bg = background_color)
    plot(1, type = "n", ann = FALSE)
    text(1, 1, "No data available", cex = 1.5)
    dev.off()
  } else {
    png(save_name, width = 1600, height = 1600, res = 600, bg = background_color)
    merged_data <- processed_data$data
    sector.data <- processed_data$sector_data
    datasets <- colnames(merged_data)[-1]
    
    # Colors
    dot_colors <- c("#fc8d62", "#8da0cb", "#e78ac3","#66c2a5", "#ffd92f","#a6d854")
    if (length(datasets) > length(dot_colors)) {
      warning("More datasets than colors provided. Colors will repeat.")
    }
    
    merged_data[,2:ncol(merged_data)] <- as.integer(merged_data[,2:ncol(merged_data)] != 0)
    merged_data[merged_data == 0] <- NA
    print(head(merged_data)) 
    circos.clear()
    
    # Compute unique Entrez IDs for each dataset
    unique_counts <- sapply(datasets, function(dataset) {
      length(unique(merged_data[merged_data[[dataset]] == 1, "entrez_id"]))
    })
    
    total_unique_ids = sum(unique_counts)
    
    # Number of datasets
    num_datasets <- length(datasets)
    
    # Each dot represents X Entrez IDs
    X = 200
    
    # Initialize the plotting with two tracks
    circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0))
    circos.initialize(factors = "a", xlim = cbind(c(0), c(total_unique_ids)))
    
    max_limit = (length(processed_data$data) - 1)  # Set the maximum limit
    
    # Adjust the data such that any value greater than the max_limit is set to max_limit
    df_adjusted <- merged_data
    df_adjusted[,2:ncol(merged_data)] <- pmin(df_adjusted[,2:ncol(merged_data)], max_limit)
    
    circos.trackPlotRegion(track.index = 1, bg.border = NA, ylim = c(0, max_limit), panel.fun = function(x, y) {
      for (i in 1:nrow(df_adjusted)) {
        current_position <- df_adjusted$entrez_id[i]
        bottom_limit <- 0
        
        # Loop through each dataset to create the segments of the stacked bar
        for(j in 2:ncol(df_adjusted)) {
          if (!is.na(df_adjusted[i, j])) { # Check for NA values
            top_limit <- bottom_limit + 1  # Each segment takes up only one unit of height
            circos.rect(xleft = current_position - 0.5, xright = current_position + 0.5,
                        ybottom = bottom_limit, ytop = top_limit,
                        col = dot_colors[j-1], border = dot_colors[j-1])
            bottom_limit <- top_limit
          } else {
            bottom_limit <- bottom_limit + 1  # Increment bottom_limit even if data is NA
          }
        }
      }
    })

    # Dots for the outer track
    circos.track(track.index = 2, ylim = c(1, num_datasets+1), bg.border = NA, panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim")
      circos.segments(rep(xlim[1], num_datasets), 1:num_datasets,
                      rep(xlim[2], num_datasets), 1:num_datasets,
                      col = "#CCCCCC")
      
      # Plotting dots
      for(i in 1:num_datasets) {
        num_dots = ceiling(unique_counts[i] / X)
        dot_positions = seq(0, unique_counts[i], length.out = num_dots)
        circos.points(dot_positions, rep(i, length(dot_positions)),
                      pch = 16, cex = .8, col = dot_colors[i])
      }
    })
    
    # Add the legend
    legend("center", legend = datasets, fill = dot_colors, cex = .6)
    dev.off()
  }
}

# Function to generate individual summary statistics
generate_individ_summary_stats <- function(file_names) {
  summary_df <- data.frame(
    Dataset = character(),
    Unique_Entrez_IDs = integer(),
    Unique_Sample_IDs = integer(),
    stringsAsFactors = FALSE
  )
  
  for(dataset_name in names(file_names)) {
    file_path <- file_names[[dataset_name]]
    if (file.exists(file_path)) {
      dataset <- read_data(file_path)
      unique_genes <- length(unique(dataset$entrez_id))
      unique_samples <- length(unique(dataset$improve_sample_id))
      
      summary_df <- rbind(summary_df, data.frame(
        Dataset = tools::toTitleCase(dataset_name),
        Unique_Entrez_IDs = unique_genes,
        Unique_Sample_IDs = unique_samples,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(summary_df)
}

# Group summary statistics function
generate_group_summary_stats <- function(file_names) {
  # Initialize an empty dataframe for summary
  summary_df <- data.frame(
    Dataset = character(),
    Samples = integer(),
    Variations = integer(),
    `Cancer Type` = character(),
    Source = character(),
    stringsAsFactors = FALSE
  )
  
  # Iterate through each dataset
  for(dataset_name in names(file_names)) {
    file_path <- file_names[[dataset_name]]
    
    # Check if file exists before reading
    if (file.exists(file_path)) {
      data <- read.csv(file_path)
      
      # Calculate unique counts
      unique_samples <- length(unique(data$improve_sample_id))
      unique_cancer_types <- length(unique(data$cancer_type))
      cancer_type_str <- if (unique_cancer_types == 1) {
        as.character(unique(data$cancer_type))
      } else {
        paste(unique_cancer_types, "cancer types")
      }
      
      unique_sources <- length(unique(data$other_id_source))
      source_str <- if(unique_sources == 1) {
        as.character(unique(data$other_id_source))
      } else if(unique_sources == 0) {
        unique_sources <- length(unique(data$id_source))
        paste(unique_sources, "sources")
      } else {
        paste(unique_sources, "sources")
      }
      
      # Append to the summary dataframe
      summary_df <- rbind(summary_df, data.frame(
        Dataset = dataset_name,
        Samples = unique_samples,
        `Cancer Type` = cancer_type_str,
        Source = source_str,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Return the summary dataframe
  return(summary_df)
}


# Function to generate group summary plot
generate_group_summary_plot <- function(all_file_names) {
  
  # Initialize a dataframe to store sample counts
  samples_df <- data.frame(
    DataType = character(),
    Samples = integer(),
    Source = character(),
    stringsAsFactors = FALSE
  )
  
  # Loop over the sources and data types to calculate unique sample counts
  for(source in names(all_file_names)) {
    for(data_type in names(all_file_names[[source]])) {
      file_path <- all_file_names[[source]][[data_type]]
      if (file.exists(file_path)) {
        data <- read.csv(file_path)
        unique_samples <- length(unique(data$improve_sample_id))
        samples_df <- rbind(samples_df, data.frame(
          DataType = data_type,
          Samples = unique_samples,
          Source = source,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Check if samples_df is empty
  if (nrow(samples_df) == 0) {
    print("No data available for plotting")
    return()
  }
  
  # Plot the data
  p <- ggplot(samples_df, aes(x = DataType, y = Samples, fill = Source)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Number of Samples by Data Type and Source",
         x = "Data Type",
         y = "Number of Samples") +
    scale_fill_manual(values = c("beataml" = "#fc8d62", "hcmi" = "#8da0cb", "depmap" = "#66c2a5", "cptac" = "#8511c1")) +
    theme(plot.background = element_rect(fill = background_color, color = background_color))
  ggsave('Fig5_Sample_Summary.png', p, height = 9, width = 12, bg = background_color)
}

# Data file names for each group
beataml_names <- list(
  transcriptomics = "beataml_transcriptomics.csv.gz",
  proteomics = "beataml_proteomics.csv.gz",
  mutations = "beataml_mutations.csv.gz"
)
hcmi_names <- list(
  transcriptomics = "hcmi_transcriptomics.csv.gz",
  mutations = "hcmi_mutations.csv.gz",
  copy_number = "hcmi_copy_number.csv.gz"
)
depmap_names <- list(
  transcriptomics = "depmap_transcriptomics.csv.gz",
  proteomics = "depmap_proteomics.csv.gz",
  miRNA = "depmap_miRNA.csv.gz",
  copy_number = "depmap_copy_number.csv.gz",
  mutations = "depmap_mutations.csv.gz"
)
cptac_names <- list(
  transcriptomics = "cptac_transcriptomics.csv.gz",
  proteomics = "cptac_proteomics.csv.gz",
  copy_number = "cptac_copy_number.csv.gz",
  mutations = "cptac_mutations.csv.gz"
)

# Combine all file names into one list
all_file_names <- list(
  beataml = beataml_names,
  hcmi = hcmi_names,
  depmap = depmap_names,
  cptac = cptac_names
)

save_summary <- function(file_group, individ_summary) {
  file_name <- paste0(file_group,"_table", ".csv")
  write.csv(individ_summary, file = file_name, row.names = FALSE)
}

# # Run individual summary stats for each data group
for (file_group_name in names(all_file_names)) {
  file_group <- all_file_names[[file_group_name]]
  individ_summary <- generate_individ_summary_stats(file_group)
  save_summary(file_group_name,individ_summary)
  processed_data <- read_and_preprocess_data(file_group)
  if (is.null(processed_data$data) || nrow(processed_data$data) != 0) {
    generate_circos_plot(processed_data,file_group_name)
  }
}


# Sample file names
samples_names <- list(
  HCMI = "hcmi_samples.csv",
  BEATAML = "beataml_samples.csv",
  DepMap = "depmap_samples.csv",
  CPTAC = "cptac_samples.csv"
)

# # Generate and print the summary
# group_summary <- generate_group_summary_stats(samples_names)

# Generate group summary plot
generate_group_summary_plot(all_file_names)
