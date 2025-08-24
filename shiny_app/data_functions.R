# Data loading and processing functions for ConTra Shiny app

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Load ConTra datasets
load_contra_data <- function(base_dir = "..") {
  data_dir <- file.path(base_dir, "data", "cleaned_datasets")
  
  message("Loading gene expression data...")
  genes <- read_csv(file.path(data_dir, "gene_counts_cleaned.csv"), 
                    show_col_types = FALSE)
  
  message("Loading lncRNA expression data...")
  lncrna <- read_csv(file.path(data_dir, "lncrna_counts_cleaned.csv"), 
                     show_col_types = FALSE)
  
  message("Loading miRNA expression data...")
  mirna <- read_csv(file.path(data_dir, "mirna_counts_cleaned.csv"), 
                    show_col_types = FALSE)
  
  message("Loading methylation data...")
  methyl <- read_csv(file.path(data_dir, "wgbs_counts_cleaned.csv"), 
                     show_col_types = FALSE)
  
  # Convert to proper format with entity names
  # Handle different file formats
  if ("...1" %in% names(genes)) {
    genes <- genes %>% rename(entity = `...1`)
  } else if (names(genes)[1] == "") {
    names(genes)[1] <- "entity"
  }
  
  if ("Geneid" %in% names(lncrna)) {
    lncrna <- lncrna %>% rename(entity = Geneid)
  } else if ("...1" %in% names(lncrna)) {
    lncrna <- lncrna %>% rename(entity = `...1`)
  }
  
  if ("Geneid" %in% names(mirna)) {
    mirna <- mirna %>% rename(entity = Geneid)
  } else if ("...1" %in% names(mirna)) {
    mirna <- mirna %>% rename(entity = `...1`)
  }
  
  if ("Geneid" %in% names(methyl)) {
    methyl <- methyl %>% rename(entity = Geneid)
  } else if ("...1" %in% names(methyl)) {
    methyl <- methyl %>% rename(entity = `...1`)
  }
  
  list(
    genes = genes,
    lncrna = lncrna,
    mirna = mirna,
    methyl = methyl
  )
}

# Get gene choices for dropdown
get_gene_choices <- function(genes_df) {
  if (nrow(genes_df) == 0 || !"entity" %in% names(genes_df)) {
    return("No genes available")
  }
  gene_ids <- genes_df$entity
  # Return first 50 genes for performance, sorted
  head(sort(gene_ids), 50)
}

# Aggregate timepoints (replicate the Python logic)
aggregate_timepoints <- function(df) {
  if (nrow(df) == 0 || !"entity" %in% names(df)) {
    return(data.frame())
  }
  
  tp_labels <- c("TP1", "TP2", "TP3", "TP4")
  
  result <- data.frame(entity = df$entity, stringsAsFactors = FALSE)
  
  for (tp in tp_labels) {
    # Find columns ending with this timepoint
    tp_cols <- names(df)[grepl(paste0(tp, "$"), names(df))]
    
    if (length(tp_cols) > 0) {
      # Calculate mean across replicates for this timepoint
      result[[paste0("T", substr(tp, 3, 3))]] <- rowMeans(df[tp_cols], na.rm = TRUE)
    }
  }
  
  result
}

# Get expression data for a gene and its regulators
get_gene_expression_data <- function(data_list, gene_id) {
  # Aggregate timepoints for all datasets
  genes_tp <- aggregate_timepoints(data_list$genes)
  lncrna_tp <- aggregate_timepoints(data_list$lncrna)
  mirna_tp <- aggregate_timepoints(data_list$mirna)
  methyl_tp <- aggregate_timepoints(data_list$methyl)
  
  # Check if datasets are valid
  if (nrow(genes_tp) == 0) {
    return(data.frame(entity = character(), type = character(), 
                      timepoint = character(), expression = numeric()))
  }
  
  # For this simplified version, we'll just get some sample regulators
  # In a full implementation, this would use interaction data
  
  records <- list()
  
  # Add gene data
  gene_row <- genes_tp[genes_tp$entity == gene_id, ]
  if (nrow(gene_row) > 0) {
    for (tp in c("T1", "T2", "T3", "T4")) {
      if (tp %in% names(gene_row) && !is.na(gene_row[[tp]])) {
        records <- append(records, list(data.frame(
          entity = gene_id,
          type = "gene",
          timepoint = tp,
          expression = as.numeric(gene_row[[tp]]),
          stringsAsFactors = FALSE
        )))
      }
    }
  }
  
  # Add some sample regulators (first few from each type)
  # miRNA regulators (first 3)
  if (nrow(mirna_tp) > 0) {
    sample_mirna <- head(mirna_tp$entity, 3)
    for (mirna_id in sample_mirna) {
      mirna_row <- mirna_tp[mirna_tp$entity == mirna_id, ]
      if (nrow(mirna_row) > 0) {
        for (tp in c("T1", "T2", "T3", "T4")) {
          if (tp %in% names(mirna_row) && !is.na(mirna_row[[tp]])) {
            records <- append(records, list(data.frame(
              entity = mirna_id,
              type = "miRNA",
              timepoint = tp,
              expression = as.numeric(mirna_row[[tp]]),
              stringsAsFactors = FALSE
            )))
          }
        }
      }
    }
  }
  
  # lncRNA regulators (first 3)
  if (nrow(lncrna_tp) > 0) {
    sample_lncrna <- head(lncrna_tp$entity, 3)
    for (lncrna_id in sample_lncrna) {
      lncrna_row <- lncrna_tp[lncrna_tp$entity == lncrna_id, ]
      if (nrow(lncrna_row) > 0) {
        for (tp in c("T1", "T2", "T3", "T4")) {
          if (tp %in% names(lncrna_row) && !is.na(lncrna_row[[tp]])) {
            records <- append(records, list(data.frame(
              entity = lncrna_id,
              type = "lncRNA",
              timepoint = tp,
              expression = as.numeric(lncrna_row[[tp]]),
              stringsAsFactors = FALSE
            )))
          }
        }
      }
    }
  }
  
  # Methylation regulators (first 3)
  if (nrow(methyl_tp) > 0) {
    sample_methyl <- head(methyl_tp$entity, 3)
    for (methyl_id in sample_methyl) {
      methyl_row <- methyl_tp[methyl_tp$entity == methyl_id, ]
      if (nrow(methyl_row) > 0) {
        for (tp in c("T1", "T2", "T3", "T4")) {
          if (tp %in% names(methyl_row) && !is.na(methyl_row[[tp]])) {
            records <- append(records, list(data.frame(
              entity = methyl_id,
              type = "methylation",
              timepoint = tp,
              expression = as.numeric(methyl_row[[tp]]),
              stringsAsFactors = FALSE
            )))
          }
        }
      }
    }
  }
  
  if (length(records) > 0) {
    result <- do.call(rbind, records)
    result$timepoint <- factor(result$timepoint, levels = c("T1", "T2", "T3", "T4"))
    return(result)
  } else {
    return(data.frame(entity = character(), type = character(), 
                      timepoint = character(), expression = numeric(),
                      stringsAsFactors = FALSE))
  }
}

# Calculate z-scores (replicate Python logic)
calculate_zscores <- function(expr_df) {
  if (nrow(expr_df) == 0) return(expr_df)
  
  # Calculate z-scores per entity
  expr_df %>%
    group_by(entity) %>%
    mutate(
      mean_expr = mean(expression, na.rm = TRUE),
      sd_expr = sd(expression, na.rm = TRUE),
      sd_expr = ifelse(sd_expr == 0 | is.na(sd_expr), 1, sd_expr),
      zscore = (expression - mean_expr) / sd_expr
    ) %>%
    ungroup() %>%
    select(-mean_expr, -sd_expr)
}

# Create facet z-score plot
create_facet_zscore_plot <- function(expr_df, gene_id, show_gene_overlay = TRUE) {
  if (nrow(expr_df) == 0) {
    return(ggplot() + 
           annotate("text", x = 1, y = 1, label = "No data available", size = 6) +
           theme_void())
  }
  
  # Calculate z-scores
  zdf <- calculate_zscores(expr_df)
  
  # Define colors matching Python version
  type_colors <- c(
    "miRNA" = "#1f77b4",
    "lncRNA" = "#2ca02c", 
    "methylation" = "#d62728",
    "gene" = "#000000"
  )
  
  # Get regulator types (exclude gene)
  regulator_types <- unique(zdf$type[zdf$type != "gene"])
  
  if (length(regulator_types) == 0) {
    return(ggplot() + 
           annotate("text", x = 1, y = 1, label = "No regulator data found", size = 6) +
           theme_void())
  }
  
  # Create the plot
  p <- ggplot()
  
  # Add regulator lines for each type
  for (rtype in regulator_types) {
    type_data <- zdf[zdf$type == rtype, ]
    if (nrow(type_data) > 0) {
      p <- p + 
        geom_line(data = type_data,
                  aes(x = timepoint, y = zscore, group = entity),
                  color = type_colors[[rtype]], 
                  alpha = 0.7, 
                  linewidth = 1) +
        geom_point(data = type_data,
                   aes(x = timepoint, y = zscore),
                   color = type_colors[[rtype]], 
                   alpha = 0.7, 
                   size = 2)
    }
  }
  
  # Add gene overlay if requested
  if (show_gene_overlay) {
    gene_data <- zdf[zdf$entity == gene_id & zdf$type == "gene", ]
    if (nrow(gene_data) > 0) {
      p <- p + 
        geom_line(data = gene_data,
                  aes(x = timepoint, y = zscore, group = entity),
                  color = type_colors[["gene"]], 
                  linewidth = 3) +
        geom_point(data = gene_data,
                   aes(x = timepoint, y = zscore),
                   color = type_colors[["gene"]], 
                   size = 3)
    }
  }
  
  # Facet by regulator type
  p <- p + 
    facet_wrap(~ type, scales = "free_y") +
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed", alpha = 0.6) +
    labs(
      title = paste("Gene", gene_id, "regulators (z-score by type)"),
      x = "Time Point",
      y = "Z-score"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      strip.text = element_text(size = 12, face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}