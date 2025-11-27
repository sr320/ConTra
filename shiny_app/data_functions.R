# Data loading and processing functions for ConTra Shiny app

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(jsonlite)

# Load ConTra datasets
load_contra_data <- function() {
  # Prefer in-app data directory, fallback to parent project structure
  candidate_dirs <- c(
    file.path("data", "cleaned_datasets"),
    file.path("..", "data", "cleaned_datasets")
  )
  data_dir <- NULL
  for (cand in candidate_dirs) {
    if (file.exists(file.path(cand, "gene_counts_cleaned.csv"))) {
      data_dir <- cand
      break
    }
  }
  if (is.null(data_dir)) {
    stop("Could not find data/cleaned_datasets directory with required CSVs in app or parent.")
  }

  # Try to locate an interactions file
  interactions <- load_interactions_df()
  
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
  } else if (names(lncrna)[1] == "") {
    names(lncrna)[1] <- "entity"
  }
  
  if ("Geneid" %in% names(mirna)) {
    mirna <- mirna %>% rename(entity = Geneid)
  } else if ("...1" %in% names(mirna)) {
    mirna <- mirna %>% rename(entity = `...1`)
  } else if ("Name" %in% names(mirna)) {
    mirna <- mirna %>% rename(entity = Name)
  } else if (names(mirna)[1] == "") {
    names(mirna)[1] <- "entity"
  }
  
  if ("Geneid" %in% names(methyl)) {
    methyl <- methyl %>% rename(entity = Geneid)
  } else if ("...1" %in% names(methyl)) {
    methyl <- methyl %>% rename(entity = `...1`)
  } else if ("CpG" %in% names(methyl)) {
    methyl <- methyl %>% rename(entity = CpG)
  } else if (names(methyl)[1] == "") {
    names(methyl)[1] <- "entity"
  }
  
  list(
    genes = genes,
    lncrna = lncrna,
    mirna = mirna,
    methyl = methyl,
    interactions = interactions
  )
}

# Attempt to locate and load multi_way_interactions.csv from common locations
load_interactions_df <- function() {
  candidate_files <- c(
    file.path("data", "multi_way_interactions.csv"),
    file.path("data", "interactions", "multi_way_interactions.csv"),
    file.path("..", "data", "multi_way_interactions.csv")
  )
  found <- NULL
  for (f in candidate_files) {
    if (file.exists(f)) {
      found <- f
      break
    }
  }
  if (is.null(found)) {
    # Search most recent in ../output/**/tables/
    out_root <- file.path("..", "output")
    if (dir.exists(out_root)) {
      matches <- list.files(out_root, pattern = "multi_way_interactions\\.csv$", recursive = TRUE, full.names = TRUE)
      if (length(matches) > 0) {
        info <- file.info(matches)
        ord <- order(info$mtime, decreasing = TRUE)
        found <- matches[ord[1]]
      }
    }
  }
  if (is.null(found)) {
    message("Interactions file not found; falling back to sample regulators.")
    return(NULL)
  }
  message("Loading interactions from: ", found)
  df <- readr::read_csv(found, show_col_types = FALSE)
  df
}

# Parse the regulator_types field which may be JSON-like or double-encoded
parse_regulator_list <- function(txt) {
  if (is.null(txt) || is.na(txt) || nchar(trimws(txt)) == 0) return(character(0))
  val <- as.character(txt)
  # If wrapped in quotes with inner list, strip outer quotes
  if (grepl('^\"\\[', val) || grepl("^'\\[", val)) {
    val <- sub('^\"', '', val)
    val <- sub('\"$', '', val)
    val <- sub("^'", '', val)
    val <- sub("'$", '', val)
  }
  # Normalize single quotes to double quotes for JSON parsing
  val_json <- gsub("'", '"', val, fixed = TRUE)
  regs <- tryCatch(jsonlite::fromJSON(val_json), error = function(e) NULL)
  if (is.null(regs)) return(character(0))
  as.character(regs)
}

extract_core_gene_id <- function(gene_name) {
  # Extract core ID like FUN_005353 from strings like MOC2A_CHLRE_FUN_005353
  m <- regmatches(gene_name, regexpr("FUN_\\d+", gene_name))
  if (length(m) == 1 && nchar(m) > 0) {
    return(m)
  }
  return(gene_name)
}

# Get gene choices for dropdown (optionally filtered to those with interactions)
get_gene_choices <- function(genes_df, interactions_df = NULL) {
  if (nrow(genes_df) == 0 || !"entity" %in% names(genes_df)) {
    return("No genes available")
  }
  gene_ids <- genes_df$entity
  
  if (!is.null(interactions_df) && "gene" %in% names(interactions_df)) {
    interaction_gene_raw <- interactions_df$gene
    interaction_gene_core <- vapply(interaction_gene_raw, extract_core_gene_id, character(1))
    valid_set <- unique(c(interaction_gene_raw, interaction_gene_core))
    gene_ids <- intersect(gene_ids, valid_set)
  }
  
  if (length(gene_ids) == 0) {
    gene_ids <- genes_df$entity
  }
  
  head(sort(unique(gene_ids)), 50)
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
  
  # Determine regulators: prefer interactions file if available, else fallback to samples
  regulators <- character(0)
  if (!is.null(data_list$interactions) &&
      all(c("gene", "regulator_types") %in% names(data_list$interactions))) {
    # Match either exact gene or core ID match
    core_ids <- vapply(data_list$interactions$gene, extract_core_gene_id, character(1))
    match_idx <- which(data_list$interactions$gene == gene_id | core_ids == gene_id)
    if (length(match_idx) > 0) {
      row <- data_list$interactions[match_idx[1], ]
      regulators <- parse_regulator_list(row$regulator_types[[1]])
    }
  }

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
  
  if (length(regulators) > 0) {
    # Use real regulators from interactions
    for (r in regulators) {
      if (startsWith(r, "mirna_")) {
        mirna_id <- sub("^mirna_", "", r)
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
      } else if (startsWith(r, "lncrna_")) {
        lncrna_id <- sub("^lncrna_", "", r)
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
      } else if (startsWith(r, "methylation_")) {
        cpg_id <- sub("^methylation_", "", r)
        methyl_row <- methyl_tp[methyl_tp$entity == cpg_id, ]
        if (nrow(methyl_row) > 0) {
          for (tp in c("T1", "T2", "T3", "T4")) {
            if (tp %in% names(methyl_row) && !is.na(methyl_row[[tp]])) {
              records <- append(records, list(data.frame(
                entity = cpg_id,
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
  } else {
    # Fallback: sample first few from each type
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