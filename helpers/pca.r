#' Prepare PCA Data from LFQ or LBQ Dataset
#'
#' This function prepares the PCA input data from either LFQ or LBQ quantitative proteomics data.
#' It reshapes the data, performs PCA, and merges sample-level annotation metadata.
#'
#' @param df A data frame containing the input proteomics data.
#' @param annotation A data frame with sample-level annotations. Must include all grouping columns.
#' @param type Character string, either `"LFQ"` or `"LBQ"`, indicating the data type.
#'
#' @return A tibble with PCA results merged with annotation metadata.
get_pca_data <- function(data, annotation, type = c("LFQ", "LBQ")) {
  type <- match.arg(type)
  
  if (type == "LFQ") {
    data_matrix <- data %>%
      dplyr::select(ID = originalRUN, Protein, Abundance = LogIntensities) %>%
      pivot_wider(names_from = Protein, values_from = Abundance) %>%
      column_to_rownames(var = "ID") %>%
      mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
    
    annotation_tmp <- annotation %>%
      select(ID = Run, Condition)
  } else if (type == "LBQ") {
    data_matrix <- data %>%
      mutate(ID = sprintf('%s_%s', Run, Channel)) %>%
      dplyr::select(ID, Protein, Abundance) %>%
      pivot_wider(names_from = Protein, values_from = Abundance) %>%
      column_to_rownames(var = "ID") %>%
      mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
    
    annotation_tmp <- annotation %>%
      mutate(ID = sprintf('Pool_%s_%s', TechRepMixture, Channel)) %>%
      select(ID, Channel, Condition, TechRepMixture) %>%
      mutate(across(everything(), ~ as.character(.)))
  }
  
  pca <- prcomp(data_matrix, center = TRUE, scale = TRUE)
  pca <- pca$x %>% as.data.frame() %>% rownames_to_column("ID")
  
  pca_data <- pca %>%
    left_join(annotation_tmp, by = "ID") %>%
    
    return(pca_data)
}

#' Generate PCA Plots Colored by Metadata
#'
#' Produces PCA scatter plots for visualizing sample clustering. Points are colored by different
#' metadata variables (e.g., Condition, BioReplicate, TechRepMixture).
#' The layout and content adapt based on whether the input is LFQ or LBQ.
#'
#' @param pca_data A tibble containing PCA output and sample-level metadata.
#' @param PCx Column name of the principal component for the x-axis (unquoted).
#' @param PCy Column name of the principal component for the y-axis (unquoted).
#' @param type Character string, either `"LFQ"` or `"LBQ"`. Determines which metadata variables
#'        are included in the plots.
#'
#' @return A patchwork-combined `ggplot` object with scatter plots colored by sample-level variables.
get_pca_plots <- function(pca_data, PCx, PCy, type = c("LFQ", "LBQ")) {
  type <- match.arg(type)
  
  # Base plot style
  base_plot <- function(data, fill_var, fill_label = NULL) {
    ggplot(data, aes(x = {{ PCx }}, y = {{ PCy }}, fill = .data[[fill_var]])) +
      geom_point(pch = 21, size = 3) +
      guides(fill = guide_legend(title.position = "top")) +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.justification.bottom = "top"
      ) +
      (if (!is.null(fill_label))
        labs(fill = fill_label)
       else
         NULL)
  }
  
  # Always plot Condition
  plots <- list(base_plot(pca_data, "Condition"))
  
  # Add TechRepMixture for LBQ
  if (type == "LBQ" && "TechRepMixture" %in% names(pca_data)) {
    plots <- append(plots, list(base_plot(
      pca_data, "TechRepMixture", "Tech. Replicate"
    )))
  }
  
  # Combine using patchwork
  plot_layout_widths <- rep(1, length(plots))
  combined_plot <- patchwork::wrap_plots(plots, widths = plot_layout_widths) +
    patchwork::plot_layout(axes = "collect")
  
  return(combined_plot)
}