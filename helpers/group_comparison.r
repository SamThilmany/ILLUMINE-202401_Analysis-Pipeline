#' Create Volcano Plots for Differentially Abundant Proteins
#'
#' Generates faceted volcano plots for each comparison in the dataset,
#' highlighting significantly differentially abundant proteins.
#'
#' @param data A data frame containing at least `log2FC`, `adj.pvalue`, and `LabelFactor` columns.
#' @param fdr_cutoff A numeric value indicating the FDR (adjusted p-value) cutoff for significance.
#' @param source A string label indicating the data source ("LFQ" or "LBQ").
#'
#' @return A `ggplot` object displaying volcano plots faceted by `LabelFactor`.
get_volcano_plots <- function(data, fdr_cutoff, source = c("LFQ", "LBQ")) {
  source <- match.arg(source)
  max_log2FC <- max(abs(data$log2FC), na.rm = TRUE)
  
  data <- data %>%
    group_by(LabelFactor) %>%
    mutate(Mean = mean(log2FC), SD = sd(log2FC)) %>%
    ungroup() %>%
    mutate(
      Source = source,
      Regulation = case_when(
        adj.pvalue < fdr_cutoff & (log2FC < (Mean - SD)) ~ "downregulated",
        adj.pvalue < fdr_cutoff &
          (log2FC > (Mean + SD)) ~ "upregulated",
        TRUE ~ "not significant"
      )
    )
  
  ggplot(data, aes(x = log2FC, y = -log10(adj.pvalue))) +
    facet_grid(Source ~ LabelFactor,
               labeller = labeller(LabelFactor = conditionLabels)) +
    scale_x_continuous(limits = c(-max_log2FC, max_log2FC),
                       expand = c(0.1, 0.1)) +
    geom_point(aes(color = Regulation), alpha = 0.1) +
    geom_hline(
      yintercept = -log10(fdr_cutoff),
      linetype = "dashed",
      color = hue_pal()(4)[1],
      linewidth = 0.25
    ) +
    geom_vline(
      aes(xintercept = (Mean - SD)),
      linetype = "dashed",
      color = hue_pal()(4)[1],
      linewidth = 0.25
    ) +
    geom_vline(
      aes(xintercept = (Mean + SD)),
      linetype = "dashed",
      color = hue_pal()(4)[1],
      linewidth = 0.25
    ) +
    scale_color_manual(
      values = c(
        "downregulated" = hue_pal()(4)[3],
        "upregulated" = hue_pal()(4)[1],
        "not significant" = "gray"
      )
    ) +
    labs(x = expression("log"[2] * " FC"),
         y = expression("-log"[10] * " (adj. p-value)")) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 6
      )
    )
}

#' Stacked volcano plots for LFQ and LBQ
#'
#' Creates a two‑panel figure: one panel for LFQ, one for LBQ.
#' Internally calls \code{get_volcano_plots()} for each data set and
#' stacks the results with \pkg{patchwork}.
#'
#' @param lfq_data A data frame containing at least `log2FC`, `adj.pvalue`, and `LabelFactor` columns.
#' @param lfq_fdr_cutoff A numeric value indicating the FDR (adjusted p-value) cutoff for significance.
#' @param lbq_data A data frame containing at least `log2FC`, `adj.pvalue`, and `LabelFactor` columns.
#' @param lbq_fdr_cutoff A numeric value indicating the FDR (adjusted p-value) cutoff for significance.
#'
#' @return A patchwork-combined `ggplot` object combining the two \code{ggplot}s.
get_combined_volcano_plot <- function(lfq_data,
                                      lfq_fdr_cutoff,
                                      lbq_data,
                                      lbq_fdr_cutoff) {
  lfq_plot <- get_volcano_plots(lfq_data, lfq_fdr_cutoff, "LFQ")
  lbq_plot <- get_volcano_plots(lbq_data, lbq_fdr_cutoff, "LBQ")
  
  (lfq_plot / lbq_plot) +
    plot_layout(heights = c(1, 1), axes = "collect")
}

#' Summarize Significant Protein Counts Across Conditions
#'
#' Calculates the number and percentage of significantly regulated proteins per condition,
#' and complements with the number of non-significant proteins.
#'
#' @param data A data frame containing at least `Protein`, `adj.pvalue`, and `LabelFactor` columns.
#' @param fdr_cutoff A numeric FDR threshold used to define significance.
#'
#' @return A data frame with columns `LabelFactor`, `Protein_Count`, `Percentage`,
#'         and `Significance`.
#'
#' @details Ensures all proteins are accounted for by computing both significant and
#' non-significant protein counts per `LabelFactor`.
get_significant_proteins_amount_data <- function(data, fdr_cutoff) {
  all <- data %>%
    group_by(LabelFactor) %>%
    summarise(Protein_Count = n_distinct(Protein))
  
  significant <- data %>%
    dplyr::filter(adj.pvalue < fdr_cutoff) %>%
    group_by(LabelFactor) %>%
    summarise(Protein_Count = n_distinct(Protein)) %>%
    left_join(all,
              by = "LabelFactor",
              suffix = c("_significant", "_all")) %>%
    mutate(Percentage = if_else(
      is.na(Protein_Count_significant),
      0,
      Protein_Count_significant
    ) / Protein_Count_all * 100) %>%
    select(LabelFactor, Protein_Count = Protein_Count_significant, Percentage)
  
  non.significant <- all %>%
    left_join(significant,
              by = "LabelFactor",
              suffix = c("_all", "_significant")) %>%
    mutate(
      Protein_Count = Protein_Count_all - if_else(
        is.na(Protein_Count_significant),
        0,
        Protein_Count_significant
      ),
      Percentage = (
        1 - if_else(
          is.na(Protein_Count_significant),
          0,
          Protein_Count_significant
        ) / Protein_Count_all
      ) * 100
    ) %>%
    select(LabelFactor, Protein_Count, Percentage)
  
  combined <- bind_rows(
    significant %>% mutate(Significance = "Yes"),
    non.significant %>% mutate(Significance = "No")
  )
}

#' Plot Significant vs Non-Significant Protein Counts as Bar Charts
#'
#' Creates a faceted bar chart showing the number of significant and non-significant proteins
#' per condition, annotated with counts and percentages.
#'
#' @param significant_proteins_amount_data A data frame as returned by `get_significant_proteins_amount_data()`,
#'        containing columns `LabelFactor`, `Protein_Count`, `Percentage`, and `Significance`.
#' @param source A string label indicating the data source ("LFQ" or "LBQ").
#'
#' @return A `ggplot` object with faceted bar charts per `LabelFactor`.
#'
#' @details Bars are annotated with both absolute counts and percentages.
#' Faceting is used to compare across multiple conditions.
get_significance_barchart <- function(significant_proteins_amount_data,
                                      source = c("LFQ", "LBQ")) {
  source <- match.arg(source)
  
  annotation <- significant_proteins_amount_data %>%
    transmute(
      LabelFactor,
      Protein_Count,
      Significance,
      nudge_y = if_else(
        Protein_Count > 0.8 * max(Protein_Count),
        -0.05 * max(Protein_Count),
        0.05 * max(Protein_Count)
      ),
      vjust = if_else(Protein_Count > 0.8 * max(Protein_Count), 1, 0),
      label1 = Protein_Count,
      label2 = round(Percentage, 2),
      Source = source
    ) %>%
    distinct()
  
  ggplot(
    significant_proteins_amount_data,
    aes(x = Significance, y = Protein_Count, fill = Significance)
  ) +
    geom_col(position = "dodge") +
    facet_grid(Source ~ LabelFactor,
               labeller = labeller(LabelFactor = conditionLabels)) +
    geom_text(
      data = annotation,
      aes(
        y = Protein_Count + nudge_y,
        vjust = vjust,
        label = sprintf('%s\n(%s %%)', comma(label1), label2)
      ),
      hjust = 0.5,
      size = 5 / .pt
    ) +
    labs(x = "Significant", y = "Protein Count") +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

#' Stacked Significant vs Non-Significant Bar Charts for LFQ and LBQ
#'
#' Creates a two‑panel figure: one panel for LFQ, one for LBQ.
#' Internally calls \code{get_significance_barchart()} for each data set and
#' stacks the results with \pkg{patchwork}.
#'
#' @param lfq_data A data frame containing at least `log2FC`, `adj.pvalue`, and `LabelFactor` columns.
#' @param lbq_data A data frame containing at least `log2FC`, `adj.pvalue`, and `LabelFactor` columns.
#'
#' @return A patchwork-combined `ggplot` object combining the two \code{ggplot}s.
get_combined_significance_barchart <- function(lfq_data, lbq_data) {
  lfq_plot <- get_significance_barchart(lfq_data, "LFQ")
  lbq_plot <- get_significance_barchart(lbq_data, "LBQ")
  
  (lfq_plot / lbq_plot) +
    plot_layout(heights = c(1, 1), axes = "collect")
}

#' Calculate fold change statistics per condition
#'
#' Computes the mean and standard deviation of log2 fold change (log2FC) values
#' for each experimental condition (LabelFactor). Prepares additional columns
#' for annotation in downstream plots.
#'
#' @param data A data.frame containing columns `LabelFactor` and `log2FC`.
#'
#' @return A data.frame summarizing mean, standard deviation, and annotation positions for each condition.
get_fold_change_stats <- function(data) {
  data %>%
    group_by(LabelFactor) %>%
    summarize(Mean = mean(log2FC), SD = sd(log2FC)) %>%
    mutate(
      LabelID = seq_along(LabelFactor),
      xpos = rep(-Inf, nrow(.)),
      ypos = rep(Inf, nrow(.)),
      hjustvar = rep(1.1, nrow(.)),
      vjustvar = rep(2, nrow(.)),
      annotation = sub("-", "\u2212", sprintf("Mean: %.2f", Mean))
    )
}

#' Plot fold change density plots with highlighted tails
#'
#' Creates density plots of log2 fold changes for each experimental condition,
#' with shaded areas indicating extreme fold changes (more than ±1 standard deviation from the mean).
#'
#' @param data A data.frame containing log2FC values and corresponding LabelFactor annotations.
#' @param stats A data.frame output from `get_fold_change_stats()`, providing means and standard deviations.
#' @param source A string label indicating the data source ("LFQ" or "LBQ").
#'
#' @return A `ggplot` object showing density curves and shaded extreme regions per condition.
get_fold_change_density_plot <- function(data, stats, source = c("LFQ", "LBQ")) {
  source <- match.arg(source)
  
  data <- data %>% mutate(Source = source)
  
  max_log2FC <- 5 * sd(data$log2FC) # Only plot the values that are in between 5 standard deviations
  
  fold_change_density_base_plot <- ggplot(data, aes(x = log2FC)) +
    facet_grid(Source ~ LabelFactor,
               labeller = labeller(LabelFactor = conditionLabels)) +
    scale_x_continuous(limits = c(-max_log2FC, max_log2FC),
                       expand = c(0.1, 0.1)) +
    geom_density() +
    geom_vline(
      data = stats,
      aes(xintercept = Mean - SD),
      color = hue_pal()(4)[1],
      linetype = "dashed"
    ) +
    geom_vline(
      data = stats,
      aes(xintercept = Mean + SD),
      color = hue_pal()(4)[1],
      linetype = "dashed"
    ) +
    geom_text(
      data = stats,
      aes(
        x = xpos,
        y = ypos,
        hjust = hjustvar,
        vjust = vjustvar,
        label = annotation
      ),
      angle = 90
    ) +
    labs(x = expression("log"[2] * " FC"), y = "Density") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 6
      )
    )
  
  fold_change_density_base_plot_data = ggplot_build(fold_change_density_base_plot)$data[[1]] %>%
    mutate(PanelInt = as.numeric(PANEL)) %>%
    left_join(stats, by = c("PanelInt" = "LabelID"))
  
  fold_change_density_plot <- fold_change_density_base_plot +
    geom_area(
      data = subset(fold_change_density_base_plot_data, x < Mean - SD),
      aes(x = x, y = density),
      fill = hue_pal()(4)[1],
      alpha = 0.5
    ) +
    geom_area(
      data = subset(fold_change_density_base_plot_data, x > Mean + SD),
      aes(x = x, y = density),
      fill = hue_pal()(4)[1],
      alpha = 0.5
    )
}

#' Stacked Density Plots for LFQ and LBQ
#'
#' Creates a two‑panel figure: one panel for LFQ, one for LBQ.
#' Internally calls \code{get_fold_change_density_plot()} for each data set and
#' stacks the results with \pkg{patchwork}.
#'
#' @param lfq_data A data.frame containing log2FC values and corresponding LabelFactor annotations.
#' @param lfq_stats A data.frame output from `get_fold_change_stats()`, providing means and standard deviations.
#' @param lbq_data A data.frame containing log2FC values and corresponding LabelFactor annotations.
#' @param lbq_stats A data.frame output from `get_fold_change_stats()`, providing means and standard deviations.
#'
#' @return A patchwork-combined `ggplot` object combining the two \code{ggplot}s.
get_combined_fold_change_density_plot <- function(lfq_data, lfq_stats, lbq_data, lbq_stats) {
  lfq_plot <- get_fold_change_density_plot(lfq_data, lfq_stats, "LFQ")
  lbq_plot <- get_fold_change_density_plot(lbq_data, lbq_stats, "LBQ")
  
  (lfq_plot / lbq_plot) +
    plot_layout(heights = c(1, 1), axes = "collect")
}

#' Filter data based on log2 fold change and adjusted p-value thresholds
#'
#' This function filters a dataset to retain only proteins with significant
#' log2 fold changes (outside ±1 standard deviation from the group mean) and an
#' adjusted p-value below a specified threshold. It also prints a summary showing
#' how many entries were removed per condition.
#'
#' @param data A data.frame containing the full dataset with at least `LabelFactor`, `log2FC`,
#'        and `adj.pvalue` columns.
#' @param stats A data.frame returned by `get_fold_change_stats()`, containing `Mean`, `SD`,
#'        and `LabelFactor` for each condition.
#' @param fdr_cutoff Numeric. The adjusted p-value threshold (e.g., 0.05) used to filter for
#'        statistical significance.
#'
#' @return A filtered data.frame containing only the proteins with significant fold changes and
#'         adjusted p-values. Also displays a summary table of filtering impact per condition.
get_log2fc_and_pvalue_filtered_data <- function(data, stats, fdr_cutoff) {
  current_summary <- data %>%
    group_by(LabelFactor) %>%
    summarize(Before = n())
  
  filtered_pvalue <- data %>%
    left_join(stats, by = "LabelFactor") %>%
    dplyr::filter(adj.pvalue < fdr_cutoff)
  
  after_pvalue_summary <- filtered_pvalue %>%
    group_by(LabelFactor) %>%
    summarize(After.pvalue = n())
  
  summary <- current_summary %>%
    left_join(after_pvalue_summary, by = "LabelFactor") %>%
    replace_na(list(Before = 0, After.pvalue = 0)) %>%
    mutate(Removed.pvalue = Before - After.pvalue)
  
  filtered_pvalue_log2fc <- filtered_pvalue %>%
    dplyr::filter(log2FC < (Mean - SD) | log2FC > (Mean + SD))
  
  after_pvalue_log2fc_summary <- filtered_pvalue_log2fc %>%
    group_by(LabelFactor) %>%
    summarize(After.pvalue.log2fc = n())
  
  summary <- summary %>%
    left_join(after_pvalue_log2fc_summary, by = "LabelFactor") %>%
    replace_na(list(After.pvalue.log2fc = 0)) %>%
    mutate(
      Removed.pvalue.log2FC = After.pvalue - After.pvalue.log2fc,
      Removed = Before - After.pvalue.log2fc
    )
  
  display(summary)
  
  return(filtered_pvalue_log2fc)
}

#' Create a Venn Diagram for Two Protein Sets
#'
#' Generates a Venn diagram comparing proteins between LFQ and LBQ datasets,
#' annotated with a custom subtitle.
#'
#' @param lfq_data A character vector of protein identifiers from the LFQ dataset.
#' @param lbq_data A character vector of protein identifiers from the LBQ dataset.
#' @param name A string to be used as the subtitle of the Venn diagram.
#'
#' @return A `ggplot` object representing the Venn diagram.
#'
#' @details
#' - The Venn diagram shows the overlap between the two datasets.
get_venn_diagram_plot_label <- function(lfq_data, lbq_data, name) {
  lfq_data <- as.character(lfq_data)
  lbq_data <- as.character(lbq_data)
  
  data_list <- list("LFQ" = lfq_data, "LBQ" = lbq_data)
  
  venn = Venn(data_list)
  data = process_data(venn, shape_id = "201")
  
  plot_venn(
    data,
    set_size = 2,
    label_size = 2,
    edge_size = 0.5
  ) +
    labs(subtitle = name) +
    theme(
      legend.position = "none",
      text = element_text(size = default_font_size, color = default_text_color),
      plot.title.position = "panel",
      plot.subtitle = element_text(
        size = 0.8 * default_font_size,
        color = default_text_color,
        hjust = 0.5,
        margin = margin(
          t = 0,
          r = 0,
          b = 5,
          l = 0
        )
      ),
      plot.margin = margin(
        t = 5,
        r = 0,
        b = 5,
        l = 0
      )
    )
}

#' Create Multiple Venn Diagrams for Protein Overlap
#'
#' Generates a collection of Venn diagrams for all proteins and per-condition protein sets
#' across LFQ and LBQ datasets.
#'
#' @param lfq_summarizationData A data frame containing `Protein` and `LabelFactor` columns for the LFQ dataset.
#' @param lbq_summarizationData A data frame containing `Protein` and `LabelFactor` columns for the LBQ dataset.
#' @param lfq_comparisonData A data frame containing `Protein` and `LabelFactor` columns for the LFQ dataset.
#' @param lbq_comparisonData A data frame containing `Protein` and `LabelFactor` columns for the LBQ dataset.
#'
#' @return A patchwork-combined `ggplot` object containing multiple Venn diagrams arranged in a grid.
#'
#' @details
#' - A global ("All") Venn diagram is created for all proteins.
#' - Additional Venn diagrams are created for each condition defined in `conditionLabels`
get_venn_diagram_plot <- function(lfq_comparisonData,
                                  lbq_comparisonData,
                                  lfq_summarizationData = NULL,
                                  lbq_summarizationData = NULL) {
  venn_diagrams <- list()
  
  if (!is.null(lfq_summarizationData) &
      !is.null(lbq_summarizationData)) {
    lfq_data_tmp <- lfq_summarizationData %>% pull(Protein) %>% unique() %>% as.character()
    lbq_data_tmp <- lbq_summarizationData %>% pull(Protein) %>% unique() %>% as.character()
  } else {
    lfq_data_tmp <- lfq_comparisonData %>% pull(Protein) %>% unique() %>% as.character()
    lbq_data_tmp <- lbq_comparisonData %>% pull(Protein) %>% unique() %>% as.character()
  }
  
  venn_diagrams[["All"]] <- get_venn_diagram_plot_label(lfq_data_tmp, lbq_data_tmp, "All")
  
  for (label in names(conditionLabels)) {
    lfq_data_tmp <- lfq_comparisonData %>% dplyr::filter(LabelFactor == label) %>% pull(Protein) %>%
      unique() %>% as.character()
    lbq_data_tmp <- lbq_comparisonData %>% dplyr::filter(LabelFactor == label) %>% pull(Protein) %>%
      unique() %>% as.character()
    
    venn_diagrams[[label]] <- get_venn_diagram_plot_label(lfq_data_tmp, lbq_data_tmp, conditionLabels[label])
  }
  
  wrap_plots(venn_diagrams, ncol = 3)
}

#' Create Correlation Scatter Plots of Fold Changes Between LFQ and LBQ
#'
#' Generates scatter plots comparing log2 fold changes between LFQ and LBQ datasets,
#' filtered by significance thresholds.
#'
#' @param lfq_data A data frame with `Protein`, `LabelFactor`, `log2FC`, and `adj.pvalue`
#'        columns for the LFQ dataset.
#' @param lbq_data A data frame with `Protein`, `LabelFactor`, `log2FC`, and `adj.pvalue`
#'        columns for the LBQ dataset.
#' @param lfq_fdr_cutoff FDR cutoff value to determine significance for LFQ data.
#' @param lbq_fdr_cutoff FDR cutoff value to determine significance for LBQ data.
#' @param state A string label indicating the data state ("Before Filtering" or "After Filtering").
#'
#' @return A `ggplot` object showing scatter plots of significant proteins' fold changes across datasets.
#'
#' @details
#' - Only proteins significant in both datasets (LFQ and LBQ) are included.
#' - A linear model (forced through origin) is fitted and displayed for each facet.
#' - Red dashed lines represent zero fold change on both axes.
get_correlation_scatter_plot <- function(lfq_data,
                                         lbq_data,
                                         lfq_fdr_cutoff,
                                         lbq_fdr_cutoff,
                                         state = c("Before", "After")) {
  state <- match.arg(state)
  
  combined_set <- inner_join(
    lfq_data,
    lbq_data,
    by = c("Protein", "LabelFactor"),
    suffix = c("_lfq", "_lbq")
  )
  
  # Filter extreme outliers (defined as values outside Q1 - 20*IQR and Q3 + 20*IQR)
  # Standard boxplot outliers use 1.5*IQR, but here a broader threshold is applied to exclude only
  # the most extreme values
  log2FC_lfq_iqr <- IQR(combined_set$log2FC_lfq, na.rm = TRUE)
  log2FC_lfq_q1 <- quantile(combined_set$log2FC_lfq, 0.25, na.rm = TRUE)
  log2FC_lfq_q3 <- quantile(combined_set$log2FC_lfq, 0.75, na.rm = TRUE)
  log2FC_lfq_lower <- log2FC_lfq_q1 - 20 * log2FC_lfq_iqr
  log2FC_lfq_upper <- log2FC_lfq_q3 + 20 * log2FC_lfq_iqr
  
  log2FC_lbq_iqr <- IQR(combined_set$log2FC_lbq, na.rm = TRUE)
  log2FC_lbq_q1 <- quantile(combined_set$log2FC_lbq, 0.25, na.rm = TRUE)
  log2FC_lbq_q3 <- quantile(combined_set$log2FC_lbq, 0.75, na.rm = TRUE)
  log2FC_lbq_lower <- log2FC_lbq_q1 - 20 * log2FC_lbq_iqr
  log2FC_lbq_upper <- log2FC_lbq_q3 + 20 * log2FC_lbq_iqr
  
  combined_set_significant <- combined_set %>%
    dplyr::filter(adj.pvalue_lfq < lfq_fdr_cutoff) %>%
    dplyr::filter(adj.pvalue_lbq < lbq_fdr_cutoff) %>%
    drop_na(log2FC_lfq, log2FC_lbq) %>%
    dplyr::filter(
      log2FC_lfq >= log2FC_lfq_lower & log2FC_lfq <= log2FC_lfq_upper,
      log2FC_lbq >= log2FC_lbq_lower &
        log2FC_lbq <= log2FC_lbq_upper
    )
  
  data_stats <- combined_set_significant %>%
    group_by(LabelFactor) %>%
    summarize(PCC = cor(log2FC_lfq, log2FC_lbq, method = 'pearson')) %>%
    mutate(
      xpos = rep(-Inf, nrow(.)),
      ypos = rep(Inf, nrow(.)),
      hjustvar = rep(-0.2, nrow(.)),
      vjustvar = rep(1.7, nrow(.)),
      annotation = sub("-", "\u2212", sprintf("italic(r): %.2f", PCC))
    )
  
  base <- ggplot(combined_set_significant, aes(x = log2FC_lfq, y = log2FC_lbq))
  
  if (state == "Before") {
    base <- base +
      facet_grid(
        "Before applying the\nFC thresholds" ~ LabelFactor,
        labeller = labeller(LabelFactor = conditionLabels)
      )
  } else if (state == "After") {
    base <- base +
      facet_grid(
        "After applying the\nFC thresholds" ~ LabelFactor,
        labeller = labeller(LabelFactor = conditionLabels)
      )
  }
  
  plot <- base +
    geom_text(
      data = data_stats,
      aes(
        x = xpos,
        y = ypos,
        hjust = hjustvar,
        vjust = vjustvar,
        label = annotation
      ),
      parse = TRUE
    ) +
    geom_point(color = hue_pal()(4)[1], alpha = 0.5) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               color = "red") +
    geom_vline(xintercept = 0,
               linetype = "dashed",
               color = "red") +
    geom_smooth(method = "lm",
                formula = y ~ 0 + x,
                color = hue_pal()(4)[3]) +
    labs(x = expression("log"[2] * " FC (LFQ)"),
         y = expression("log"[2] * " FC (LBQ)")) +
    theme(
      panel.border = element_rect(color = default_text_color, fill = NA),
      axis.ticks = element_line(color = default_text_color)
    )
}

#' Stacked Correlation Scatter Plots
#'
#' Creates a two‑panel figure: one panel for the correlation before filtering, one for after filtering.
#' Internally calls \code{get_correlation_scatter_plot()} for each data set and
#' stacks the results with \pkg{patchwork}.
#'
#' @param lfq_data_before A data frame with `Protein`, `LabelFactor`, `log2FC`, and `adj.pvalue`
#'        columns for the LFQ dataset.
#' @param lbq_data_before A data frame with `Protein`, `LabelFactor`, `log2FC`, and `adj.pvalue`
#'        columns for the LBQ dataset.
#' @param lfq_fdr_cutoff_before FDR cutoff value to determine significance for LFQ data.
#' @param lbq_fdr_cutoff_before FDR cutoff value to determine significance for LBQ data.
#' @param lfq_data_after A data frame with `Protein`, `LabelFactor`, `log2FC`, and `adj.pvalue`
#'        columns for the LFQ dataset.
#' @param lbq_data_after A data frame with `Protein`, `LabelFactor`, `log2FC`, and `adj.pvalue`
#'        columns for the LBQ dataset.
#' @param lfq_fdr_cutoff_after FDR cutoff value to determine significance for LFQ data.
#' @param lbq_fdr_cutoff_after FDR cutoff value to determine significance for LBQ data.
#'
#' @return A patchwork-combined `ggplot` object combining the two \code{ggplot}s.
get_combined_correlation_scatter_plot <- function(lfq_data_before,
                                                  lbq_data_before,
                                                  lfq_fdr_cutoff_before,
                                                  lbq_fdr_cutoff_before,
                                                  lfq_data_after,
                                                  lbq_data_after,
                                                  lfq_fdr_cutoff_after,
                                                  lbq_fdr_cutoff_after) {
  before_plot <- get_correlation_scatter_plot(
    lfq_data_before,
    lbq_data_before,
    lfq_fdr_cutoff_before,
    lbq_fdr_cutoff_before,
    state = "Before"
  )
  after_plot <- get_correlation_scatter_plot(
    lfq_data_after,
    lbq_data_after,
    lfq_fdr_cutoff_after,
    lbq_fdr_cutoff_after,
    state = "After"
  )
  
  (before_plot / after_plot) +
    plot_layout(heights = c(1, 1), axes = "collect")
}

#' Generate a Venn diagram showing overlap of significant proteins across conditions
#'
#' This function creates a Venn diagram to visualize the overlap of significant proteins between
#' different experimental conditions. The conditions compared are:
#' "EE_vs_DMSO", "EE+LNG_vs_DMSO", "LNG_vs_DMSO", and "S-23_vs_DMSO".
#'
#' @param data A data.frame containing protein information with columns `Protein`, `Label`, and `adj.pvalue`.
#'
#' @return A `ggplot` object showing the Venn diagram of significant protein overlaps between conditions.
get_condition_overlap_venn_diagram_plot <- function(data) {
  set1 <- data %>% dplyr::filter(Label == "EE_vs_DMSO") %>% pull(Protein) %>%
    unique() %>% as.character()
  set2 <- data %>% dplyr::filter(Label == "LNG_vs_DMSO") %>% pull(Protein) %>%
    unique() %>% as.character()
  set3 <- data %>% dplyr::filter(Label == "EE+LNG_vs_DMSO") %>% pull(Protein) %>%
    unique() %>% as.character()
  set4 <- data %>% dplyr::filter(Label == "S-23_vs_DMSO") %>% pull(Protein) %>%
    unique() %>% as.character()
  
  dap_sets <- list(
    "EE" = set1,
    "LNG" = set2,
    "EE+LNG" = set3,
    "S-23" = set4
  )
  
  ggVennDiagram(
    dap_sets,
    set_size = 2,
    label_size = 2,
    edge_size = 0.5
  ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = default_font_size, color = default_text_color),
      legend.position = "none"
    )
}