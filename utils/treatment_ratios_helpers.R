library(ggplot2, quietly = TRUE)
library(ggrepel, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(stringr, quietly = TRUE)

source("../utils/helpers.R")
source("../utils/renamer.R")

# Helper function to get a data frame of abundance ratios
get_ratio_df <- function(df, quan_method, node_name, gda_accession_list) {
  ratios <- unique(
    gsub("log2_|pValue_", "", colnames(df)[grep("log2_|pValue_", colnames(df))])
  )

  ratio_df_list <- list()
  for (ratio in ratios) {
    ratio_df_list[[ratio]] <- df %>%
      select(
        -matches("log2_|pValue_|cv_|found_in_"),
        all_of(c(
          paste0("log2_", ratio),
          paste0("pValue_", ratio)
        ))
      ) %>%
      setNames(c(
        "Accession",
        "Gene.Symbol",
        "Entrez.Gene.ID",
        "log2.Ratio",
        "p.Value"
      )) %>%
      na.omit() %>%
      mutate(
        Comparison = ratio,
        Control = str_match(ratio, "([A-za-z0-9.]+)_([A-Za-z0-9.]+)")[, 3],
        Treatment = str_match(ratio, "([A-za-z0-9.]+)_([A-Za-z0-9.]+)")[, 2],
        Neg.log10.p.Value = -log10(p.Value),
        Abs.log2.Ratio = abs(log2.Ratio),
        Diff.Expressed = ifelse(
          log2.Ratio >= log2_fc_threshold & Neg.log10.p.Value >= -log10(p_value_threshold),
          "up",
          ifelse(
            log2.Ratio <= -log2_fc_threshold & Neg.log10.p.Value >= -log10(p_value_threshold),
            "down",
            "not significant"
          )
        ),
        Diff.Expressed.Label = ifelse(
          Diff.Expressed != "not significant",
          Accession,
          NA
        ),
        Diff.Expressed.Label.Accession.List = ifelse(
          Diff.Expressed.Label %in% gda_accession_list,
          Diff.Expressed.Label,
          NA
        ),
        Quan.Method = quan_method,
        Node = node_name
      )
  }

  ratio_df <- do.call(rbind, ratio_df_list)

  return(ratio_df)
}

# Helper function to filter the ratio_df for regulated proteins
get_ratio_df_regulated <- function(ratio_df) {
  ratio_df_regulated <- ratio_df %>%
    filter(Diff.Expressed != "not significant")

  return(ratio_df_regulated)
}

# Helper function to filter the ratio_df for inconsistent regulations
get_ratio_df_consistent <- function(ratio_df) {
  ratio_df_consistent <- ratio_df %>%
    group_by(Accession) %>%
    filter(!any(Diff.Expressed != "not significant" & Treatment == "DMSO")) %>%
    ungroup() %>%
    filter(Diff.Expressed != "not significant") %>%
    group_by(Accession, Treatment, Quan.Method) %>%
    filter(all(Diff.Expressed == "up") | all(Diff.Expressed == "down")) %>%
    ungroup()

  return(ratio_df_consistent)
}

# Helper function to draw volcano plots
plot_volcano <- function(df, treatment, control, node_name, show_labels = FALSE, match_gene_list = FALSE) {
  labels <- if (match_gene_list) df$Diff.Expressed.Label.Accession.List else df$Diff.Expressed.Label

  plot <- ggplot(df, aes(x = log2.Ratio, y = Neg.log10.p.Value, col = df$Diff.Expressed, label = labels)) +
    geom_vline(xintercept = c(-log2_fc_threshold, log2_fc_threshold), col = "#231f2050", linetype = "dashed") +
    geom_hline(yintercept = -log10(p_value_threshold), col = "#231f2050", linetype = "dashed") +
    geom_point(size = 1, alpha = 0.5) +
    scale_color_manual(
      values = c(
        "down" = "#00AFBB",
        "not significant" = "gray",
        "up" = "#BB0C00"
      ),
      labels = c(
        "Down-regulated",
        "Not significant",
        "Up-regulated"
      )
    ) +
    coord_cartesian(
      xlim = c(-max(df$Abs.log2.Ratio), max(df$Abs.log2.Ratio)),
      ylim = c(0, max(df$Neg.log10.p.Value))
    ) +
    labs(
      title = sprintf("%s vs. %s (%s)", condition_renamer(treatment), condition_renamer(control), node_name),
      x = expression("log"[2] * "FC"),
      y = expression("-log"[10] * "p-value")
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(linetype = "solid", fill = NA),
      legend.position = "none"
    )

  if (show_labels || match_gene_list) {
    plot <- plot + geom_text_repel(max.overlaps = Inf, size = default_font_size / (1.5 * .pt))
  }

  return(plot)
}