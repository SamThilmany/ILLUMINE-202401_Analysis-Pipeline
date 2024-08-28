library(ggplot2, quietly = TRUE)
library(ggrepel, quietly = TRUE)
library(dplyr, quietly = TRUE)

source("../utils/helpers.R")
source("../utils/renamer.R")

# Helper function to draw volcano plots
plot_volcano <- function(df, treatment, control, show_labels = FALSE, match_gene_list = FALSE) {
  labels <- if (match_gene_list) df$Diff.Expressed.Label.Accession.List else df$Diff.Expressed.Label

  plot <- ggplot(df, aes(x = log2.Ratio, y = Neg.log10.p.Value, col = df$Diff.Expressed, label = labels)) +
    geom_vline(xintercept = c(-log2_fc_threshold, log2_fc_threshold), col = "gray", linetype = "dashed") +
    geom_hline(yintercept = -log10(p_value_threshold), col = "gray", linetype = "dashed") +
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
      title = sprintf("%s vs. %s", condition_renamer(treatment), condition_renamer(control)),
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

# Helper function to draw venn plots for the different treatments
plot_venn_treatments <- function(df, name, treatments, match_gene_list = FALSE, print_intersections = FALSE) {
  df_all <- df %>%
    filter(
      !Treatment %in% c("CTRL", "DMSO")
    )

  if (length(treatments) == 1 && unlist(head(treatments, 1)) == "all") {
    treatments <- unique(df_all$Treatment)
  }

  comparison <- list()
  for (treatment in treatments) {
    comparison[[condition_renamer(treatment)]] <- df_all %>%
      filter(Treatment == treatment) %>%
      pull(ifelse(match_gene_list, Diff.Expressed.Label.Accession.List, Accession)) %>%
      unique() %>%
      na.omit()
  }

  name_string <- ifelse(match_gene_list, "%s in\ndepression-related diseases\n(%s)", "%s\n(%s)")
  plot <- plot_base_venn(comparison, sprintf(name_string, name, paste(sapply(treatments, condition_renamer), collapse = " vs. ")), print_intersections)

  return(plot)
}

# Helper function to draw venn plots for treatment comparisons
plot_venn_up_down_treatments <- function(df, name, treatments, match_gene_list = FALSE, print_intersections = FALSE) {
  comparison <- list()
  for (treatment in treatments) {
    comparison[[sprintf("%s ↑", treatment)]] <- df %>%
      filter(Treatment == treatment) %>%
      filter(Diff.Expressed == "up") %>%
      pull(Accession) %>%
      unique() %>%
      na.omit()

    comparison[[sprintf("%s ↓", treatment)]] <- df %>%
      filter(Treatment == treatment) %>%
      filter(Diff.Expressed == "down") %>%
      pull(Accession) %>%
      unique() %>%
      na.omit()
  }

  plot <- plot_base_venn(comparison, sprintf("%s (%s vs. %s)", name, condition_renamer(unlist(head(treatments, 1))), condition_renamer(unlist(tail(treatments, 1)))), print_intersections)

  return(plot)
}

# Helper function to draw venn plots for the different quantification methods
plot_venn_quan_method <- function(df, name, treatment, print_intersections) {
  quan_methods <- unique(df$Quan.Method)

  comparison <- list()
  for (quan_method in quan_methods) {
    comparison[[quan_method]] <- df

    if (treatment != "all") {
      comparison[[quan_method]] <- comparison[[quan_method]] %>%
        filter(Treatment == treatment)
    }

    comparison[[quan_method]] <- comparison[[quan_method]] %>%
      filter(Quan.Method == quan_method) %>%
      pull(Accession) %>%
      unique() %>%
      na.omit()
  }

  plot <- plot_base_venn(comparison, sprintf("%s (%s)", name, condition_renamer(treatment)), print_intersections)

  return(plot)
}

# Helper function to draw venn plots for the up- and down-regulated proteins of different quantification methods
plot_venn_up_down_quan_method <- function(df, treatment, print_intersections) {
  quan_methods <- unique(df$Quan.Method)

  comparison <- list()
  for (quan_method in quan_methods) {
    comparison[[sprintf("%s ↑", quan_method)]] <- df %>%
      filter(Treatment == treatment) %>%
      filter(Diff.Expressed == "up") %>%
      filter(Quan.Method == quan_method) %>%
      pull(Accession) %>%
      unique() %>%
      na.omit()

    comparison[[sprintf("%s ↓", quan_method)]] <- df %>%
      filter(Treatment == treatment) %>%
      filter(Diff.Expressed == "down") %>%
      filter(Quan.Method == quan_method) %>%
      pull(Accession) %>%
      unique() %>%
      na.omit()
  }

  plot <- plot_base_venn(comparison, sprintf("Regulated proteins (%s)", condition_renamer(treatment)), print_intersections)

  return(plot)
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
    group_by(Accession, Treatment) %>%
    filter(all(Diff.Expressed == "up") | all(Diff.Expressed == "down")) %>%
    ungroup()

  return(ratio_df_consistent)
}