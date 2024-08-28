library(dplyr, quietly = TRUE)
library(stringr, quietly = TRUE)

source("../utils/helpers.R")
source("../utils/renamer.R")
source("../utils/treatment_ratios_helpers.R")

setClass(
  Class = "TreatmentRatios",
  slots = c(
    protein_df = "data.frame",
    ratio_df = "data.frame",
    ratio_df_regulated = "data.frame",
    ratio_df_consistent = "data.frame",
    geneDiseaseAssociations = "data.frame"
  ),
  prototype = list(
    protein_df = data.frame(),
    ratio_df = data.frame(),
    ratio_df_regulated = data.frame(),
    ratio_df_consistent = data.frame(),
    geneDiseaseAssociations = data.frame()
  )
)

setMethod(
  f = "initialize",
  signature = "TreatmentRatios",
  definition = function(
    .Object,
    protein_df = data.frame(),
    quan_method = NA_character_,
    node_name = NA_character_,
    gene_disease_associations = data.frame()
  ) {
    .Object@protein_df <- protein_df
    .Object@geneDiseaseAssociations <- gene_disease_associations
    .Object@ratio_df <- get_ratio_df(
      .Object@protein_df,
      quan_method,
      node_name,
      unique(.Object@geneDiseaseAssociations$Accession)
    )
    .Object@ratio_df_regulated <- get_ratio_df_regulated(.Object@ratio_df)
    .Object@ratio_df_consistent <- get_ratio_df_consistent(.Object@ratio_df)

    return(.Object)
  }
)

setGeneric(
  name = "plotVolcanos",
  def = function(.Object, ...) {
    standardGeneric("plotVolcanos")
  }
)

setMethod(
  f = "plotVolcanos",
  signature = "TreatmentRatios",
  definition = function(.Object, directory = "./", node_name, node_slug, show_labels = FALSE, match_gene_list = FALSE) {
    ratios <- unique(.Object@ratio_df$Comparison)

    for (ratio in ratios) {
      df <- .Object@ratio_df %>% filter(Comparison == ratio)
      treatment <- unique(df$Treatment)
      control <- unique(df$Control)

      plot <- plot_volcano(
        df = df,
        treatment = treatment,
        control = control,
        node_name = node_name,
        show_labels = show_labels,
        match_gene_list = match_gene_list
      )

      save_and_show_plot(sprintf("volcano_%s_vs_%s_%s%s%s", treatment, control, node_slug, ifelse(show_labels, "_labeled", ""), ifelse(match_gene_list, "_inGeneList", "")), directory, plot)

      if (show_labels || match_gene_list) {
        up_regulated <- df %>%
          filter(Diff.Expressed == "up") %>%
          pull(ifelse(match_gene_list, Diff.Expressed.Label.Accession.List, Accession)) %>%
          unique() %>%
          na.omit()

        down_regulated <- df %>%
          filter(Diff.Expressed == "down") %>%
          pull(ifelse(match_gene_list, Diff.Expressed.Label.Accession.List, Accession)) %>%
          unique() %>%
          na.omit()

        cat("Up-regulated", fill = TRUE)
        cat(paste(up_regulated, collapse = ", "), fill = TRUE)
        cat(sprintf("Number of up-regulated proteins: %s", length(up_regulated)), fill = TRUE)

        cat("", fill = TRUE)

        cat("Down-regulated", fill = TRUE)
        cat(paste(down_regulated, collapse = ", "), fill = TRUE)
        cat(sprintf("Number of down-regulated proteins: %s", length(down_regulated)), fill = TRUE)
      }
    }
  }
)

setGeneric(
  name = "printStats",
  def = function(.Object, ...) {
    standardGeneric("printStats")
  }
)

setMethod(
  f = "printStats",
  signature = "TreatmentRatios",
  definition = function(.Object, node_name) {
    df <- .Object@ratio_df %>%
      filter(Treatment %in% c("EE", "LNG", "EE.LNG", "S.23"))

    log2_positive <- df %>%
      filter(log2.Ratio > 0) %>%
      select(Accession, log2.Ratio)

    log2_negative <- df %>%
      filter(log2.Ratio < 0) %>%
      select(Accession, log2.Ratio)

    log2_up <- df %>%
      filter(Diff.Expressed == "up") %>%
      select(Accession, log2.Ratio)

    log2_down <- df %>%
      filter(Diff.Expressed == "down") %>%
      select(Accession, log2.Ratio)

    log2_positive_string <- sprintf(
      "Positive log 2 FC in %s: %s (IQR: %s); %s features",
      node_name,
      round(mean(log2_positive$log2.Ratio), 2),
      round(IQR(log2_positive$log2.Ratio), 2),
      length(unique(log2_positive$Accession))
    )

    log2_negative_string <- sprintf(
      "Negative log 2 FC in %s: %s (IQR: %s); %s features",
      node_name,
      round(mean(log2_negative$log2.Ratio), 2),
      round(IQR(log2_negative$log2.Ratio), 2),
      length(unique(log2_negative$Accession))
    )

    log2_up_string <- sprintf(
      "Up-regulated log 2 FC in %s: %s (IQR: %s); %s features",
      node_name,
      round(mean(log2_up$log2.Ratio), 2),
      round(IQR(log2_up$log2.Ratio), 2),
      length(log2_up$Accession)
    )

    log2_down_string <- sprintf(
      "Down-regulated log 2 FC in %s: %s (IQR: %s); %s features",
      node_name,
      round(mean(log2_down$log2.Ratio), 2),
      round(IQR(log2_down$log2.Ratio), 2),
      length(log2_down$Accession)
    )

    cat(log2_positive_string, fill = TRUE)
    cat(log2_negative_string, fill = TRUE)
    cat("\n")
    cat(log2_up_string, fill = TRUE)
    cat(log2_down_string, fill = TRUE)
  }
)