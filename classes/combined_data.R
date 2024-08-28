library(dplyr, quietly = TRUE)
library(stringr, quietly = TRUE)

source("../utils/helpers.R")
source("../utils/renamer.R")
source("../utils/combined_data_helpers.R")

setClass(
  Class = "CombinedData",
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
  signature = "CombinedData",
  definition = function(
    .Object,
    ...,
    gene_disease_associations = data.frame()
  ) {
    instances <- list(...)

    protein_df_list <- lapply(instances, function(instance) {
      instance@protein_df
    })

    ratio_df_list <- lapply(instances, function(instance) {
      instance@ratio_df
    })

    .Object@protein_df <- bind_rows(protein_df_list)
    .Object@ratio_df <- bind_rows(ratio_df_list)
    .Object@ratio_df_regulated <- get_ratio_df_regulated(.Object@ratio_df)
    .Object@ratio_df_consistent <- get_ratio_df_consistent(.Object@ratio_df)
    .Object@geneDiseaseAssociations <- gene_disease_associations

    write.csv(
      .Object@ratio_df_regulated,
      file.path(
        biological_analysis_directory,
        "ratio_df_regulated.csv"
      ),
      row.names = FALSE
    )

    write.csv(
      .Object@ratio_df_consistent,
      file.path(
        biological_analysis_directory,
        "ratio_df_consistent.csv"
      ),
      row.names = FALSE
    )

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
  signature = "CombinedData",
  definition = function(.Object, directory = "./", show_labels = FALSE, match_gene_list = FALSE) {
    ratios <- unique(.Object@ratio_df$Comparison)

    for (ratio in ratios) {
      df <- .Object@ratio_df %>% filter(Comparison == ratio)
      treatment <- unique(df$Treatment)
      control <- unique(df$Control)

      if (control == "CTRL") next

      plot <- plot_volcano(
        df = df,
        treatment = treatment,
        control = control,
        show_labels = show_labels,
        match_gene_list = match_gene_list
      )

      save_and_show_plot(sprintf("volcano_%s_vs_%s_%s%s", treatment, control, ifelse(show_labels, "_labeled", ""), ifelse(match_gene_list, "_inGeneList", "")), directory, plot)

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
  name = "plotVennTreatments",
  def = function(.Object, ...) {
    standardGeneric("plotVennTreatments")
  }
)

setMethod(
  f = "plotVennTreatments",
  signature = "CombinedData",
  definition = function(.Object, directory = "./", treatments = "all", set = "consistent", match_gene_list = FALSE, print_intersections = FALSE) {
    if (! set %in% c("all", "regulated", "consistent")) {
      stop(sprintf("%s is not a valid argument for `set` in the `plotVennTreatments` method.", set))
    }

    if (is.character(treatments) && length(treatments) == 1) {
      treatments <- list(treatments)
    }

    if (is.character(treatments) && length(treatments) > 1) {
      treatments <- as.list(treatments)
    }

    if (set == "all") {
      df <- .Object@ratio_df
      slug <- "all_proteins"
      name <- "Quantified proteins"
    } else if (set == "regulated") {
      df <- .Object@ratio_df_regulated
      slug <- "regulated_proteins"
      name <- "Regulated proteins"
    } else {
      df <- .Object@ratio_df_consistent
      slug <- "consistent_proteins"
      name <- "Consistent proteins"
    }

    plot <- plot_venn_treatments(df, name, treatments, match_gene_list, print_intersections)

    save_and_show_plot(sprintf("venn_plot_%s_treatments_%s%s", slug, set, ifelse(match_gene_list, "_inGeneList", "")), directory, plot)
  }
)

setGeneric(
  name = "plotVennUpDownTreatment",
  def = function(.Object, ...) {
    standardGeneric("plotVennUpDownTreatment")
  }
)

setMethod(
  f = "plotVennUpDownTreatment",
  signature = "CombinedData",
  definition = function(.Object, directory = "./", treatments, set = "consistent", match_gene_list = FALSE, print_intersections = FALSE) {
    if (is.character(treatments) && length(treatments) == 1) {
      treatments <- list(treatments)
    }

    if (is.character(treatments) && length(treatments) > 1) {
      treatments <- as.list(treatments)
    }

    if (length(treatments) != 2) {
      stop(sprintf("The `treatments` parameter in the `plotVennUpDownTreatment` method needs to contain two treatments, but %s were given.", length(treatments)))
    }

    if (! set %in% c("all", "regulated", "consistent")) {
      stop(sprintf("%s is not a valid argument for `set` in the `plotVennUpDownTreatment` method.", set))
    }

    if (set == "all") {
      df <- .Object@ratio_df
      slug <- "all_proteins"
      name <- "Quantified proteins"
    } else if (set == "regulated") {
      df <- .Object@ratio_df_regulated
      slug <- "regulated_proteins"
      name <- "Regulated proteins"
    } else {
      df <- .Object@ratio_df_consistent
      slug <- "consistent_proteins"
      name <- "Consistent proteins"
    }

    plot <- plot_venn_up_down_treatments(df, name, treatments, match_gene_list, print_intersections)

    save_and_show_plot(sprintf("venn_plot_%s_direction_treatment_%s-vs-%s_%s%s", slug, unlist(head(treatments, 1)), unlist(tail(treatments, 1)), set, ifelse(match_gene_list, "_inGeneList", "")), directory, plot)
  }
)

setGeneric(
  name = "plotVennQuanMethods",
  def = function(.Object, ...) {
    standardGeneric("plotVennQuanMethods")
  }
)

setMethod(
  f = "plotVennQuanMethods",
  signature = "CombinedData",
  definition = function(.Object, directory = "./", treatments, set = "consistent", print_intersections = FALSE) {
    if (is.character(treatments) && length(treatments) == 1) {
      treatments <- list(treatments)
    }

    if (is.character(treatments) && length(treatments) > 1) {
      treatments <- as.list(treatments)
    }

    if (! set %in% c("all", "regulated", "consistent")) {
      stop(sprintf("%s is not a valid argument for `set` in the `plotVennQuanMethods` method.", set))
    }

    if (set == "all") {
      df <- .Object@ratio_df
      slug <- "all_proteins"
      name <- "Quantified proteins"
    } else if (set == "regulated") {
      df <- .Object@ratio_df_regulated
      slug <- "regulated_proteins"
      name <- "Regulated proteins"
    } else {
      df <- .Object@ratio_df_consistent
      slug <- "consistent_proteins"
      name <- "Consistent proteins"
    }

    for (treatment in treatments) {
      plot <- plot_venn_quan_method(df, name, treatment, print_intersections)

      save_and_show_plot(sprintf("venn_plot_%s_quan_method_%s_%s", slug, treatment, set), directory, plot)
    }
  }
)

setGeneric(
  name = "plotVennUpDownQuanMethods",
  def = function(.Object, ...) {
    standardGeneric("plotVennUpDownQuanMethods")
  }
)

setMethod(
  f = "plotVennUpDownQuanMethods",
  signature = "CombinedData",
  definition = function(.Object, directory = "./", treatments, set = "consistent", print_intersections = FALSE) {
    if (! set %in% c("regulated", "consistent")) {
      stop(sprintf("%s is not a valid argument for `set` in the `plotVennQuanMethods` method.", set))
    }

    df <- if (set == "regulated") {
      .Object@ratio_df_regulated
    } else {
      .Object@ratio_df_consistent
    }

    for (treatment in treatments) {
      plot <- plot_venn_up_down_quan_method(df, treatment, print_intersections)

      save_and_show_plot(sprintf("venn_plot_regulated_direction_quan_method_%s_%s", treatment, set), directory, plot)
    }
  }
)
