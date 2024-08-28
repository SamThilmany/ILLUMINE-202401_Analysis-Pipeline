library(dplyr, quietly = TRUE)

source("../classes/treatment_ratios.R")
source("../utils/kegg_analysis_helpers.R")

setClass(
  Class = "KeggAnalysis",
  contains = "TreatmentRatios",
  slots = list(
    enrichment = "enrichResult"
  )
)

setMethod(
  f = "initialize",
  signature = "KeggAnalysis",
  definition = function(
    .Object,
    inherit,
    ratio = NA_character_,
    regulation = "up"
  ) {
    if (!is(inherit, "TreatmentRatios")) {
      stop("`inherit` must be an instance of the `TreatmentRatios` class")
    }

    .Object@protein_df <- inherit@protein_df
    .Object@ratio_df <- inherit@ratio_df
    .Object@ratio_df_regulated <- inherit@ratio_df_regulated
    .Object@ratio_df_consistent <- inherit@ratio_df_consistent
    .Object@geneDiseaseAssociations <- inherit@geneDiseaseAssociations

    if (!tolower(regulation) %in% c("both", "up", "down")) {
      stop(sprintf("%s is not a valid argument for the `regulation` parameter in the `KeggAnalysis` initialization method.", regulation))
    }

    regulated_gene_list <- prepare_gene_list(.Object@ratio_df_consistent, ratio, regulation)

    background_gene_list <- .Object@protein_df %>%
      pull(Entrez.Gene.ID) %>%
      na.omit(Entrez.Gene.ID) %>%
      unique()

    background_gene_list <- unlist(strsplit(background_gene_list, ";\\s*"))

    .Object@enrichment <- get_kegg_enrichment(regulated_gene_list, background_gene_list)

    return(.Object)
  }
)

setGeneric(
  name = "plotKEGG",
  def = function(.Object, ...) {
    standardGeneric("plotKEGG")
  }
)

setMethod(
  f = "plotKEGG",
  signature = "KeggAnalysis",
  definition = function(.Object, directory = "./", regulation, treatment, control, node_name, node_slug) {
    plot_kegg(
      .Object@enrichment,
      directory,
      regulation,
      treatment,
      control,
      node_name,
      node_slug
    )
  }
)