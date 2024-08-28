library(dplyr, quietly = TRUE)

source("../classes/treatment_ratios.R")
source("../utils/go_analysis_helpers.R")

setClass(
  Class = "GoAnalysisSingle",
  contains = "TreatmentRatios",
  slots = list(
    enrichment = "enrichResult"
  )
)

setMethod(
  f = "initialize",
  signature = "GoAnalysisSingle",
  definition = function(
    .Object,
    inherit,
    ratio = NA_character_,
    regulation = "up",
    ontology = "bp"
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
      stop(sprintf("%s is not a valid argument for the `regulation` parameter in the `GoAnalysisSingle` initialization method.", regulation))
    }

    if (!tolower(ontology) %in% c("bp", "mf", "cc", "all")) {
      stop(sprintf("%s is not a valid argument for the `ontology` parameter in the `GoAnalysisSingle` initialization method.", ontology))
    }

    regulated_gene_list <- prepare_gene_list(.Object@ratio_df_consistent, ratio, regulation)

    background_gene_list <- .Object@protein_df %>%
      pull(Gene.Symbol) %>%
      na.omit(Gene.Symbol) %>%
      unique()

    background_gene_list <- unlist(strsplit(background_gene_list, ";\\s*"))

    .Object@enrichment <- get_go_enrichment(regulated_gene_list, background_gene_list, ontology)

    return(.Object)
  }
)

setGeneric(
  name = "plotGO",
  def = function(.Object, ...) {
    standardGeneric("plotGO")
  }
)

setMethod(
  f = "plotGO",
  signature = "GoAnalysisSingle",
  definition = function(.Object, directory = "./", ontology = "bp", regulation, treatment, control, node_name, node_slug) {
    plot_go(
      .Object@enrichment,
      directory,
      ontology,
      regulation,
      treatment,
      control,
      node_name,
      node_slug
    )
  }
)