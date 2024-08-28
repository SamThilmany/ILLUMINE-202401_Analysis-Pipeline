library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)

source("../classes/combined_data.R")
source("../utils/disgenet_analysis_helpers.R")

setClass(
  Class = "DisgenetAnalysisCombined",
  contains = "CombinedData",
  slots = list(
    enrichment = "DataGeNET.DGN",
    enrichment_df = "data.frame"
  )
)

setMethod(
  f = "initialize",
  signature = "CombinedData",
  definition = function(
    .Object,
    inherit,
    ratio = NA_character_
  ) {
    if (!is(inherit, "CombinedData")) {
      stop("`inherit` must be an instance of the `CombinedData` class")
    }

    .Object@protein_df <- inherit@protein_df
    .Object@ratio_df <- inherit@ratio_df
    .Object@ratio_df_regulated <- inherit@ratio_df_regulated
    .Object@ratio_df_consistent <- inherit@ratio_df_consistent
    .Object@geneDiseaseAssociations <- inherit@geneDiseaseAssociations

    regulated_gene_list <- prepare_gene_list(.Object@ratio_df_consistent, ratio)
    .Object@enrichment <- get_disgenet_enrichment(regulated_gene_list)

    .Object@enrichment_df <- clean_enrichment_data(.Object@enrichment)

    return(.Object)
  }
)

setGeneric(
  name = "plotDiseaseClasses",
  def = function(.Object, ...) {
    standardGeneric("plotDiseaseClasses")
  }
)

setMethod(
  f = "plotDiseaseClasses",
  signature = "DisgenetAnalysisCombined",
  definition = function(.Object, directory = "./", treatment, control) {
    plot <- plot_disease_classes(.Object@enrichment_df, treatment, control)

    save_and_show_plot(sprintf("disease_classes_%s-vs-%s", treatment, control), directory, plot)
  }
)