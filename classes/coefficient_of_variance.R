source("../utils/coefficient_of_variance_helpers.R")

setClass(
  Class = "CoefficientOfVariance",
  slots = list(
    protein_df = "data.frame"
  ),
  prototype = list(
    protein_df = data.frame()
  )
)

setMethod(
  f = "initialize",
  signature = "CoefficientOfVariance",
  definition = function(.Object, protein_df = data.frame()) {
    .Object@protein_df <- get_protein_cv(protein_df)

    return(.Object)
  }
)

setGeneric(
  name = "plotCV",
  def = function(.Object, ...) {
    standardGeneric("plotCV")
  }
)

setMethod(
  f = "plotCV",
  signature = "CoefficientOfVariance",
  definition = function(.Object, directory = "./", node_name, node_slug) {
    plot <- plot_cv(
      .Object@protein_df,
      node_name,
      node_slug
    )

    save_and_show_plot(sprintf("protein_cv_density_%s", node_slug), directory, plot)
  }
)