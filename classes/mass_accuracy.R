library(ggplot2, quietly = TRUE)

source("../utils/helpers.R")
source("../utils/mass_accuracy_helpers.R")

setClass(
  Class = "MassAccuracy",
  slots = c(
    psmDataFrame = "data.frame",
    psmMassDiff_df = "data.frame",
    psmMassDiff_df_quantile = "data.frame"
  ),
  prototype = list(
    psmDataFrame = data.frame(),
    psmMassDiff_df = data.frame(),
    psmMassDiff_df_quantile = data.frame()
  )
)

setMethod(
  f = "initialize",
  signature = "MassAccuracy",
  definition = function(.Object, psm_data_frame = data.frame()) {
    .Object@psmDataFrame <- psm_data_frame
    .Object@psmMassDiff_df <- get_psm_mass_diff(.Object@psmDataFrame)
    .Object@psmMassDiff_df_quantile <- get_df_stats(.Object@psmMassDiff_df, column_name = "deltaM")

    return(.Object)
  }
)

setGeneric(
  name = "plotViolin",
  def = function(.Object, ...) {
    standardGeneric("plotViolin")
  }
)

setMethod(
  f = "plotViolin",
  signature = "MassAccuracy",
  definition = function(.Object, directory = "./") {
    plot <- plot_mass_diff_violin(.Object@psmMassDiff_df)

    save_and_show_plot("psm_massDiff_violin", directory, plot)
  }
)

setGeneric(
  name = "plotDensity",
  def = function(.Object, ...) {
    standardGeneric("plotDensity")
  }
)

setMethod(
  f = "plotDensity",
  signature = "MassAccuracy",
  definition = function(.Object, directory = "./") {
    plot <- plot_mass_diff_density(
      df = .Object@psmMassDiff_df,
      df_stats = .Object@psmMassDiff_df_quantile
    )

    save_and_show_plot("psm_massDiff_density", directory, plot)
  }
)

setGeneric(
  name = "plotScatter",
  def = function(.Object, ...) {
    standardGeneric("plotScatter")
  }
)

setMethod(
  f = "plotScatter",
  signature = "MassAccuracy",
  definition = function(.Object, directory = "./") {
    plot <- plot_mass_diff_scatter(.Object@psmMassDiff_df)

    save_and_show_plot("psm_deltaM_vs_mass", directory, plot)
  }
)