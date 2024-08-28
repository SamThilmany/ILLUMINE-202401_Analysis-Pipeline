library(ggplot2, quietly = TRUE)

source("../utils/helpers.R")
source("../utils/ion_inject_time_helpers.R")

setClass(
  Class = "IonInjectTime",
  slots = c(
    ms2PsmScans = "data.frame",
    ms2ScanChunks = "list",
    ms2PsmIonInjectTime_df = "data.frame",
    ms2PsmIonInjectTime_df_quantile = "data.frame"
  ),
  prototype = list(
    ms2PsmScans = data.frame(),
    ms2ScanChunks = list(),
    ms2PsmIonInjectTime_df = data.frame(),
    ms2PsmIonInjectTime_df_quantile = data.frame()
  )
)

setMethod(
  f = "initialize",
  signature = "IonInjectTime",
  definition = function(.Object, ms2_scans = data.frame()) {
    .Object@ms2PsmScans <- get_ms2_psm_scans(ms2_scans)
    .Object@ms2ScanChunks <- get_ms2_scan_chunks(.Object@ms2PsmScans, subset = FALSE)
    .Object@ms2PsmIonInjectTime_df <- get_ion_inject_time_values(.Object@ms2ScanChunks)
    .Object@ms2PsmIonInjectTime_df_quantile <- get_df_stats(.Object@ms2PsmIonInjectTime_df, group_by_column = "File", column_name = "Ion.Inject.Time")

    return(.Object)
  }
)

setGeneric(
  name = "plot",
  def = function(.Object, ...) {
    standardGeneric("plot")
  }
)

setMethod(
  f = "plot",
  signature = "IonInjectTime",
  definition = function(.Object, directory = "./") {
    plot <- plot_ion_injection_time(
      df = .Object@ms2PsmIonInjectTime_df,
      df_stats = .Object@ms2PsmIonInjectTime_df_quantile
    )

    save_and_show_plot("psm_ionInjectTime_density_plot", directory, plot, plot_height = 2 * default_height)
  }
)