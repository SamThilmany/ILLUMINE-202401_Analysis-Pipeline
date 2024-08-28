source("../utils/helpers.R")
source("../utils/agc_fill_helpers.R")

setClass(
  Class = "AgcFill",
  slots = c(
    rawFilesPathList = "list",
    ms2PsmScans = "data.frame",
    ms2ScanChunks = "list",
    ms2PsmAgc_df = "data.frame",
    ms2PsmAgc_df_quantile = "data.frame"
  ),
  prototype = list(
    rawFilesPathList = list(),
    ms2PsmScans = data.frame(),
    ms2ScanChunks = list(),
    ms2PsmAgc_df = data.frame(),
    ms2PsmAgc_df_quantile = data.frame()
  )
)

setMethod(
  f = "initialize",
  signature = "AgcFill",
  definition = function(.Object, raw_files_path_list = list(), ms2_scans = data.frame()) {
    .Object@rawFilesPathList <- raw_files_path_list
    .Object@ms2PsmScans <- get_ms2_psm_scans(ms2_scans)
    .Object@ms2ScanChunks <- get_ms2_scan_chunks(.Object@ms2PsmScans, subset = TRUE)
    .Object@ms2PsmAgc_df <- read_ms2_agc_fill(.Object@ms2ScanChunks, .Object@rawFilesPathList)
    .Object@ms2PsmAgc_df_quantile <- get_df_stats(.Object@ms2PsmAgc_df, group_by_column = "File", column_name = "AGC.Fill")

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
  signature = "AgcFill",
  definition = function(.Object, directory = "./") {
    plot <- plot_agc_fill(
      df = .Object@ms2PsmAgc_df,
      df_stats = .Object@ms2PsmAgc_df_quantile
    )

    save_and_show_plot("psm_agcFill_density_plot", directory, plot, plot_height = 2 * default_height)
  }
)