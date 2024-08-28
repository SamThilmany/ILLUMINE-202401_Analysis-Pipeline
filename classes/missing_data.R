source("../utils/helpers.R")
source("../utils/missing_data_helpers.R")

setClass(
  Class = "MissingData",
  slots = c(
    dataFrame = "data.frame",
    missingValues_file_df = "data.frame",
    missingValues_channel_df = "data.frame",
    missingValues_replicate_df = "data.frame"
  ),
  prototype = list(
    dataFrame = data.frame(),
    missingValues_file_df = data.frame(),
    missingValues_channel_df = data.frame(),
    missingValues_replicate_df = data.frame()
  )
)

setMethod(
  f = "initialize",
  signature = "MissingData",
  definition = function(.Object, data_frame = data.frame(), by_channel = FALSE, by_replicate = FALSE) {
    .Object@dataFrame <- data_frame

    if (!is.logical(by_channel) || length(by_channel) != 1) {
      stop(sprintf("%s is not a valid argument for `by_channel` in the `performAnalysis` method.", by_channel))
    }

    if (!is.logical(by_replicate) || length(by_replicate) != 1) {
      stop(sprintf("%s is not a valid argument for `by_replicate` in the `performAnalysis` method.", by_replicate))
    }

    .Object@missingValues_file_df <- get_files_missing_values(.Object@dataFrame)
    .Object@missingValues_channel_df <- if (by_channel) get_channels_missing_values(.Object@dataFrame) else data.frame()
    .Object@missingValues_replicate_df <- if (by_replicate) get_replicate_missing_values(.Object@dataFrame) else data.frame()

    return(.Object)
  }
)

setGeneric(
  name = "plotByFile",
  def = function(.Object, ...) {
    standardGeneric("plotByFile")
  }
)

setMethod(
  f = "plotByFile",
  signature = "MissingData",
  definition = function(.Object, directory = "./", node_name, node_slug) {
    plot <- plot_missing_data(
      df = .Object@missingValues_file_df,
      x_col = "File",
      x_name = "File",
      node_name = node_name
    )

    save_and_show_plot(sprintf("missing_values_%s_byFiles", node_slug), directory, plot)
  }
)

setGeneric(
  name = "plotByChannel",
  def = function(.Object, ...) {
    standardGeneric("plotByChannel")
  }
)

setMethod(
  f = "plotByChannel",
  signature = "MissingData",
  definition = function(.Object, directory = "./", node_name, node_slug) {
    plot <- plot_missing_data(
      df = .Object@missingValues_channel_df,
      x_col = "Channel",
      x_name = "Channel",
      node_name = node_name
    )

    save_and_show_plot(sprintf("missing_values_%s_byChannels", node_slug), directory, plot)
  }
)

setGeneric(
  name = "plotByReplicate",
  def = function(.Object, ...) {
    standardGeneric("plotByReplicate")
  }
)

setMethod(
  f = "plotByReplicate",
  signature = "MissingData",
  definition = function(.Object, directory = "./", node_name, node_slug) {
    plot <- plot_missing_data(
      df = .Object@missingValues_replicate_df,
      x_col = "Tech.Replicate",
      x_name = "Replicate",
      node_name = node_name
    )

    save_and_show_plot(sprintf("missing_values_%s_byReplicates", node_slug), directory, plot)
  }
)