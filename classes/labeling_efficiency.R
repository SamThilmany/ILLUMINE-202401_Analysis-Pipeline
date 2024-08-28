library(ggplot2, quietly = TRUE)

source("../utils/helpers.R")
source("../utils/labeling_efficiency_helpers.R")

setClass(
  Class = "LabelingEfficiency",
  slots = c(
    dataFrame = "data.frame",
    dataFrameLabeled = "data.frame",
    reporterIonIntensityDist = "data.frame",
    labelingEfficiency = "data.frame"
  ),
  prototype = list(
    dataFrame = data.frame(),
    dataFrameLabeled = data.frame(),
    reporterIonIntensityDist = data.frame(),
    labelingEfficiency = data.frame()
  )
)

setMethod(
  f = "initialize",
  signature = "LabelingEfficiency",
  definition = function(.Object, data_frame = data.frame()) {
    .Object@dataFrame <- data_frame
    .Object@dataFrameLabeled <- get_labeled_df(.Object@dataFrame)
    .Object@reporterIonIntensityDist <- get_reporter_ion_intensity_dist(.Object@dataFrame)
    .Object@labelingEfficiency <-  get_labeling_efficiency(.Object@dataFrame, .Object@dataFrameLabeled)

    return(.Object)
  }
)

setGeneric(
  name = "plotRiiRatio",
  def = function(.Object, ...) {
    standardGeneric("plotRiiRatio")
  }
)

setMethod(
  f = "plotRiiRatio",
  signature = "LabelingEfficiency",
  definition = function(.Object, directory = "./", node_name, node_slug) {
    plot <- ggplot(.Object@reporterIonIntensityDist, aes(x = log2.RII.Ratios)) +
      facet_wrap(vars(Channel), nrow = 2) +
      geom_vline(aes(xintercept = median(log2.RII.Ratios)), color = "Blue", alpha = 0.25, linetype = "dashed", linewidth = 0.5) +
      geom_density(fill = "Blue", alpha = 0.2) +
      labs(
        title = sprintf("Reporter Ion Intensity (%s)", node_name),
        x = expression("log"[2] * "(RII/Mean RII)"),
        y = "Density"
      ) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.minor = element_blank()
      )

    save_and_show_plot(sprintf("RII_distribution_%s", node_slug), directory, plot)
  }
)

setGeneric(
  name = "plotLabelingEfficiency",
  def = function(.Object, ...) {
    standardGeneric("plotLabelingEfficiency")
  }
)

setMethod(
  f = "plotLabelingEfficiency",
  signature = "LabelingEfficiency",
  definition = function(.Object, directory = "./", node_name, node_slug) {
    plot <- plot_labeling_efficiency(.Object@labelingEfficiency, node_name)

    save_and_show_plot(sprintf("labeling_efficiency_%s", node_slug), directory, plot)
  }
)