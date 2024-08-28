source("../utils/helpers.R")
source("../utils/feature_identification_helpers.R")

setClass(
  Class = "FeatureIdentification",
  slots = c(
    a_name = "character",
    ms2DataFrame_a = "data.frame",
    psmDataFrame_a = "data.frame",
    proteinDataFrame_a = "data.frame",
    proteinDataFrame_a_complete = "data.frame",

    b_name = "character",
    ms2DataFrame_b = "data.frame",
    psmDataFrame_b = "data.frame",
    proteinDataFrame_b = "data.frame",
    proteinDataFrame_b_complete = "data.frame",

    identQuantProteinNo = "data.frame",
    psmMs2IdentRate = "data.frame"
  ),
  prototype = list(
    a_name = NA_character_,
    ms2DataFrame_a = data.frame(),
    psmDataFrame_a = data.frame(),
    proteinDataFrame_a = data.frame(),
    proteinDataFrame_a_complete = data.frame(),

    b_name = NA_character_,
    ms2DataFrame_b = data.frame(),
    psmDataFrame_b = data.frame(),
    proteinDataFrame_b = data.frame(),
    proteinDataFrame_b_complete = data.frame(),

    identQuantProteinNo = data.frame(),
    psmMs2IdentRate = data.frame()
  )
)

setMethod(
  f = "initialize",
  signature = "FeatureIdentification",
  definition = function(.Object,
    a_name = NA_character_,
    ms2_data_frame_a = data.frame(),
    psm_data_frame_a = data.frame(),
    protein_data_frame_a = data.frame(),
    protein_data_frame_a_complete = data.frame(),

    b_name = NA_character_,
    ms2_data_frame_b = data.frame(),
    psm_data_frame_b = data.frame(),
    protein_data_frame_b = data.frame(),
    protein_data_frame_b_complete = data.frame()
  ) {
    .Object@a_name <- a_name
    .Object@ms2DataFrame_a <- ms2_data_frame_a
    .Object@psmDataFrame_a <- psm_data_frame_a
    .Object@proteinDataFrame_a <- protein_data_frame_a
    .Object@proteinDataFrame_a_complete <- protein_data_frame_a_complete

    .Object@b_name <- b_name
    .Object@ms2DataFrame_b <- ms2_data_frame_b
    .Object@psmDataFrame_b <- psm_data_frame_b
    .Object@proteinDataFrame_b <- protein_data_frame_b
    .Object@proteinDataFrame_b_complete <- protein_data_frame_b_complete

    .Object@identQuantProteinNo <- get_ident_quant_protein_no(.Object@proteinDataFrame_a, .Object@proteinDataFrame_a_complete, .Object@proteinDataFrame_b, .Object@proteinDataFrame_b_complete, .Object@a_name, .Object@b_name)
    .Object@psmMs2IdentRate <- get_psm_ms2_ident_rate(.Object@ms2DataFrame_a, .Object@psmDataFrame_a, .Object@ms2DataFrame_b, .Object@psmDataFrame_b, .Object@a_name, .Object@b_name)

    return(.Object)
  }
)

setGeneric(
  name = "plotFeatures",
  def = function(.Object, ...) {
    standardGeneric("plotFeatures")
  }
)

setMethod(
  f = "plotFeatures",
  signature = "FeatureIdentification",
  definition = function(.Object, directory = "./") {
    plot <- plot_features(.Object@identQuantProteinNo)

    save_and_show_plot("feature_ident_quant", directory, plot)
  }
)

setGeneric(
  name = "plotIdRate",
  def = function(.Object, ...) {
    standardGeneric("plotIdRate")
  }
)

setMethod(
  f = "plotIdRate",
  signature = "FeatureIdentification",
  definition = function(.Object, plot_mode = "absolute", directory = "./") {
    if (!plot_mode %in% c("absolute", "relative")) {
      stop(sprintf("%s is not a valid argument for `plot_mode` in the `plotIdRate` method.", plot_mode))
    }

    plot <- plot_feature_id_rate(.Object@psmMs2IdentRate, plot_mode)

    if (plot_mode == "absolute") {
      save_and_show_plot("ms2_identification", directory, plot)
    }

    if (plot_mode == "relative") {
      save_and_show_plot("ms2_identification_percentage", directory, plot)
    }
  }
)