#' Prepare abundance‐density data for visualisation
#'
#' Combines the raw and normalised abundance tables into a single long data frame that is ready to be plotted.
#'
#' @param data_beforeNormalization A data frame **before** normalisation. Must contain at least the columns
#'        given in \code{sampleColName_before} and \code{abundanceColName_before}.
#' @param data_afterNormalization  A data frame **after** normalisation. Must contain at least the columns
#'        given in \code{sampleColName_after} and \code{abundanceColName_after}.
#' @param sampleColName_before Column representing the sample or run in \code{data_beforeNormalization} (unquoted).
#' @param sampleColName_after  Column representing the sample or run in \code{data_afterNormalization} (unquoted).
#' @param sampleLevels A character vector giving the desired order of the samples on the x‑axis.
#' @param abundanceColName_before Column representing the abundance in \code{data_beforeNormalization}.
#' @param abundanceColName_after  Column representing the abundance in \code{data_afterNormalization}.
#' @param source A string label indicating the data source ("LFQ" or "LBQ").
#'
#' @return A tibble with the columns `Sample`, `Abundance, `Normalization`, and `Source`.
#'
#' @details Rows with missing abundance values are dropped.
get_abundance_density_data <- function(data_beforeNormalization,
                                       data_afterNormalization,
                                       sampleColName_before,
                                       sampleColName_after,
                                       sampleLevels,
                                       abundanceColName_before,
                                       abundanceColName_after,
                                       source = c("LFQ", "LBQ")) {
  source <- match.arg(source)
  
  data_before <- data_beforeNormalization %>%
    select(Sample = {{ sampleColName_before }}, Abundance = {{ abundanceColName_before }}) %>%
    mutate(Sample = factor(Sample, levels = sampleLevels),
           Normalization = "Before Normalization")
  
  data_after <- data_afterNormalization %>%
    select(Sample = {{ sampleColName_after }}, Abundance = {{ abundanceColName_after}}) %>%
    mutate(Sample = factor(Sample, levels = sampleLevels),
           Normalization = "After Normalization")
  
  abundance_density_data <- bind_rows(data_before, data_after) %>%
    mutate(Source = source,
           Normalization = factor(
             Normalization,
             levels = c("Before Normalization", "After Normalization")
           ),) %>%
    drop_na(Abundance)
  
  return(abundance_density_data)
}

#' Density‐heatmap of abundances per sample
#'
#' Produces a 2D density raster (essentially a heatmap) of log‑intensity
#' distribution per run, facetted by data source and by normalisation state.
#'
#' @param abundance_density_data A data frame created by \code{get_abundance_density_data()}.
#'
#' @return A `ggplot` object showing the abundance density per sample and source.
get_abundance_density_plot <- function(abundance_density_data) {
  ggplot(abundance_density_data, aes(x = Sample, y = Abundance)) +
    facet_grid(Source ~ Normalization) +
    stat_density(aes(fill = after_stat(density)),
                 geom = "raster",
                 position = "identity") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Run", y = "Abundance (arb. unit)", fill = "Density") +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 6
      ),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.spacing = unit(1, "lines")
    )
}

#' Stacked abundance‑density plots for LFQ and LBQ
#'
#' Creates a two‑panel figure: one panel for LFQ, one for LBQ.
#' Internally calls \code{get_abundance_density_plot()} for each data set and
#' stacks the results with \pkg{patchwork}.
#'
#' @param lfq_abundance_density_data A data frame produced by
#'   \code{get_abundance_density_data(source = "LFQ")}.
#' @param lbq_abundance_density_data A data frame produced by
#'   \code{get_abundance_density_data(source = "LBQ")}.
#'
#' @return A patchwork-combined `ggplot` object combining the two \code{ggplot}s.
get_complete_abundance_density_plot <- function(lfq_abundance_density_data,
                                                lbq_abundance_density_data) {
  lfq_plot <- get_abundance_density_plot(lfq_abundance_density_data)
  lbq_plot <- get_abundance_density_plot(lbq_abundance_density_data)
  
  (lfq_plot / lbq_plot) +
    plot_layout(heights = c(1, 1), axes = "collect")
}