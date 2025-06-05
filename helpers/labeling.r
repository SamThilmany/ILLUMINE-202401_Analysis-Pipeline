#' Calculate TMT Labeling Efficiency
#'
#' Computes the percentage of peptide-spectrum matches (PSMs) that contain TMTpro modifications,
#' as a measure of labeling efficiency.
#'
#' @param psm_originalData A data frame containing PSM-level data with an `Modifications` column
#'        indicating any observed modifications.
#'
#' @return A character string summarizing the labeling efficiency as a percentage
#'         (e.g., `"The labeling efficiency is 97.45 %"`).
#'
#' @details This function assumes that TMT-labeled peptides are annotated with `"TMTpro"` in the
#'          `Modifications` column.
get_labeling_efficiency <- function(psm_originalData) {
  all_psm <- nrow(psm_originalData)
  labeled_psm <- nrow(psm_originalData %>% dplyr::filter(str_detect(Modifications, "TMTpro")))
  label_efficiency <- labeled_psm / all_psm
  
  return(sprintf("The labeling efficiency is %.2f %%", label_efficiency * 100))
}

#' Get Reporter Ion Intensity Distribution
#'
#' Prepares the reporter ion intensity data for plotting by selecting relevant columns and removing
#' missing values.
#'
#' @param reporterIonIntensityData A data frame containing at least the `Channel` and `log2Intensity` columns.
#'
#' @return A tibble with non-missing reporter ion intensities.
get_reporterIon_distribution <- function(reporterIonIntensityData) {
  reporterIonIntensity_dist <- reporterIonIntensityData %>%
    dplyr::filter(str_detect(Modifications, "TMTpro")) %>%
    dplyr::select(matches("^Abundance")) %>%
    drop_na() %>%
    pivot_longer(
      cols = matches("^Abundance"),
      names_to = "Channel",
      values_to = "Abundance"
    ) %>%
    mutate(
      Channel = gsub("Abundance\\.", "", Channel),
      RII.Ratios = Abundance / mean(Abundance),
      log2.RII.Ratios = log2(RII.Ratios)
    )
  
  return(reporterIonIntensity_dist)
}

#' Compute Median Reporter Ion Ratios and Annotation Data
#'
#' Calculates the median of log2-transformed reporter ion intensity ratios per channel and prepares
#' annotation data for plotting.
#'
#' @param reporterIonIntensity_dist A data frame containing at least the columns `Channel` and
#'        `log2.RII.Ratios`, typically output from `get_reporterIon_distribution()`
#'
#' @return A data frame with per-channel median values and additional columns for plot annotation:
#' \describe{
#'   \item{median}{The median log2 ratio for each channel.}
#'   \item{xpos, ypos}{Coordinates for placing the annotation text.}
#'   \item{hjustvar, vjustvar}{Justification variables for text placement.}
#'   \item{annotation}{Formatted label string for `geom_text`, e.g., `"~x == 0.42"` for parsed display.}
#' }
get_reporterIon_stats <- function(reporterIonIntensity_dist) {
  reporterIonIntensity_stats <- reporterIonIntensity_dist %>%
    group_by(Channel) %>%
    summarise(median = median(log2.RII.Ratios)) %>%
    mutate(
      xpos = rep(-Inf, nrow(.)),
      ypos = rep(Inf, nrow(.)),
      hjustvar = rep(-0.15, nrow(.)),
      vjustvar = rep(2, nrow(.)),
      annotation = sprintf("%s == %.2f", expression(tilde(x)), median)
    )
  
  return(reporterIonIntensity_stats)
}

#' Plot Reporter Ion Intensity Density with Annotation
#'
#' Generates a density plot for the reporter ion intensities by channel, with annotation.
#'
#' @param reporterIonIntensity_dist A data frame created by `get_reporterIon_distribution()`.
#' @param reporterIonIntensity_stats A data frame with `xpos`, `ypos`, and `annotation` columns
#'        for plot annotation.
#'
#' @return A `ggplot` object showing the density distribution with annotation.
get_reporterIon_density_plots <- function(reporterIonIntensity_dist,
                                          reporterIonIntensity_stats) {
  reporterIon_density_plots <- ggplot(reporterIonIntensity_dist, aes(x = log2.RII.Ratios)) +
    facet_wrap(vars(Channel), nrow = 2) +
    geom_density(fill = hue_pal()(4)[3], alpha = 0.2) +
    geom_vline(
      aes(xintercept = median(log2.RII.Ratios)),
      color = hue_pal()(4)[1],
      linetype = "dashed",
      linewidth = 0.5
    ) +
    geom_text(
      data = reporterIonIntensity_stats,
      aes(
        x = xpos,
        y = ypos,
        hjust = 1.1,
        vjust = 2,
        label = annotation
      ),
      parse = TRUE,
      angle = 90
    ) +
    labs(x = expression("log"[2] * " (RII/Mean RII)"), y = "Density") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      panel.grid.minor = element_blank()
    )
  
  return(reporterIon_density_plots)
}