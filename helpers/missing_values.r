#' Prepare Standardized Missingness Data
#'
#' Extracts relevant columns from a peptide dataset and annotates with a source label (e.g., "LFQ", "LBQ").
#'
#' @param data A data frame with peptide data.
#' @param sample_col Column representing the sample or run (unquoted).
#' @param peptide_col Column representing the peptide sequence (unquoted).
#' @param source A string label indicating the data source ("LFQ" or "LBQ").
#'
#' @return A standardized data frame with columns: Sample, Peptide, Missing, Source.
get_missing_data_source <- function(data,
                                    sample_col,
                                    protein_col,
                                    peptide_col,
                                    source = c("LFQ", "LBQ")) {
  source <- match.arg(source)
  data %>%
    select(Sample = {{ sample_col }},
           Protein = {{ protein_col }},
           Peptide = {{ peptide_col }},
           Missing) %>%
    mutate(Source = source)
}

#' Combine LFQ and LBQ Missingness Data
#'
#' @param lfq_data LFQ peptide-level data.
#' @param lbq_data LBQ peptide-level data.
#'
#' @return A combined data frame for both LFQ and LBQ sources.
get_missing_data <- function(lfq_data, lbq_data) {
  bind_rows(
    get_missing_data_source(lfq_data, Run, ProteinName, PeptideModifiedSequence, "LFQ"),
    get_missing_data_source(lbq_data, Channel, ProteinName, PeptideSequence, "LBQ")
  )
}

#' Count Missing and Detected Peptides
#'
#' @param data A data frame with a logical `Missing` column.
#' @return A named vector with counts for "Missing" and "Detected".
count_missing <- function(data) {
  c(Missing = sum(data$Missing),
    Detected = sum(!data$Missing))
}

#' Generate Missingness Summary Table
#'
#' Calculates total missing/detected counts and percentages for LFQ and LBQ data.
#'
#' @param lfq_data LFQ peptide data with `Missing` column.
#' @param lbq_data LBQ peptide data with `Missing` column.
#'
#' @return A tidy summary table for plotting.
get_missing_summary_data <- function(lfq_data, lbq_data) {
  data.frame(
    Method = c("LFQ", "LBQ"),
    Missing = c(count_missing(lfq_data)["Missing"], count_missing(lbq_data)["Missing"]),
    Detected = c(count_missing(lfq_data)["Detected"], count_missing(lbq_data)["Detected"])
  ) %>%
    mutate(MissingPercentage = Missing / (Missing + Detected)) %>%
    pivot_longer(
      cols = c(Missing, Detected),
      names_to = "Status",
      values_to = "Count"
    )
}

#' Compute Per-Sample Missingness
#'
#' @param missing_data Data frame with Sample, Missing, Source columns.
#'
#' @return Data frame with missing percentages and label alignment info.
get_missingness_data <- function(missing_data) {
  missing_data %>%
    group_by(Sample, Source) %>%
    summarize(MissingPercent = mean(Missing), .groups = "drop") %>%
    group_by(Source) %>%
    mutate(
      hjust = if_else(MissingPercent > 0.5 * max(MissingPercent), 1, 0),
      nudge_y = if_else(
        MissingPercent > 0.5 * max(MissingPercent),
        -0.05 * MissingPercent,
        0.05 * MissingPercent
      )
    )
}

#' Bar Plot: Missingness per Sample
#'
#' @param data Data frame with Sample, MissingPercent, hjust, Source.

#' @return A `ggplot` object showing a bar chart of missing peptides per sample and source.
get_missingness_bar_plot <- function(data) {
  ggplot(data, aes(x = Sample, y = MissingPercent)) +
    facet_wrap( ~ factor(Source, levels = c("LFQ", "LBQ")), scales = "free") +
    geom_col(fill = "gray90") +
    scale_y_continuous(labels = label_percent(), position = "left") +
    geom_text(
      aes(
        hjust = hjust,
        y = MissingPercent + nudge_y,
        label = if_else(
          MissingPercent < 0.001,
          sprintf(' %.2f ‰', MissingPercent * 1000),
          sprintf(' %.2f %%', MissingPercent * 100)
        )
      ),
      vjust = 0.5,
      size = 4 / .pt,
      angle = 90
    ) +
    labs(x = NULL, y = "Missingness") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      )
    )
}

#' Heatmap: Missingness by Peptide
#'
#' @param data Data frame with Sample, Peptide, Missing, Source.

#' @return A `ggplot` object showing a heatmap showing missing peptides per sample and source.
get_missingness_heatmap_plot <- function(data) {
  ggplot(data, aes(
    x = Sample,
    y = Peptide,
    fill = Missing,
    color = Missing
  )) +
    facet_wrap( ~ factor(Source, levels = c("LFQ", "LBQ")), scales = "free") +
    geom_raster() +
    labs(x = "Run / Channel", y = "Peptide Sequence") +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 6
      ),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text = element_blank(),
      legend.position = "none",
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      )
    ) +
    plot_layout(tag_level = "new")
}

#' Bar Plot: Overall Missingness Summary
#'
#' @param data Data frame with Method, Status, Count, and MissingPercentage.

#' @return `ggplot` object showing a bar chart of missing peptides per source.
get_missingness_summary_plot <- function(data) {
  annotation <- data %>%
    group_by(Method) %>%
    mutate(yPos = 0.1 * max(Count)) %>%
    ungroup() %>%
    transmute(Method,
              Status,
              yPos,
              label1 = Count,
              label2 = round(
                if_else(Status == "Missing", MissingPercentage, 1 - MissingPercentage) * 100,
                2
              )) %>%
    distinct()
  
  ggplot(data, aes(x = Status, y = Count, fill = Status)) +
    facet_wrap( ~ factor(Method, levels = c("LFQ", "LBQ")),
                scales = "free_y",
                ncol = 1) +
    geom_bar(stat = "identity") +
    geom_text(
      data = annotation,
      aes(
        y = yPos,
        label = sprintf('%s\n(%s %%)', comma(label1), label2)
      ),
      hjust = 0,
      vjust = 0.5,
      size = 5 / .pt,
      angle = 90
    ) +
    labs(x = NULL, y = NULL) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 6),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      )
    )
}

#' Generate Complete Missingness Visualization
#'
#' @param per_sample_df Data for per-sample missingness bar plot.
#' @param peptide_df Long-format peptide missingness data.
#' @param summary_df Summary table of missing/detected counts.
#'
#' @return A patchwork-combined `ggplot` object with all subplots.
get_complete_missingness_plot <- function(per_sample_df, peptide_df, summary_df) {
  left_plot <- (
    get_missingness_bar_plot(per_sample_df) /
      get_missingness_heatmap_plot(peptide_df)
  ) +
    plot_layout(heights = c(0.5, 1)) &
    theme(plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    ))
  
  right_plot <- get_missingness_summary_plot(summary_df)
  
  (left_plot | right_plot) +
    plot_layout(widths = c(3, 1)) +
    plot_annotation(tag_levels = c("A"))
}