library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(ggplot2, quietly = TRUE)

source("../utils/helpers.R")
source("../utils/renamer.R")

get_protein_cv <- function(df) {
  df <- df %>%
    select(matches("^Gene"), matches("^CV")) %>%
    pivot_longer(
      cols = matches("^CV"),
      names_to = "Condition",
      values_to = "CV"
    ) %>%
    mutate(Condition = gsub("cv_", "", Condition)) %>%
    drop_na()

  return(df)
}

plot_cv <- function(df, node_name, node_slug) {
  label_renamer <- function(labels) {
    sapply(labels, condition_renamer)
  }

  y_limits <- boxplot.stats(df$CV)$stats[c(1, 5)]

  plot <- ggplot(df, aes(x = Condition, y = CV, fill = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    coord_cartesian(ylim = y_limits + c(-0.05, 0.05) * diff(y_limits) / 2) +
    scale_x_discrete(labels = label_renamer) +
    labs(
      title = sprintf("Abundance Coeff. of Variance (%s)", node_name),
      x = "Condition",
      y = "CV [%]"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
  
  return(plot)
}