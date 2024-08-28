library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(scales, quietly = TRUE)

source("../utils/settings.R")
source("../utils/helpers.R")

# Helper function to get the missing values per file
get_files_missing_values <- function(psm_df) {
  files_missing_values_df <- psm_df %>%
    select(matches("File"), matches("Abundance")) %>%
    pivot_longer(cols = matches("Abundance"), names_to = "Channel", values_to = "Abundance") %>%
    setNames(c("File", "Channel", "Abundance")) %>%
    group_by(File) %>%
    summarize(
      missing_values_abs = sum(is.na(Abundance)),
      missing_values_rel = missing_values_abs / n() * 100
    )

  return(files_missing_values_df)
}

# Helper function to get the missing values per channel
get_channels_missing_values <- function(psm_df) {
  channels_missing_values_df <- psm_df %>%
    select(matches("Abundance")) %>%
    pivot_longer(cols = matches("Abundance"), names_to = "Channel", values_to = "Abundance") %>%
    mutate(Channel = gsub("Abundance\\.", "TMT-", Channel)) %>%
    group_by(Channel) %>%
    summarize(
      missing_values_abs = sum(is.na(Abundance)),
      missing_values_rel = missing_values_abs / n() * 100
    )

  return(channels_missing_values_df)
}

# Helper function to get the missing values per technical replicate
get_replicate_missing_values <- function(psm_df) {
  replicate_missing_values_df <- psm_df %>%
    select(matches("Tech.Replicate"), matches("Abundance")) %>%
    pivot_longer(cols = matches("Abundance"), names_to = "Channel", values_to = "Abundance") %>%
    setNames(c("Tech.Replicate", "Channel", "Abundance")) %>%
    group_by(Tech.Replicate) %>%
    summarize(
      missing_values_abs = sum(is.na(Abundance)),
      missing_values_rel = missing_values_abs / n() * 100
    )

  return(replicate_missing_values_df)
}

# Helper function to plot missing data
plot_missing_data <- function(df, x_col, x_name, node_name) {
  plot <- ggplot(df, aes(x = !!sym(x_col), y = missing_values_rel)) +
    geom_bar(stat = "identity", fill = hue_pal()(1)) +
    labs(
      title = sprintf("Missing Values (%s)", node_name),
      x = x_name,
      y = "Missing Values [%]"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  return(plot)
}