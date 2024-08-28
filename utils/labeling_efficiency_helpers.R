library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(ggplot2, quietly = TRUE)

# Helper function to get the labeled PSMs
get_labeled_df <- function(df) {
  labeled_df <- df %>%
    filter(grepl("TMTpro", Modifications))

  return(labeled_df)
}

# Helper function to get the distribution of the reporter ion intensity
get_reporter_ion_intensity_dist <- function(labeled_df) {
  reporter_ion_intensity_dist <- labeled_df %>%
    select(matches("^Abundance")) %>%
    drop_na() %>%
    pivot_longer(cols = matches("^Abundance"), names_to = "Channel", values_to = "Abundance") %>%
    mutate(
      Channel = gsub("Abundance\\.", "", Channel),
      RII.Ratios = Abundance/mean(Abundance),
      log2.RII.Ratios = log2(RII.Ratios)
    )

  return(reporter_ion_intensity_dist)
}

# Helper function to get the labeling efficiency
get_labeling_efficiency <- function(df, labeled_df) {
  labeling_efficiency <- nrow(labeled_df) / nrow(df)

  labeling_efficiency_df <- data.frame(
    Type = c("All", "TMTpro-labeled"),
    PSM = c(100, 100 * labeling_efficiency)
  )

  return(labeling_efficiency_df)
}

# Helper function to plot the labeling efficiency
plot_labeling_efficiency <- function(df, node_name) {
  plot <- ggplot(df, aes(x = Type, y = PSM, fill = Type)) +
    geom_bar(stat = "identity", width = 0.75) +
    geom_text(aes(label = format(PSM, big.mark = ",", digits = 3)), angle = 90, vjust = 0.5, hjust = 1.25) +
    labs(
      title = sprintf("Labeling Efficiency (%s)", node_name),
      x = "",
      y = "N° PSM [% from total]"
    ) +
    theme(
      legend.position = "none"
    )

  return(plot)
}