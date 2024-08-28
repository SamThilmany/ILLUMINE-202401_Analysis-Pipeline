library(dplyr, quietly = TRUE)
library(ggplot2, quietly = TRUE)

# Helper function to get the mass difference
get_psm_mass_diff <- function(psm_df) {
  psm_mass_diff <- psm_df %>%
    select(matches("File"), matches("mz"), matches("deltaM")) %>%
    setNames(c("File", "mz", "deltaM"))

  return(psm_mass_diff)
}

# Helper function to plot a mass difference violin plot
plot_mass_diff_violin <- function(df) {
  plot <- ggplot(df, aes(x = File, y = deltaM, fill = File)) +
    geom_violin() +
    labs(
      title = "Relative mass diff.",
      x = "File",
      y = "DeltaM [ppm]"
    ) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  return(plot)
}

# Helper function to plot a mass difference density plot
plot_mass_diff_density <- function(df, df_stats) {
  plot <- ggplot(df, aes(x = deltaM)) +
    geom_density(fill = "Blue", alpha = 0.2, aes(y = after_stat(scaled))) +
    geom_vline(data = df_stats, aes(xintercept = perc5), color = "Blue", alpha = 0.25, linetype = "dashed", linewidth = 0.5) +
    geom_vline(data = df_stats, aes(xintercept = mean), color = "Red", alpha = 0.25, linetype = "dashed", linewidth = 0.5) +
    geom_vline(data = df_stats, aes(xintercept = perc95), color = "Blue", alpha = 0.25, linetype = "dashed", linewidth = 0.5) +
    labs(
      title = "Relative mass diff.",
      x = "DeltaM [ppm]",
      y = "Density"
    ) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  return(plot)
}

# Helper function to plot a mass difference scatter plot
plot_mass_diff_scatter <- function(df) {
  plot <- ggplot(df, aes(x = mz, y = deltaM)) +
    geom_point(alpha = 0.01) +
    geom_smooth(color = "Green", linetype = "dotted", linewidth = 0.5, method = "lm", formula = "y ~ x") +
    labs(
      title = "Relative mass diff.",
      x = expression(italic("m/z")),
      y = "DeltaM [ppm]"
    )

  return(plot)
}