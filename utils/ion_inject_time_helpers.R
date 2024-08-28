library(dplyr, quietly = TRUE)
library(doParallel, quietly = TRUE)
library(foreach, quietly = TRUE)
library(ggplot2, quietly = TRUE)

source("../utils/settings.R")
source("../utils/helpers.R")

# Helper function to get the ion injection time
get_ion_inject_time_values <- function(ms2_scan_chunks) {
  ms2_psm_ion_inject_time <- list()

  for (file in names(ms2_scan_chunks)) {
    ms2_scan_chunk <- ms2_scan_chunks[[file]]
    ms2_psm_ion_inject_time[[file]] <- ms2_scan_chunk %>%
      select(matches("^Scan"), matches("^Ion"))
  }

  ms2_psm_ion_inject_time_df <- do.call(rbind, Map(cbind, ms2_psm_ion_inject_time, File = as.numeric(names(ms2_psm_ion_inject_time))))

  return(ms2_psm_ion_inject_time_df)
}

# Helper function to plot the ion injection time
plot_ion_injection_time <- function(df, df_stats) {
  plot <- ggplot(df, aes(x = Ion.Inject.Time)) +
    facet_wrap(vars(File), ncol = 4) +
    geom_density(fill = "Blue", alpha = 0.2, aes(y = after_stat(scaled))) +
    geom_vline(data = df_stats, aes(xintercept = perc5), color = "Blue", alpha = 0.25, linetype = "dashed", linewidth = 0.5) +
    geom_vline(data = df_stats, aes(xintercept = perc95), color = "Blue", alpha = 0.25, linetype = "dashed", linewidth = 0.5) +
    geom_vline(data = df_stats, aes(xintercept = mean), color = "Red", alpha = 0.25, linetype = "dashed", linewidth = 0.5) +
    geom_text(data = df_stats, aes(x = +Inf, y = +Inf, label = sprintf("%.2f %%", mean), hjust = 1, vjust = 1)) +
    labs(
      title = paste0("Ion Inject. Time"),
      x = "Ion Inject. Time [ms]",
      y = "Density"
    ) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.minor = element_blank()
    )

  return(plot)
}