library(dplyr, quietly = TRUE)
library(doParallel, quietly = TRUE)
library(foreach, quietly = TRUE)
library(ggplot2, quietly = TRUE)

source("../utils/settings.R")
source("../utils/helpers.R")

# Helper function to retrieve the AGC fill of a spectrum
get_agc_fill <- function(spectrum) {
  getAgcFill <- rawrr::makeAccessor(key = "AGC Fill:", returnType = "double")

  agc_fill <- getAgcFill(spectrum) * 100

  return(agc_fill)
}

get_agc_fill_values <- function(spectrum_set) {
  agc_fill_values <- sapply(spectrum_set, get_agc_fill)

  return(agc_fill_values)
}

# Function to loop through the selected spectra and read the AGC fill value
read_ms2_agc_fill <- function(ms2_scan_chunks, data_raw_list) {
  # Flatten ms2_scan_chunks to allow parallel processing of each scan
  flat_chunks <- do.call(rbind, lapply(names(ms2_scan_chunks), function(file) {
    data.frame(File = file, Scan = ms2_scan_chunks[[file]]$Scan, stringsAsFactors = FALSE)
  }))

  # Group scans by file
  flat_chunks_grouped <- flat_chunks %>%
    group_by(File) %>%
    summarize(Scans = list(Scan))

  # Define the chunk size (number of scans per chunk)
  chunk_size <- 200

  # Create a parallel backend
  cl <- makeCluster(default_cores)
  registerDoParallel(cl)

  # Export the function and other necessary objects to each worker
  clusterExport(cl, c("get_agc_fill", "get_agc_fill_values"))

  # Function to process a chunk of scans
  process_scan_chunk <- function(file, scan_chunk, data_raw_list) {
    raw_file <- data_raw_list[[file]]
    spectrum_set <- rawrr::readSpectrum(rawfile = raw_file, scan = scan_chunk)
    agc_fill_values <- get_agc_fill_values(spectrum_set)

    data.frame(
      File = file,
      Scan = scan_chunk,
      AGC.Fill = agc_fill_values,
      stringsAsFactors = FALSE
    )
  }

  # Parallel processing of each scan
  results <- foreach(group = seq_len(nrow(flat_chunks_grouped)), .combine = rbind, .packages = c("rawrr", "dplyr")) %dopar% {
    file <- flat_chunks_grouped$File[[group]]
    scan_list <- flat_chunks_grouped$Scans[[group]]
    scan_chunks <- split(scan_list, ceiling(seq_along(scan_list) / chunk_size))

    group_results <- lapply(scan_chunks, function(scan_chunk) {
      process_scan_chunk(file, scan_chunk, data_raw_list)
    })

    do.call(rbind, group_results)
  }

  # Stop the parallel backend
  stopCluster(cl)

  return(results)
}

# Helper function to plot the data
plot_agc_fill <- function(df, df_stats) {
  plot <- ggplot(df, aes(x = AGC.Fill)) +
    facet_wrap(vars(File), ncol = 4) +
    geom_density(fill = "Blue", alpha = 0.2, aes(y = after_stat(scaled))) +
    geom_vline(data = df_stats, aes(xintercept = perc5), color = "Blue", alpha = 0.25, linetype = "dashed", linewidth = 0.5) +
    geom_vline(data = df_stats, aes(xintercept = perc95), color = "Blue", alpha = 0.25, linetype = "dashed", linewidth = 0.5) +
    geom_vline(data = df_stats, aes(xintercept = mean), color = "Red", alpha = 0.25, linetype = "dashed", linewidth = 0.5) +
    geom_text(data = df_stats, aes(x = +Inf, y = +Inf, label = sprintf("%.2f %%", mean), hjust = 1, vjust = 1)) +
    labs(
      title = paste0("AGC Fill"),
      x = "AGC Fill [%]",
      y = "Density"
    ) +
    xlim(0, 100) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.minor = element_blank()
    )

  return(plot)
}