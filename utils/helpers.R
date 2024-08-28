library(IRdisplay, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(ggplot2, quietly = TRUE)

source("../utils/settings.R")

# Slect all the MS2 spectra that correspond to at least one PSM
get_ms2_psm_scans <- function(ms2_df) {
  ms2_psm_scans <- ms2_df %>%
    select(matches("File"), matches("^Scan"), matches("PSM"), matches("^Ion")) %>%
    setNames(c("File.ID", "Scan", "PSM", "Ion.Inject.Time")) %>%
    na.omit() %>%
    filter(PSM > 0)

  return(ms2_psm_scans)
}

# Split the list of scans into chunks according to the File.ID
get_ms2_scan_chunks <- function(ms2_psm_scans, subset = FALSE) {
  if (!is.logical(subset) || length(subset) != 1) {
    stop(sprintf("%s is not a valid argument for `subset` in the `get_ms2_scan_chunks` method.", subset))
  }

  ms2_scan_chunks <- list()

  for (file in unique(ms2_psm_scans$File.ID)) {
    file_scans <- ms2_psm_scans %>%
      filter(File.ID == file)

    if (subset && !is.null(random_subset)) {
      ms2_scan_chunks[[file]] <- file_scans %>%
        sample_n(min(random_subset, nrow(file_scans)), replace = FALSE)
    } else {
      ms2_scan_chunks[[file]] <- file_scans
    }
  }

  return(ms2_scan_chunks)
}

# Helper function to calculate df statistics
get_df_stats <- function(df, group_by_column = NULL, column_name) {
  df_quantile <- df

  if (is.character(group_by_column)) {
    df_quantile <- df_quantile %>%
      group_by(!!sym(group_by_column))
  }

  df_quantile <- df_quantile %>%
    summarise(
      perc5 = quantile(.data[[column_name]], 0.05),
      median = median(.data[[column_name]]),
      mean = mean(.data[[column_name]]),
      perc95 = quantile(.data[[column_name]], 0.95)
    )

  return(df_quantile)
}

# Helper function to save and show a plot
save_and_show_plot <- function(filename, dir, plot, show_width = 500, plot_width = default_width, plot_height = default_height, suppress_warnings = TRUE) {
  save_and_show <- function() {
    ggsave(filename = paste0(filename, ".png"), path = dir, plot = plot, dpi = default_dpi, width = plot_width, height = plot_height, unit = default_unit, create.dir = TRUE)
    display_png(file = paste0(dir, "/", filename, ".png"), width = show_width)
  }

  if (suppress_warnings) {
    suppressWarnings(save_and_show())
  } else {
    save_and_show()
  }
}

# Helper function to retrieve gene-disease-associations from DISGENET
get_gene_disease_associations <- function(diseases_of_interest) {
  diseases_of_interest <- paste0("MESH_", diseases_of_interest)

  gene_disease_associations <- disgenet2r::disease2gene(
    disease = diseases_of_interest,
    database = "CURATED",
    score = c(0.25, 1),
    warnings = FALSE
  )

  gene_disease_associations <- gene_disease_associations@qresult %>%
    select(matches("gene_symbol"), matches("uniprotids")) %>%
    setNames(c("Gene.Symbol", "Accession")) %>%
    unnest_longer(Accession)
}

# Helper function to plot Venn plots
plot_base_venn <- function(comparison, title, print_intersections = FALSE) {
  comparison_no_names <- comparison
  names(comparison_no_names) <- c(rep("", length(comparison)))

  venn_no_names <- ggVennDiagram::Venn(comparison_no_names)
  venn <- ggVennDiagram::Venn(comparison)
  data_no_names <- ggVennDiagram::process_data(venn_no_names)
  data <- ggVennDiagram::process_data(venn)

  plot_data <- if (length(comparison) == 4) {
    data_no_names
  } else {
    data
  }

  plot <- ggVennDiagram::plot_venn(
    plot_data,
    set_size = default_font_size / .pt,
    label_size = default_font_size / (1.5 * .pt),
    edge_size = 0.5,
    label_percent_digit = 1
  ) +
    ggtitle(title) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0.5, size = default_font_size * 14 / 11),
      legend.title = element_text(size = default_font_size),
      text = element_text(size = default_font_size)
    ) +
    scale_x_continuous(expand = expansion(mult = 0.1)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
    labs(fill = "Count")

  if (length(comparison) == 4) {
    plot <- plot +
      geom_text(
        aes(ifelse(X < 0.5, X - (X * 0.4), X + ((1 - X) * 0.4)), Y, label = name),
        data = ggVennDiagram::venn_setlabel(data)
      )
  }

  if (print_intersections) {
    processed_regions <- ggVennDiagram::process_region_data(venn) %>% as.data.frame()
    ids <- processed_regions %>% pull(id)

    for (tmp_id in ids) {
      region <- processed_regions %>%
        filter(id == tmp_id) %>%
        select(name, item)

      treatment <- region$name
      items <- region$item %>% unlist() %>% paste(collapse = ", ")

      if (items != "") {
        cat(treatment, fill = TRUE)
        cat(items, fill = TRUE)
        cat("\n", fill = TRUE)
      }
    }
  }

  print(summary(comparison))

  return(plot)
}