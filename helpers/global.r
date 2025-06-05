#' Save and Display a ggplot Object in Both PNG and PDF Formats
#'
#' Saves a ggplot object to both PNG and PDF formats and displays the PNG version inline.
#' Useful for analysis workflows where both high-quality vector output (PDF) and inline visualization (PNG)
#' are desired.
#'
#' @param filename A character string specifying the base filename (without extension) for the saved plot.
#' @param plot A `ggplot` object to be saved and displayed.
#' @param dir A character string specifying the directory to save the plot in.
#'        Default is `"./graphs/"`.
#' @param show_width Width (in pixels) for displaying the image inline.
#'        Default is proportional to `plot_width`.
#' @param plot_width Width of the plot in units defined by `default_unit`. Used for saving the image.
#'        Default is `default_width`.
#' @param plot_height Height of the plot in units defined by `default_unit`. Used for saving the image.
#'        Default is `default_height`.
#' @param suppress_warnings Logical. If `TRUE`, suppresses warnings during saving and rendering.
#'        Default is `TRUE`.
#'
#' @details
#' Requires global variables: `default_width`, `default_height`, `default_dpi`, and `default_unit`.
#' Also depends on `display_png()` for rendering the image inline.
#' PNGs are used for inline display; PDFs provide publication-quality output.
#'
#' @return Invisibly returns `NULL`. The function is called for its side effects.
save_and_show_plot <- function(filename,
                               plot,
                               dir = "./graphs/",
                               show_width = 500 / default_width * plot_width,
                               plot_width = default_width,
                               plot_height = default_height,
                               suppress_warnings = TRUE) {
  save_and_show <- function() {
    # Save PNG
    ggsave(
      filename = paste0(filename, ".png"),
      path = dir,
      plot = plot,
      dpi = default_dpi,
      width = plot_width,
      height = plot_height,
      unit = default_unit,
      device = "png",
      create.dir = TRUE
    )
    
    # Save PDF
    ggsave(
      filename = paste0(filename, ".pdf"),
      path = dir,
      plot = plot,
      width = plot_width,
      height = plot_height,
      unit = default_unit,
      device = "pdf",
      create.dir = TRUE
    )
    
    # Show PNG inline
    display_png(file = file.path(dir, paste0(filename, ".png")), width = show_width)
  }
  
  if (suppress_warnings) {
    suppressWarnings(save_and_show())
  } else {
    save_and_show()
  }
  
  invisible(NULL)
}

#' Save a data table as a CSV file and display it
#'
#' This function saves a given data.frame to a CSV file in the specified directory
#' and then displays the table-head inline.
#'
#' @param filename A character string specifying the name of the CSV file (without extension).
#' @param table A data.frame containing the table to be saved and displayed.
#' @param dir A character string specifying the directory where the CSV file will be saved.
#'        Default is "./results/".
#'
#' @return Invisibly returns the input `table`. The head of the table is also displayed in the
#'         console/viewer and saved to a CSV file.
save_and_show_table <- function(filename, table, dir = "./results/") {
  file_path <- file.path(dir, paste0(filename, ".csv"))
  
  write.csv(table, file_path)
  display(head(table))
  
  invisible(NULL)
}
