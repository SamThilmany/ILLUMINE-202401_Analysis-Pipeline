#' Generate and save LBQ annotation file
#'
#' Creates a tab-delimited annotation file for LBQ (Label-Based Quantification) analysis,
#' detailing experimental design including channels, conditions, fractions, technical replicates,
#' and pseudo-biological replicates across multiple TMT experiments.
#'
#' This function supports a fixed mapping of TMT channels to conditions and replicates,
#' iterates across specified fractions and technical replicates, and outputs a well-structured
#' annotation table suitable for input into downstream proteomics analysis pipelines.
#'
#' @param filename The base name (without extension) for the output annotation file.
#' @param dir The directory where the annotation file should be saved.
#'
#' @return No return value. Writes a `.tsv` annotation file to the specified directory.
create_lbq_annotation_file <- function(filename, dir) {
  fractions <- 1:8
  techReplicates <- 1:3
  conditions <- list(
    "EE" = c("126", "127N", "127C"),
    "LNG" = c("128N", "128C", "129N"),
    "EE+LNG" = c("129C", "130N", "130C"),
    "S-23" = c("131N", "131C", "132N"),
    "DMSO" = c("132C", "133N"),
    "CTRL" = c("133C", "134N")
  )
  
  rows <- list()
  row_index <- 1
  
  for (techReplicate in techReplicates) {
    for (fraction in fractions) {
      for (condition in names(conditions)) {
        labels <- conditions[[condition]]
        for (i in seq_along(labels)) {
          rows[[row_index]] <- list(
            Run = sprintf("Fraction-%s_Repl-%s.raw", fraction, techReplicate),
            Channel = labels[i],
            Condition = condition,
            Mixture = "Pool",
            TechRepMixture = techReplicate,
            Fraction = fraction,
            BioReplicate = condition
          )
          row_index <- row_index + 1
        }
      }
    }
  }
  
  lbqAnnotation <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  
  # Save to file
  write.table(
    lbqAnnotation,
    file.path(dir, filename),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}