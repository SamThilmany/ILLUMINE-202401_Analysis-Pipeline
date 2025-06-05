#' Load and filter proteomics quantification data
#'
#' Reads a tab-delimited text file containing protein quantification data and applies
#' filtering steps depending on the data type (`LFQ` or `LBQ`). For `LFQ` data,
#' missing precursor abundances are set to zero. For `LBQ` data, entries with high
#' isolation interference (>50%) are filtered out. In both cases, only proteins
#' identified in a single protein group are retained. Moreover, contaminants are
#' removed from the data.
#'
#' @param table Character. Path to the tab-delimited file containing the quantification data.
#' @param contaminants Character. Path to the tab-delimited file containing contaminants.
#'        Download `.tsv` file from UniProt.
#' @param type Character. Type of quantification data. Must be either `"LFQ"` (Label-Free Quantification)
#'        or `"LBQ"` (Label-Based Quantification).
#'
#' @return A data.frame with the filtered and preprocessed quantification data.
load_data <- function(table, contaminants, type = c("LFQ", "LBQ")) {
  type <- match.arg(type)
  
  data <- read.table(table, sep = "\t", header = TRUE)
  contaminants <- read.table(contaminants, sep = "\t", header = TRUE) %>% pull(Entry)
  
  cat(sprintf("Initial PSM count: %s\n", nrow(data)))
  
  if (type == "LFQ") {
    data <- data %>%
      mutate(Precursor.Abundance = if_else(is.na(Precursor.Abundance), 0, Precursor.Abundance))
  }
  
  if (type == "LBQ") {
    data <- data %>%
      dplyr::filter(
        !is.na(Isolation.Interference.in.Percent) &
          Isolation.Interference.in.Percent < 50
      )
    
    cat(sprintf(
      "After filtering for an isolation interference < 50%%: %s\n",
      nrow(data)
    ))
  }
  
  data <- data %>%
    dplyr::filter(Number.of.Protein.Groups == 1)
  
  cat(sprintf(
    "After filtering for PSMs with only one matching protein group: %s\n",
    nrow(data)
  ))
  
  data <- data %>%
    dplyr::filter(!(Master.Protein.Accessions %in% contaminants))
  
  cat(sprintf("After removing contaminants: %s\n", nrow(data)))
  
  return(data)
}