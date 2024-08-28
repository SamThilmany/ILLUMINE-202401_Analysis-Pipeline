library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(stringr, quietly = TRUE)

source("../utils/renamer.R")

# Preprocess the protein data frame
preprocess_protein_df <- function(df) {
  df <- df %>%
    select(
      matches("^Accession"),
      matches("^Gene"),
      matches("^Entrez"),
      matches("Grouped..CV"),
      matches("Ratio..log2"),
      matches("Ratio.Adj..P.Value")
    ) %>%
    setNames(c(
      "Accession",
      "Gene.Symbol",
      "Entrez.Gene.ID",
      sapply(colnames(.)[grep("Group", colnames(.))], abundance_group_column_renamer),
      sapply(colnames(.)[grep("Ratio", colnames(.))], abundance_ratio_column_renamer)
    )) %>%
    filter(
      !stringr::str_detect(Gene.Symbol, "KRT")
    )

  return(df)
}

# Preprocess the peptide data frame
preprocess_peptide_df <- function(df) {
  return(df)
}

# Preprocess the psm data frame
preprocess_psm_df <- function(df, method, add_tech_repl = FALSE, add_fraction = FALSE) {
  if (!method %in% c("lfq", "lbq")) {
    stop(sprintf("%s is not a valid argument for `method` in the `preprocess_psm_df` method.", method))
  }

  if (!is.logical(add_tech_repl) || length(add_tech_repl) != 1) {
    stop(sprintf("%s is not a valid argument for `add_tech_repl` in the `preprocess_psm_df` method.", add_tech_repl))
  }

  if (!is.logical(add_fraction) || length(add_fraction) != 1) {
    stop(sprintf("%s is not a valid argument for `add_fraction` in the `preprocess_psm_df` method.", add_fraction))
  }

  # Define initial column patterns and new names
  column_patterns <- c(
    "^File",
    "^Charge",
    "^m.z",
    "^DeltaM",
    "^Isolation",
    "^Ion",
    "^Modifications"
  )
  new_names <- c(
    "File.ID",
    "Charge",
    "mz",
    "deltaM",
    "Isolation.Interference",
    "Ion.Inject.Time",
    "Modifications"
  )

  # Extend patterns and new names based on method
  if (method == "lfq") {
    column_patterns <- c(column_patterns, "^Precursor")
    new_names <- c(new_names, "Precursor.Abundance")
  } else if (method == "lbq") {
    column_patterns <- c(column_patterns, "^Abundance")
  }

  # Select matched columns
  df <- df %>%
    select(matches(column_patterns))

  # Check the presence of the columns and set new names accordingly
  present_columns <- colnames(df)
  if (method == "lbq") {
    new_names <- c(
      new_names,
      gsub("\\.\\.", ".", present_columns[-c(1:length(new_names))])
    )
  }

  # Rename the columns
  df <- setNames(df, new_names)

  # Continue with the pipeline
  df <- df %>%
    mutate(
      File.ID = sub("F", "", File.ID),
      Main.ID = as.integer(sub("\\..*", "", File.ID)),
      Sub.ID = ifelse(is.na(sub("^[^.]*\\.", "", File.ID)), "", sub("^[^.]*\\.", "", File.ID))
    ) %>%
    arrange(Main.ID, as.numeric(Sub.ID)) %>%
    mutate(File.ID = factor(File.ID, levels = unique(File.ID))) %>%
    select(-Main.ID, -Sub.ID)

  if (add_tech_repl) {
    df <- df %>%
      mutate(Tech.Replicate = as.numeric(sub("(\\d)\\.\\d", "\\1", File.ID)))
  }

  if (add_fraction) {
    df <- df %>%
      mutate(Fraction = as.numeric(sub("\\d\\.(\\d)", "\\1", File.ID)))
  }

  return(df)
}

# Preprocess the ms2 data frame
preprocess_ms2_df <- function(df, add_tech_repl = FALSE, add_fraction = FALSE) {
  if (!is.logical(add_tech_repl) || length(add_tech_repl) != 1) {
    stop(sprintf("%s is not a valid argument for `add_tech_repl` in the `preprocess_ms2_df` method.", add_tech_repl))
  }

  if (!is.logical(add_fraction) || length(add_fraction) != 1) {
    stop(sprintf("%s is not a valid argument for `add_fraction` in the `preprocess_ms2_df` method.", add_fraction))
  }

  df <- df %>%
    select(
      matches("File"),
      matches("Scan"),
      matches("PSM"),
      matches("^Ion")
    ) %>%
    setNames(c(
      "File.ID",
      "Scan",
      "PSMs",
      "Ion.Inject.Time"
    )) %>%
    mutate(
      File.ID = sub("F", "", File.ID),
      Main.ID = as.integer(sub("\\..*", "", File.ID)),
      Sub.ID = ifelse(is.na(sub("^[^.]*\\.", "", File.ID)), "", sub("^[^.]*\\.", "", File.ID))
    ) %>%
    arrange(Main.ID, as.numeric(Sub.ID)) %>%
    mutate(File.ID = factor(File.ID, levels = unique(File.ID))) %>%
    select(-Main.ID, -Sub.ID)

  if (add_tech_repl) {
    df <- df %>%
      mutate(Tech.Replicate = as.numeric(sub("(\\d)\\.\\d", "\\1", File.ID)))
  }

  if (add_fraction) {
    df <- df %>%
      mutate(Fraction = as.numeric(sub("\\d\\.(\\d)", "\\1", File.ID)))
  }

  return(df)
}

# Helper function to filter a protein data frame for features quantified in every sample group
filter_for_quantified_features <- function(protein_df) {
  protein_df <- protein_df %>%
    filter(complete.cases(select(., matches("log2|pValue"))))

  return(protein_df)
}