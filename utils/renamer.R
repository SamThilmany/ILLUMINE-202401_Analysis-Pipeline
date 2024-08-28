# Rename condition slugs and names
condition_renamer <- function(condition, output = "name") {
  if (!output %in% c("name", "slug")) {
    stop(sprintf("%s is not a valid argument for `output` in the `condition_renamer` method.", output))
  }

  # Some conditions were named differently by the cooperation partners
  # this is changed here to be consistent in every analysis
  renamed_condition <- switch(condition,
    "ctrl" = "CTRL",
    "EtL" = "EE.LNG",
    "S23" = "S.23",
    condition
  )

  if (output == "name") {
    renamed_condition <- switch(condition,
      "EE.LNG" = "EE + LNG",
      "S.23" = "S-23",
      condition
    )
  }

  return(renamed_condition)
}

# Rename abundance ratio columns
abundance_ratio_column_renamer <- function(column_name) {
  regex <- "Abundance\\.Ratio\\.{1,2}(log2|Adj\\.+P\\.Value)\\.{3,4}(.+)\\.{5}(.+)\\.{1}"

  matches <- stringr::str_match(column_name, regex)

  if (is.na(matches[1]) || is.na(matches[2]) || is.na(matches[3]) || is.na(matches[4])) {
    return(column_name)
  }

  if (grepl("Adj\\.+P\\.Value", matches[2])) {
    replacement_string <- "pValue_%s_%s"
  } else if (grepl("log2", matches[2])) {
    replacement_string <- "log2_%s_%s"
  }

  new_name <- sprintf(
    replacement_string,
    condition_renamer(matches[3], output = "slug"),
    condition_renamer(matches[4], output = "slug")
  )

  return(new_name)
}

# Rename abundances grouped cv columns
abundance_group_column_renamer <- function(column_name) {
  regex <- "Abundances\\.{2}Grouped\\.{2}(CV)\\.{6}(.+)"

  matches <- stringr::str_match(column_name, regex)

  if (is.na(matches[1]) || is.na(matches[2]) || is.na(matches[3])) {
    return(column_name)
  }

  new_name <- sprintf(
    "cv_%s",
    condition_renamer(matches[3], output = "slug")
  )
}