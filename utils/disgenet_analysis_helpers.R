library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(stringr, quietly = TRUE)
library(ggplot2, quietly = TRUE)

source("../utils/settings.R")
source("../utils/renamer.R")

# Helper function to prepare a vector of regulated genes
prepare_gene_list <- function(ratio_df_consistent, ratio) {
  regulated_gene_list <- ratio_df_consistent %>%
    filter(Comparison == ratio) %>%
    arrange(desc(Abs.log2.Ratio)) %>%
    pull(Gene.Symbol) %>%
    na.omit(Gene.Symbol) %>%
    unique()

  regulated_gene_list <- unlist(strsplit(regulated_gene_list, ";\\s*"))

  return(regulated_gene_list)
}

# Helper function to perform DISGENET enrichment
get_disgenet_enrichment <- function(regulated_gene_list) {
  # The number of genes for a GDA analysis is limited to 100
  regulated_gene_list <- head(regulated_gene_list, n = 100)

  disgenet_enrichment <- disgenet2r::gene2disease(
    gene = regulated_gene_list,
    vocabulary = "HGNC",
    database = "CURATED",
    score = c(0.75, 1),
    warnings = FALSE
  )

  return(disgenet_enrichment)
}

# Helper function to clean the DISGENET enrichment data
clean_enrichment_data <- function(enrichment) {
  enrichment_df <- enrichment@qresult %>%
    as.data.frame() %>%
    select(gene_symbol, disease_name, diseaseUMLSCUI, diseaseClasses_MSH) %>%
    unnest(cols = diseaseClasses_MSH, keep_empty = TRUE)

  return(enrichment_df)
}

# Helper function to shorten disease classes names
shorten_disease_name <- function(name, max_length = 30) {
  # Extract the MeSH code at the end (e.g., (C15))
  mesh_code <- str_extract(name, "\\([A-Z]\\d+\\)")
  # Remove the MeSH code from the name
  name_without_mesh_code <- str_remove(name, "\\([A-Z]\\d+\\)$")

  # Check if the length of the name without the MeSH code is greater than max_length
  shortened_name <- ifelse(
    str_length(name_without_mesh_code) > (max_length - str_length(mesh_code) - 3),
    str_sub(name_without_mesh_code, 1, max_length - str_length(mesh_code) - 3),
    name_without_mesh_code
  )

  # Add '...' and the MeSH code back if it was shortened
  result <- ifelse(
    str_length(name_without_mesh_code) > (max_length - str_length(mesh_code) - 3),
    str_trim(paste0(shortened_name, "... ", mesh_code)),
    paste0(name_without_mesh_code, " ", mesh_code)
  )

  return(result)
}

# Helper function to plot the a bar plot of the disease classes
plot_disease_classes <- function(enrichment_df, treatment, control) {
  df <- as.data.frame(table(enrichment_df$diseaseClasses_MSH), stringsAsFactors = FALSE) %>%
    setNames(c("Class.Name", "Freq")) %>%
    mutate(
      Shortened.Class.Name = shorten_disease_name(Class.Name, 20),
      Color = ifelse(str_detect(Class.Name, "\\(F03\\)$"), "highlight", "normal") # MeSH F04: Mental Disorders
    ) %>%
    drop_na()

  df$Shortened.Class.Name <- factor(df$Shortened.Class.Name, levels = df$Shortened.Class.Name[order(df$Freq)])

  plot <- ggplot(df, aes(x = Shortened.Class.Name, y = Freq, fill = Color)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
      title = sprintf("Disease Classes (%s vs. %s)", condition_renamer(treatment), condition_renamer(control)),
      x = "",
      y = "Number of GDAs"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  return(plot)
}