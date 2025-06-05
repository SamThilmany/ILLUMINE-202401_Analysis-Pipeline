#' Retrieve and format top differentially abundant proteins for a specific condition
#'
#' Filters and formats the top 100 proteins (by adjusted p-value) for a selected condition label.
#' Merges protein-level statistical data with associated metadata, and formats key numeric columns
#' for clean table presentation.
#'
#' @param data A data.frame containing statistical results including columns `Protein`, `LabelFactor`,
#'        `log2FC`, and `adj.pvalue`.
#' @param protein_metadata A data.frame containing additional metadata for each protein, to be
#'        merged by `Protein`. (provided by Proteome Discoverer)
#' @param label Character. The comparison label to filter by. Must be one of: `"EE_vs_DMSO"`,
#'        `"LNG_vs_DMSO"`, `"EE+LNG_vs_DMSO"`, or `"S-23_vs_DMSO"`.
#'
#' @return A data.frame with the top 100 proteins for the specified label, including formatted
#'         log2 fold changes, adjusted p-values, q-values, XCorr scores, and peptide information.
get_formatted_top_proteins_per_label <- function(data,
                                                 protein_metadata,
                                                 label = c(
                                                   "EE_vs_DMSO",
                                                   "LNG_vs_DMSO",
                                                   "EE+LNG_vs_DMSO",
                                                   "S-23_vs_DMSO"
                                                 )) {
  label <- match.arg(label)
  
  formatted_data <- data %>%
    select(Protein, LabelFactor, log2FC, adj.pvalue) %>%
    left_join(protein_metadata, by = "Protein") %>%
    dplyr::filter(LabelFactor == label) %>%
    select(-LabelFactor) %>%
    arrange(adj.pvalue) %>%
    slice_head(n = 200) %>%
    mutate(
      log2FC = format(round(log2FC, 4), nsmall = 4),
      adj.pvalue = formatC(adj.pvalue, format = "e", digits = 2),
      q_value = formatC(q_value, format = "e", digits = 2),
      XCorr = format(round(XCorr, 2), nsmall = 2)
    ) %>%
    select(
      Protein,
      Gene,
      Coverage,
      PSMs,
      Peptides,
      Unique_Peptides,
      q_value,
      XCorr,
      log2FC,
      adj.pvalue
    ) %>%
    rename(
      "log2 FC" = log2FC,
      "adj. p-value" = adj.pvalue,
      "q-value" = q_value,
      "Unique Peptides" = Unique_Peptides
    )
}

#' Map protein identifiers to Entrez Gene IDs
#'
#' Converts a list of UniProt protein identifiers into Entrez Gene IDs using the org.Hs.eg.db database.
#'
#' @param data A data.frame containing a `Protein` column or a character vector of protein IDs.
#'
#' @return A character vector of Entrez Gene IDs.
get_genes_from_proteins <- function(data) {
  if (!is.character(data)) {
    proteins <- data %>%
      pull(Protein) %>%
      unique()
  } else {
    proteins <- data
  }
  
  clusterProfiler::bitr(
    proteins,
    fromType = "UNIPROT",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  ) %>% pull(ENTREZID)
}

#' Get the list of intersecting genes from specific condition comparisons
#'
#' Identifies proteins that are significantly regulated (based on a p-value cutoff) in specific combinations
#' of experimental conditions and converts them to Entrez Gene IDs.
#'
#' @param data A data.frame containing differential abundance results with `Protein`, `Label`, and `adj.pvalue`.
#' @param fdr_threshold Numeric threshold for adjusted p-value significance.
#' @param region Character string specifying the region name (e.g., "EE", "EE+LNG", "EE/LNG", etc.).
#'
#' @return A character vector of Entrez Gene IDs from the selected region.
get_intersection_genes <- function(data,
                                   fdr_threshold,
                                   region = c(
                                     "EE",
                                     "LNG",
                                     "EE+LNG",
                                     "S-23",
                                     "EE/LNG",
                                     "EE/EE+LNG",
                                     "EE/S-23",
                                     "LNG/EE+LNG",
                                     "LNG/S-23",
                                     "EE+LNG/S-23",
                                     "EE/LNG/EE+LNG",
                                     "EE/LNG/S-23",
                                     "EE/EE+LNG/S-23",
                                     "LNG/EE+LNG/S-23",
                                     "EE/LNG/EE+LNG/S-23"
                                   )) {
  region <- match.arg(region)
  
  significant <- data %>%
    dplyr::filter(adj.pvalue < fdr_threshold)
  
  set1 <- significant %>%
    dplyr::filter(Label == "EE_vs_DMSO") %>%
    pull(Protein) %>%
    unique() %>%
    as.character()
  set2 <- significant %>%
    dplyr::filter(Label == "LNG_vs_DMSO") %>%
    pull(Protein) %>%
    unique() %>%
    as.character()
  set3 <- significant %>%
    dplyr::filter(Label == "EE+LNG_vs_DMSO") %>%
    pull(Protein) %>%
    unique() %>%
    as.character()
  set4 <- significant %>%
    dplyr::filter(Label == "S-23_vs_DMSO") %>%
    pull(Protein) %>%
    unique() %>%
    as.character()
  
  dap_significant_sets <- list(
    "EE" = set1,
    "LNG" = set2,
    "EE+LNG" = set3,
    "S-23" = set4
  )
  
  dap_significant_sets_venn_obj <- Venn(dap_significant_sets)
  dap_significant_sets_region_data <- process_data(dap_significant_sets_venn_obj)
  
  intersection_proteins <- ggVennDiagram:::get_shape_regionlabel(dap_significant_sets_region_data) %>%
    dplyr::filter(name == region) %>%
    pull(item) %>%
    unlist() %>%
    unique()
  
  get_genes_from_proteins(intersection_proteins)
}

#' Perform Gene Ontology (GO) enrichment analysis
#'
#' Conducts GO enrichment analysis using a list of genes and a background universe,
#' and returns a formatted summary table.
#'
#' @param gene_list A character vector of Entrez Gene IDs to test.
#' @param background_list A character vector of Entrez Gene IDs as background.
#' @param ont Character string specifying which ontology to use: "BP", "CC", or "MF".
#'
#' @return A data.frame summarizing the GO enrichment results.
get_go_enrichment <- function(gene_list,
                              background_list,
                              ont = c("BP", "CC", "MF")) {
  ont <- match.arg(ont)
  
  go_enrichment <- clusterProfiler::enrichGO(
    gene = gene_list,
    universe = background_list,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  ) %>% clusterProfiler::simplify(cutoff = 0.7,
                                  by = "p.adjust",
                                  select_fun = min)
  
  go_enrichment <- go_enrichment %>%
    as.data.frame() %>%
    dplyr::select(
      "ID",
      "Description",
      "GeneRatio",
      "BgRatio",
      "FoldEnrichment",
      "Count",
      "p.adjust",
      "qvalue"
    ) %>%
    mutate(Label = sprintf(
      "%s (%s)",
      if_else(
        str_length(Description) > default_label_characters,
        paste0(str_sub(
          Description, 1, default_label_characters
        ), "..."),
        Description
      ),
      ID
    )) %>%
    unique() %>%
    arrange(p.adjust) %>%
    rowwise() %>%
    mutate(
      log10p_adjusted = -log10(p.adjust),
      calcGeneRatio = eval(parse(text = GeneRatio)),
      calcBgRatio = eval(parse(text = BgRatio))
    )
  
  return(go_enrichment)
}

#' Perform Kyoto Encyclopedia of Genes and Genomes (KEGG) enrichment analysis
#'
#' Conducts KEGG enrichment analysis using a list of genes and a background universe,
#' and returns a formatted summary table.
#'
#' @param gene_list A character vector of Entrez Gene IDs to test.
#' @param background_list A character vector of Entrez Gene IDs as background.
#'
#' @return A data.frame summarizing the KEGG enrichment results.
get_kegg_enrichment <- function(gene_list, background_list) {
  kegg_enrichment <- clusterProfiler::enrichKEGG(
    gene = gene_list,
    universe = background_list,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  kegg_enrichment <- kegg_enrichment %>%
    as.data.frame() %>%
    dplyr::select("ID",
                  "Description",
                  "GeneRatio",
                  "BgRatio",
                  "FoldEnrichment",
                  "Count",
                  "p.adjust") %>%
    mutate(Label = sprintf(
      "%s (%s)",
      if_else(
        str_length(Description) > default_label_characters,
        paste0(str_sub(
          Description, 1, default_label_characters
        ), "..."),
        Description
      ),
      ID
    )) %>%
    unique() %>%
    arrange(p.adjust) %>%
    rowwise() %>%
    mutate(
      log10p_adjusted = -log10(p.adjust),
      calcGeneRatio = eval(parse(text = GeneRatio)),
      calcBgRatio = eval(parse(text = BgRatio))
    )
  
  return(kegg_enrichment)
}

#' Perform Disease Ontology Semantic and Enrichment analysis (DOSE)
#'
#' Conducts DOSE enrichment analysis using a list of genes and a background universe,
#' and returns a formatted summary table.
#'
#' @param gene_list A character vector of Entrez Gene IDs to test.
#' @param background_list A character vector of Entrez Gene IDs as background.
#'
#' @return A data.frame summarizing the DOSE enrichment results.
get_dose_enrichment <- function(gene_list, background_list) {
  do_enrichment <- clusterProfiler::enrichDO(
    gene = gene_list,
    universe = background_list,
    ont = "HDO",
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  do_enrichment <- do_enrichment %>%
    as.data.frame() %>%
    dplyr::select("ID",
                  "Description",
                  "GeneRatio",
                  "BgRatio",
                  "FoldEnrichment",
                  "Count",
                  "p.adjust") %>%
    mutate(Label = sprintf(
      "%s (%s)",
      if_else(
        str_length(Description) > default_label_characters,
        paste0(str_sub(
          Description, 1, default_label_characters
        ), "..."),
        Description
      ),
      ID
    )) %>%
    unique() %>%
    arrange(p.adjust) %>%
    rowwise() %>%
    mutate(
      log10p_adjusted = -log10(p.adjust),
      calcGeneRatio = eval(parse(text = GeneRatio)),
      calcBgRatio = eval(parse(text = BgRatio))
    )
  
  return(do_enrichment)
}

#' Perform DisGeNET enrichment analysis
#'
#' Conducts DisGeNET enrichment analysis using a list of genes and a background universe,
#' and returns a formatted summary table.
#'
#' @param gene_list A character vector of Entrez Gene IDs to test.
#' @param background_list A character vector of Entrez Gene IDs as background.
#'
#' @return A data.frame summarizing the DISGENET enrichment results.
get_disgenet_enrichment <- function(gene_list, background_list) {
  dgn_enrichment <- DOSE::enrichDGN(
    gene = gene_list,
    universe = background_list,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  dgn_enrichment <- dgn_enrichment %>%
    as.data.frame() %>%
    dplyr::select("ID",
                  "Description",
                  "GeneRatio",
                  "BgRatio",
                  "FoldEnrichment",
                  "Count",
                  "p.adjust") %>%
    mutate(Label = sprintf(
      "%s (%s)",
      if_else(
        str_length(Description) > default_label_characters,
        paste0(str_sub(
          Description, 1, default_label_characters
        ), "..."),
        Description
      ),
      ID
    )) %>%
    unique() %>%
    arrange(p.adjust) %>%
    rowwise() %>%
    mutate(
      log10p_adjusted = -log10(p.adjust),
      calcGeneRatio = eval(parse(text = GeneRatio)),
      calcBgRatio = eval(parse(text = BgRatio))
    )
  
  return(dgn_enrichment)
}

#' Create a table from a GO, KEGG, or DISGENET enrichment result
#'
#' Formats and summarizes the results of an enrichment analysis by selecting relevant columns.
#'
#' @param enrichmentResult A data.frame produced by `get_go_enrichment()`, `get_kegg_enrichment()`,
#'        or `get_disgenet_enrichment()`.
#'
#' @return A summarized and sorted data.frame.
get_enrichment_table <- function(enrichmentResult) {
  enrichmentResult %>%
    as.data.frame() %>%
    dplyr::select(
      "ID",
      "Description",
      "Label",
      "GeneRatio",
      "calcGeneRatio",
      "BgRatio",
      "FoldEnrichment",
      "Count",
      "p.adjust"
    ) %>%
    unique() %>%
    arrange(p.adjust)
}

#' Plot the results of a GO, KEGG, or DISGENET enrichment analysis
#'
#' Creates a dot plot showing the top enriched terms, colored by p-value significance and sized by gene count.
#'
#' @param enrichmentResult A data.frame produced by `get_go_enrichment()`, `get_kegg_enrichment()`,
#'        or `get_disgenet_enrichment().
#' @param title A character string for the plot subtitle.
#'
#' @return A `ggplot` object showing the enrichment results.
get_enrichment_plot <- function(enrichmentResult, title) {
  enrichmentTable <- enrichmentResult %>%
    as.data.frame() %>%
    dplyr::select(
      "Label",
      "calcGeneRatio",
      "calcBgRatio",
      "FoldEnrichment",
      "Count",
      "log10p_adjusted"
    ) %>%
    arrange(desc(log10p_adjusted)) %>%
    slice_head(n = 15)
  
  ggplot(
    enrichmentTable,
    aes(
      x = calcGeneRatio,
      y = reorder(Label, calcGeneRatio),
      size = Count,
      color = log10p_adjusted
    )
  ) +
    geom_point() +
    scale_color_gradient(low = hue_pal()(4)[3], high = hue_pal()(4)[1]) +
    labs(
      subtitle = title,
      x = "Gene Ratio",
      y = NULL,
      size = "Gene Count",
      color = expression("-log"[10] * " (adj. p-value)")
    ) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      axis.ticks = element_line(color = "black"),
      plot.title.position = "panel",
      plot.subtitle = element_text(
        size = 0.8 * default_font_size,
        color = default_text_color,
        hjust = 0.5,
        margin = margin(
          t = 0,
          r = 0,
          b = 5,
          l = 0
        )
      ),
      legend.box = "vertical"
    )
}

#' Perform and visualize enrichment analysis for a set of genes
#'
#' Executes one of several supported enrichment analyses (GO, KEGG, DOSE, or DISGENET)
#' for a given list of genes. Optionally allows specifying a background gene list for
#' analyses that support it. Outputs both a tabular summary and a dot plot visualization
#' of enrichment results, saving them to disk with filenames derived from a given slug.
#'
#' @param type The type of enrichment analysis to perform. One of: "GO", "KEGG", "DOSE", or "DISGENET".
#' @param genes_of_interest A character vector of gene identifiers (ENTREZ IDs) to be analyzed.
#' @param background_genes An optional character vector of background genes used for enrichment calculation.
#'        Required for GO, KEGG, and DOSE enrichment.
#' @param name A string used as the title in the resulting dot plot visualization.
#' @param slug A short identifier used to generate filenames for saving the output table and plot.
#'
#' @return A data.frame.
perform_enrichment_analysis <- function(type = c("GO", "KEGG", "DOSE", "DISGENET"),
                                        genes_of_interest,
                                        background_genes = NULL,
                                        name,
                                        slug) {
  type <- match.arg(type)
  
  if (type == "GO") {
    enrichment <- get_go_enrichment(genes_of_interest, background_genes, ont = "BP")
  } else if (type == "KEGG") {
    enrichment <- get_kegg_enrichment(genes_of_interest, background_genes)
  } else if (type == "DOSE") {
    enrichment <- get_dose_enrichment(genes_of_interest, background_genes)
  } else if (type == "DISGENET") {
    enrichment <- get_disgenet_enrichment(genes_of_interest)
  }
  
  if (nrow(enrichment %>% as.data.frame()) == 0) {
    return(NULL)
  }
  
  enrichment_tab <- get_enrichment_table(enrichment)
  save_and_show_table(
    sprintf("%s_%s_tab", tolower(type), slug),
    enrichment_tab %>% select(-Label, -calcGeneRatio)
  )
  
  enrichment_dotplot <- get_enrichment_plot(enrichment, title = name)
  save_and_show_plot(sprintf("%s_%s_dotplot", tolower(type), slug),
                     enrichment_dotplot,
                     plot_width = 120)
  
  cat(sprintf(
    "The analysis resulted in %s enriched entries.",
    nrow(enrichment %>% as.data.frame())
  ))
  
  return(enrichment_tab %>% mutate(Name = name, Slug = slug))
}

#' Generate a bubble plot of enriched terms across contrasts
#'
#' This function combines multiple enrichment result data frames into a single
#' dataset and visualizes the top enriched terms (by frequency and adjusted p-value)
#' across all provided contrasts. Optionally, specific terms can be highlighted in the heatmap.
#'
#' @param results A list of enrichment result data frames, each containing at least
#'        the columns: `Name`, `Description`, `p.adjust`, and `Label`.
#'        `Name` should indicate the contrast or condition.
#' @param highlighted_terms Optional vector of `Description` strings to be highlighted
#'        in the heatmap (rendered in red using markdown styling).
#' @param top_n Integer indicating the number of top enriched terms to display, selected
#'        based on frequency of occurrence and mean adjusted p-value. Default is 30.
#'
#' @return A `ggplot` object representing a bubble plot of enriched terms, where:
#'         - y-axis shows enriched terms (with optional highlighting),
#'         - x-axis shows contrasts,
#'         - point color reflects adjusted p-value,
#'         - point size reflects gene ratio.
get_enrichment_bubbleplot <- function(results,
                                      highlighted_terms = NULL,
                                      top_n = 30) {
  data <- bind_rows(results)
  contrast_levels <- unique(data$Name)
  
  all_results <- data %>%
    group_by(Description) %>%
    mutate(
      occurences = n(),
      mean.p.adjust = mean(p.adjust),
      Name = factor(Name, levels = contrast_levels)
    ) %>%
    ungroup()
  
  top_terms <- all_results %>%
    distinct(Description, Label, occurences, mean.p.adjust) %>%
    arrange(desc(occurences), mean.p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(FormattedLabel = if_else(
      Description %in% highlighted_terms,
      paste0("<b style='color:darkred;'>", Label, "</b>"),
      Label
    ))
  
  filtered_df <- all_results %>%
    dplyr::filter(Description %in% top_terms$Description) %>%
    left_join(top_terms %>% select(Description, FormattedLabel), by = "Description") %>%
    mutate(FormattedLabel = factor(FormattedLabel, levels = rev(top_terms$FormattedLabel)))
  
  cat(unique(filtered_df$Description), sep = "\n")
  
  export <- filtered_df %>%
    select(c(ID, Description, GeneRatio, p.adjust, Name)) %>%
    pivot_wider(names_from = Name,
                values_from = c(GeneRatio, p.adjust))
  
  save_and_show_table("go_bubble", export)
  
  ggplot(filtered_df,
         aes(Name, FormattedLabel, color = p.adjust, size = calcGeneRatio)) +
    geom_point() +
    labs(
      x = "",
      y = "",
      size = "Gene Ratio",
      color = "adj. p-Value"
    ) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 6
      ),
      axis.text.y = ggtext::element_markdown(size = 6)
    )
}
