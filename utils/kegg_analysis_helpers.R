library(AnnotationDbi, quietly = TRUE)
library(org.Hs.eg.db, quietly = TRUE)
library(clusterProfiler, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(stringr, quietly = TRUE)

source("../utils/helpers.R")
source("../utils/renamer.R")

# Helper function to prepare a vector of regulated genes
prepare_gene_list <- function(ratio_df_consistent, ratio, regulation) {
  regulated_gene_list <- ratio_df_consistent %>%
    filter(Comparison == ratio) %>%
    filter(if (tolower(regulation) == "both") Diff.Expressed != "not significant" else Diff.Expressed == tolower(regulation)) %>%
    pull(Entrez.Gene.ID) %>%
    na.omit(Entrez.Gene.ID) %>%
    unique()

  regulated_gene_list <- unlist(strsplit(regulated_gene_list, ";\\s*"))

  return(regulated_gene_list)
}

# Helper function to perform GO analysis
get_kegg_enrichment <- function(regulated_gene_list, background_gene_list) {
  enrichment <- enrichKEGG(
    universe = background_gene_list,
    gene = regulated_gene_list,
    organism = "hsa",
    keyType = "kegg",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )

  return(enrichment)
}

# Helper function to draw GO plots
plot_kegg <- function(enrichment, directory = "./", regulation, treatment, control, node_name, node_slug) {
  regulation_label <- if (regulation == "both") {
    "Up- and down-regulated"
  } else {
    sprintf("%s-regulated", str_to_title(regulation))
  }

  create_plot <- function(plot_function, plot_name_prefix, ...) {
    plot <- plot_function(enrichment, ...) +
      labs(
        title = sprintf(
          "KEGG Enrichment\n(%s vs. %s; %s; %s)",
          condition_renamer(treatment), condition_renamer(control), regulation_label, node_name
        )
      ) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = default_font_size),
        axis.text.y = element_text(size = default_font_size),
        plot.title = element_text(hjust = 0.5, size = default_font_size * 14 / 11),
        legend.title = element_text(size = default_font_size),
        text = element_text(size = default_font_size)
      )

    save_and_show_plot(sprintf("kegg_%s_%s-vs-%s_%s_%s", plot_name_prefix, treatment, control, regulation, node_slug), directory, plot, show_width = 750, plot_width = 1.5 * default_width, plot_height = 2 * default_height)
  }

  create_plot(dotplot, "dotplot", showCategory = 10, font.size = default_font_size)

  create_plot(cnetplot, "cnetplot", showCategory = 10, font.size = default_font_size, cex.params = list(gene_label = 0.5, category_label = 0.75))
}