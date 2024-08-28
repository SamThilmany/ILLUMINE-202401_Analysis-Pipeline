library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(ggplot2, quietly = TRUE)

source("../utils/helpers.R")

# Helper function to get the number of identified and quantified proteins
get_ident_quant_protein_no <- function(protein_data_frame_a, protein_data_frame_a_complete, protein_data_frame_b, protein_data_frame_b_complete, a_name, b_name) {
  ident_quant_protein_no <- data.frame(
    "Node" = c(a_name, b_name),
    "Ident." = c(nrow(protein_data_frame_a), nrow(protein_data_frame_b)),
    "Quant." = c(nrow(protein_data_frame_a_complete), nrow(protein_data_frame_b_complete))
  ) %>%
    pivot_longer(
      cols = c("Ident.", "Quant."),
      names_to = "Type",
      values_to = "Count"
    )

  return(ident_quant_protein_no)
}

# Helper function to get the identification rate
get_psm_ms2_ident_rate <- function(ms2_data_frame_a, psm_data_frame_a, ms2_data_frame_b, psm_data_frame_b, a_name, b_name) {
  psm_ms2_ident_rate <- data.frame(
    "Node" = c(a_name, b_name),
    "PSM" = c(nrow(psm_data_frame_a), nrow(psm_data_frame_b)),
    "MS2" = c(nrow(ms2_data_frame_a), nrow(ms2_data_frame_b))
  ) %>%
    mutate("ID.Percentage" = PSM / MS2 * 100) %>%
    pivot_longer(cols = c(PSM, MS2)) %>%
    setNames(c("Node", "ID.Percentage", "Type", "Count")) %>%
    mutate(ID.Percentage = ifelse(Type == "MS2", 100, ID.Percentage)) %>%
    mutate(Label = ifelse(Type == "MS2", "", sprintf("%.2f %%", ID.Percentage)))

  return(psm_ms2_ident_rate)
}

# Helper function to plot features
plot_features <- function(df) {
  plot <- ggplot(df, aes(x = Type, y = Count, fill = Type)) +
    geom_bar(stat = "identity") +
    facet_grid(cols = vars(Node)) +
    geom_text(aes(label = format(Count, big.mark = ",")), angle = 90, vjust = 0.5, hjust = 1.25) +
    labs(
      title = "Identified and quantified proteins",
      x = "",
      y = "N° Proteins"
    ) +
    theme(
      legend.position = "none"
    )

  return(plot)
}

# Helper function to plot feature ID rate
plot_feature_id_rate <- function(df, plot_mode) {
  if (plot_mode == "absolute") {
    plot <- ggplot(df, aes(x = Type, y = Count, fill = Type))
  }

  if (plot_mode == "relative") {
    plot <- ggplot(df, aes(x = Type, y = ID.Percentage, fill = Type))
  }

  plot <- plot +
    geom_bar(stat = "identity") +
    scale_x_discrete(labels = c("MS2" = expression(MS^2), PSM = "PSM")) +
    facet_grid(cols = vars(Node)) +
    labs(
      title = expression(MS^2 ~ scan ~ identification),
      x = ""
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(hjust = 0.5, vjust = 0)
    )

  if (plot_mode == "absolute") {
    plot <- plot +
      scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
      geom_text(aes(label = format(Count, big.mark = ",")), angle = 90, vjust = 0.5, hjust = 1.25) +
      labs(
        y = "Scans / IDs"
      )
  }

  if (plot_mode == "relative") {
    plot <- plot +
      geom_hline(yintercept = 40, color = "Red", alpha = 0.5, linetype = "dashed", linewidth = 0.5) +
      annotate("text", 1.5, 40, label = "40 %\nID rate", color = "Red", alpha = 0.75) +
      geom_text(aes(label = Label), angle = 90, vjust = 0.5, hjust = 1.25) +
      labs(
        y = "ID rate"
      )
  }

  return(plot)
}