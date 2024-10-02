source("../utils/helpers.R")
source("../utils/data_preprocessing.R")

setClass(
  Class = "DataLoader",
  slots = c(
    lfqRawDataPaths = "list",

    lfqChimerys_proteins = "data.frame",
    lfqSequestHT_proteins = "data.frame",

    lfqChimerys_proteins_complete = "data.frame",
    lfqSequestHT_proteins_complete = "data.frame",

    lfqChimerys_peptides = "data.frame",
    lfqSequestHT_peptides = "data.frame",

    lfqChimerys_psms = "data.frame",
    lfqSequestHT_psms = "data.frame",

    lfqChimerys_ms2 = "data.frame",
    lfqSequestHT_ms2 = "data.frame",

    lbqRawDataPaths = "list",

    lbqChimerys_proteins = "data.frame",
    lbqSequestHT_proteins = "data.frame",

    lbqChimerys_proteins_complete = "data.frame",
    lbqSequestHT_proteins_complete = "data.frame",

    lbqChimerys_peptides = "data.frame",
    lbqSequestHT_peptides = "data.frame",

    lbqChimerys_psms = "data.frame",
    lbqSequestHT_psms = "data.frame",

    lbqChimerys_ms2 = "data.frame",
    lbqSequestHT_ms2 = "data.frame",

    geneDiseaseAssociations = "data.frame"
  ),
  prototype = list(
    lfqRawDataPaths = list(),

    lfqChimerys_proteins = data.frame(),
    lfqSequestHT_proteins = data.frame(),

    lfqChimerys_proteins_complete = data.frame(),
    lfqSequestHT_proteins_complete = data.frame(),

    lfqChimerys_peptides = data.frame(),
    lfqSequestHT_peptides = data.frame(),

    lfqChimerys_psms = data.frame(),
    lfqSequestHT_psms = data.frame(),

    lfqChimerys_ms2 = data.frame(),
    lfqSequestHT_ms2 = data.frame(),

    lbqRawDataPaths = list(),

    lbqChimerys_proteins = data.frame(),
    lbqSequestHT_proteins = data.frame(),

    lbqChimerys_proteins_complete = data.frame(),
    lbqSequestHT_proteins_complete = data.frame(),

    lbqChimerys_peptides = data.frame(),
    lbqSequestHT_peptides = data.frame(),

    lbqChimerys_psms = data.frame(),
    lbqSequestHT_psms = data.frame(),

    lbqChimerys_ms2 = data.frame(),
    lbqSequestHT_ms2 = data.frame(),

    geneDiseaseAssociations = data.frame()
  )
)

setMethod(
  f = "initialize",
  signature = "DataLoader",
  definition = function(.Object,
    lfq_raw_data_paths = list(),

    lfq_chimerys_proteins = NULL,
    lfq_sequest_ht_proteins = NULL,

    lfq_chimerys_peptides = NULL,
    lfq_sequest_ht_peptides = NULL,

    lfq_chimerys_psms = NULL,
    lfq_sequest_ht_psms = NULL,

    lfq_chimerys_ms2 = NULL,
    lfq_sequest_ht_ms2 = NULL,

    lbq_raw_data_paths = list(),

    lbq_chimerys_proteins = NULL,
    lbq_sequest_ht_proteins = NULL,

    lbq_chimerys_peptides = NULL,
    lbq_sequest_ht_peptides = NULL,

    lbq_chimerys_psms = NULL,
    lbq_sequest_ht_psms = NULL,

    lbq_chimerys_ms2 = NULL,
    lbq_sequest_ht_ms2 = NULL
  ) {
    .Object@lfqRawDataPaths <- lfq_raw_data_paths

    .Object@lfqChimerys_proteins <- if (!is.null(lfq_chimerys_proteins) && file.exists(lfq_chimerys_proteins)) read.table(lfq_chimerys_proteins, sep = "\t", header = TRUE) else data.frame()
    .Object@lfqSequestHT_proteins <- if (!is.null(lfq_sequest_ht_proteins) && file.exists(lfq_sequest_ht_proteins)) read.table(lfq_sequest_ht_proteins, sep = "\t", header = TRUE) else data.frame()

    .Object@lfqChimerys_peptides <- if (!is.null(lfq_chimerys_peptides) && file.exists(lfq_chimerys_peptides)) read.table(lfq_chimerys_peptides, sep = "\t", header = TRUE) else data.frame()
    .Object@lfqSequestHT_peptides <- if (!is.null(lfq_sequest_ht_peptides) && file.exists(lfq_sequest_ht_peptides)) read.table(lfq_sequest_ht_peptides, sep = "\t", header = TRUE) else data.frame()

    .Object@lfqChimerys_psms <- if (!is.null(lfq_chimerys_psms) && file.exists(lfq_chimerys_psms)) read.table(lfq_chimerys_psms, sep = "\t", header = TRUE) else data.frame()
    .Object@lfqSequestHT_psms <- if (!is.null(lfq_sequest_ht_psms) && file.exists(lfq_sequest_ht_psms)) read.table(lfq_sequest_ht_psms, sep = "\t", header = TRUE) else data.frame()

    .Object@lfqChimerys_ms2 <- if (!is.null(lfq_chimerys_ms2) && file.exists(lfq_chimerys_ms2)) read.table(lfq_chimerys_ms2, sep = "\t", header = TRUE) else data.frame()
    .Object@lfqSequestHT_ms2 <- if (!is.null(lfq_sequest_ht_ms2) && file.exists(lfq_sequest_ht_ms2)) read.table(lfq_sequest_ht_ms2, sep = "\t", header = TRUE) else data.frame()

    .Object@lbqRawDataPaths <- lbq_raw_data_paths

    .Object@lbqChimerys_proteins <- if (!is.null(lbq_chimerys_proteins) && file.exists(lbq_chimerys_proteins)) read.table(lbq_chimerys_proteins, sep = "\t", header = TRUE) else data.frame()
    .Object@lbqSequestHT_proteins <- if (!is.null(lbq_sequest_ht_proteins) && file.exists(lbq_sequest_ht_proteins)) read.table(lbq_sequest_ht_proteins, sep = "\t", header = TRUE) else data.frame()

    .Object@lbqChimerys_peptides <- if (!is.null(lbq_chimerys_peptides) && file.exists(lbq_chimerys_peptides)) read.table(lbq_chimerys_peptides, sep = "\t", header = TRUE) else data.frame()
    .Object@lbqSequestHT_peptides <- if (!is.null(lbq_sequest_ht_peptides) && file.exists(lbq_sequest_ht_peptides)) read.table(lbq_sequest_ht_peptides, sep = "\t", header = TRUE) else data.frame()

    .Object@lbqChimerys_psms <- if (!is.null(lbq_chimerys_psms) && file.exists(lbq_chimerys_psms)) read.table(lbq_chimerys_psms, sep = "\t", header = TRUE) else data.frame()
    .Object@lbqSequestHT_psms <- if (!is.null(lbq_sequest_ht_psms) && file.exists(lbq_sequest_ht_psms)) read.table(lbq_sequest_ht_psms, sep = "\t", header = TRUE) else data.frame()

    .Object@lbqChimerys_ms2 <- if (!is.null(lbq_chimerys_ms2) && file.exists(lbq_chimerys_ms2)) read.table(lbq_chimerys_ms2, sep = "\t", header = TRUE) else data.frame()
    .Object@lbqSequestHT_ms2 <- if (!is.null(lbq_sequest_ht_ms2) && file.exists(lbq_sequest_ht_ms2)) read.table(lbq_sequest_ht_ms2, sep = "\t", header = TRUE) else data.frame()

    .Object@geneDiseaseAssociations <- get_gene_disease_associations(diseases_of_interest)

    return(.Object)
  }
)

setGeneric(
  name = "preprocessData",
  def = function(.Object) {
    standardGeneric("preprocessData")
  }
)

setMethod(
  f = "preprocessData",
  signature = "DataLoader",
  definition = function(.Object) {
    .Object@lfqChimerys_proteins <- preprocess_protein_df(.Object@lfqChimerys_proteins)
    .Object@lfqSequestHT_proteins <- preprocess_protein_df(.Object@lfqSequestHT_proteins)

    .Object@lfqChimerys_proteins_complete <- filter_for_quantified_features(.Object@lfqChimerys_proteins, "LFQ Chimerys")
    .Object@lfqSequestHT_proteins_complete <- filter_for_quantified_features(.Object@lfqSequestHT_proteins, "LFQ SequestHT")

    .Object@lfqChimerys_peptides <- preprocess_peptide_df(.Object@lfqChimerys_peptides)
    .Object@lfqSequestHT_peptides <- preprocess_peptide_df(.Object@lfqSequestHT_peptides)

    .Object@lfqChimerys_psms <- preprocess_psm_df(.Object@lfqChimerys_psms, method = "lfq", add_tech_repl = FALSE, add_fraction = FALSE)
    .Object@lfqSequestHT_psms <- preprocess_psm_df(.Object@lfqSequestHT_psms, method = "lfq", add_tech_repl = FALSE, add_fraction = FALSE)

    .Object@lfqChimerys_ms2 <- preprocess_ms2_df(.Object@lfqChimerys_ms2, add_tech_repl = FALSE, add_fraction = FALSE)
    .Object@lfqSequestHT_ms2 <- preprocess_ms2_df(.Object@lfqSequestHT_ms2, add_tech_repl = FALSE, add_fraction = FALSE)

    .Object@lbqChimerys_proteins <- preprocess_protein_df(.Object@lbqChimerys_proteins)
    .Object@lbqSequestHT_proteins <- preprocess_protein_df(.Object@lbqSequestHT_proteins)

    .Object@lbqChimerys_proteins_complete <- filter_for_quantified_features(.Object@lbqChimerys_proteins, "LBQ Chimerys")
    .Object@lbqSequestHT_proteins_complete <- filter_for_quantified_features(.Object@lbqSequestHT_proteins, "LBQ SequestHT")

    .Object@lbqChimerys_peptides <- preprocess_peptide_df(.Object@lbqChimerys_peptides)
    .Object@lbqSequestHT_peptides <- preprocess_peptide_df(.Object@lbqSequestHT_peptides)

    .Object@lbqChimerys_psms <- preprocess_psm_df(.Object@lbqChimerys_psms, method = "lbq", add_tech_repl = TRUE, add_fraction = TRUE)
    .Object@lbqSequestHT_psms <- preprocess_psm_df(.Object@lbqSequestHT_psms, method = "lbq", add_tech_repl = TRUE, add_fraction = TRUE)

    .Object@lbqChimerys_ms2 <- preprocess_ms2_df(.Object@lbqChimerys_ms2, add_tech_repl = TRUE, add_fraction = TRUE)
    .Object@lbqSequestHT_ms2 <- preprocess_ms2_df(.Object@lbqSequestHT_ms2, add_tech_repl = TRUE, add_fraction = TRUE)

    return(.Object)
  }
)