lfq_raw_paths <- list(
  `1` = "../data/LFQ/raw/sample_1.raw",
  `2` = "../data/LFQ/raw/sample_2.raw",
  `3` = "../data/LFQ/raw/sample_3.raw",
  `4` = "../data/LFQ/raw/sample_4.raw",
  `5` = "../data/LFQ/raw/sample_5.raw",
  `6` = "../data/LFQ/raw/sample_6.raw",
  `7` = "../data/LFQ/raw/sample_7.raw",
  `8` = "../data/LFQ/raw/sample_8.raw",
  `9` = "../data/LFQ/raw/sample_9.raw",
  `10` = "../data/LFQ/raw/sample_10.raw",
  `11` = "../data/LFQ/raw/sample_11.raw",
  `12` = "../data/LFQ/raw/sample_12.raw",
  `13` = "../data/LFQ/raw/sample_13.raw",
  `14` = "../data/LFQ/raw/sample_14.raw",
  `15` = "../data/LFQ/raw/sample_15.raw",
  `16` = "../data/LFQ/raw/sample_16.raw"
)

lfq_chimerys_proteins_path <- "../data/LFQ/Chimerys/protein.txt"
lfq_sequest_ht_proteins_path <- "../data/LFQ/SequestHT/protein.txt"

lfq_chimerys_peptides_path <- "../data/LFQ/Chimerys/peptide.txt"
lfq_sequest_ht_peptides_path <- "../data/LFQ/SequestHT/peptide.txt"

lfq_chimerys_psms_path <- "../data/LFQ/Chimerys/psm.txt"
lfq_sequest_ht_psms_path <- "../data/LFQ/SequestHT/psm.txt"

lfq_chimerys_ms2_path <- "../data/LFQ/Chimerys/ms2.txt"
lfq_sequest_ht_ms2_path <- "../data/LFQ/SequestHT/ms2.txt"

lbq_raw_paths <- list(
  `1.1` = "../data/LBQ/raw/sample_1_1.raw",
  `2.1` = "../data/LBQ/raw/sample_2_1.raw",
  `3.1` = "../data/LBQ/raw/sample_3_1.raw",
  `1.2` = "../data/LBQ/raw/sample_1_2.raw",
  `2.2` = "../data/LBQ/raw/sample_2_2.raw",
  `3.2` = "../data/LBQ/raw/sample_3_2.raw",
  `1.3` = "../data/LBQ/raw/sample_1_3.raw",
  `2.3` = "../data/LBQ/raw/sample_2_3.raw",
  `3.3` = "../data/LBQ/raw/sample_3_3.raw",
  `1.4` = "../data/LBQ/raw/sample_1_4.raw",
  `2.4` = "../data/LBQ/raw/sample_2_4.raw",
  `3.4` = "../data/LBQ/raw/sample_3_4.raw",
  `1.5` = "../data/LBQ/raw/sample_1_5.raw",
  `2.5` = "../data/LBQ/raw/sample_2_5.raw",
  `3.5` = "../data/LBQ/raw/sample_3_5.raw",
  `1.6` = "../data/LBQ/raw/sample_1_6.raw",
  `2.6` = "../data/LBQ/raw/sample_2_6.raw",
  `3.6` = "../data/LBQ/raw/sample_3_6.raw",
  `1.7` = "../data/LBQ/raw/sample_1_7.raw",
  `2.7` = "../data/LBQ/raw/sample_2_7.raw",
  `3.7` = "../data/LBQ/raw/sample_3_7.raw",
  `1.8` = "../data/LBQ/raw/sample_1_8.raw",
  `2.8` = "../data/LBQ/raw/sample_2_8.raw",
  `3.8` = "../data/LBQ/raw/sample_3_8.raw"
)

lbq_chimerys_proteins_path <- "../data/LBQ/Chimerys/protein.txt"
lbq_sequest_ht_proteins_path <- "../data/LBQ/SequestHT/protein.txt"

lbq_chimerys_peptides_path <- "../data/LBQ/Chimerys/peptide.txt"
lbq_sequest_ht_peptides_path <- "../data/LBQ/SequestHT/peptide.txt"

lbq_chimerys_psms_path <- "../data/LBQ/Chimerys/psm.txt"
lbq_sequest_ht_psms_path <- "../data/LBQ/SequestHT/psm.txt"

lbq_chimerys_ms2_path <- "../data/LBQ/Chimerys/ms2.txt"
lbq_sequest_ht_ms2_path <- "../data/LBQ/SequestHT/ms2.txt"

data <- new("DataLoader",
  lfq_raw_data_paths = lfq_raw_paths,

  lfq_chimerys_proteins = lfq_chimerys_proteins_path,
  lfq_sequest_ht_proteins = lfq_sequest_ht_proteins_path,

  lfq_chimerys_peptides = lfq_chimerys_peptides_path,
  lfq_sequest_ht_peptides = lfq_sequest_ht_peptides_path,

  lfq_chimerys_psms = lfq_chimerys_psms_path,
  lfq_sequest_ht_psms = lfq_sequest_ht_psms_path,

  lfq_chimerys_ms2 = lfq_chimerys_ms2_path,
  lfq_sequest_ht_ms2 = lfq_sequest_ht_ms2_path,

  lbq_raw_data_paths = lbq_raw_paths,

  lbq_chimerys_proteins = lbq_chimerys_proteins_path,
  lbq_sequest_ht_proteins = lbq_sequest_ht_proteins_path,

  lbq_chimerys_peptides = lbq_chimerys_peptides_path,
  lbq_sequest_ht_peptides = lbq_sequest_ht_peptides_path,

  lbq_chimerys_psms = lbq_chimerys_psms_path,
  lbq_sequest_ht_psms = lbq_sequest_ht_psms_path,

  lbq_chimerys_ms2 = lbq_chimerys_ms2_path,
  lbq_sequest_ht_ms2 = lbq_sequest_ht_ms2_path
)

data <- preprocessData(data)