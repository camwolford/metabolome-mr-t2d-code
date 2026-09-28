# Shared helpers for the analysis and figure scripts.

# Data locations are read from environment variables, e.g. METABOLOME_MR_INPUT_DIR.
data_paths <- function(required) {
  labels <- c(
    METABOLOME_MR_INPUT_DIR = "input_dir",
    METABOLOME_MR_FULL_METABOLITE_GWAS_DIR = "full_metabolite_gwas_dir",
    METABOLOME_MR_WORK_DIR = "work_dir",
    METABOLOME_MR_OUTPUT_DIR = "output_dir",
    METABOLOME_MR_EUR_PANEL_DIR = "eur_panel_dir"
  )
  stats::setNames(Sys.getenv(required), labels[required])
}

# Allele "T" can be read back from TSV as logical TRUE; restore it.
normalise_allele <- function(values) {
  values <- as.character(values)
  values[values %in% c("TRUE", "True", "true")] <- "T"
  toupper(values)
}

# Bonferroni threshold over the metabolite-outcome pairs tested (583 type 2 diabetes, 593 fasting glucose, 593 HbA1c).
significance_threshold <- 0.05 / (583 + 593 + 593)

harmonised_directory <- function(root, design, outcome) {
  suffix <- if (identical(design, "liberal")) "_Liberal" else ""
  file.path(root, "02_instrument_selection", design, paste0("Harmonised_", outcome, "_IVs", suffix))
}

# Map each metabolite to its harmonised instrument file in one directory.
instrument_files <- function(directory) {
  files <- list.files(directory, pattern = "[.]tsv$", full.names = TRUE)
  metabolites <- vapply(files, function(f) as.character(readr::read_tsv(f, n_max = 1, show_col_types = FALSE)$Metabolite), "")
  stats::setNames(files, metabolites)
}

# Complete metabolite GWAS (Surendran et al. 2022): one METAL file per Metabolon COMP_ID, e.g.
# INTERVAL_MRC-Epi_M17769.tbl.gz, with MarkerName chr:pos:A1:A2. Read in chunks, keeping rows where keep() is TRUE.
read_metabolite_gwas <- function(comp_id, keep) {
  path <- file.path(Sys.getenv("METABOLOME_MR_FULL_METABOLITE_GWAS_DIR"), sprintf("INTERVAL_MRC-Epi_%s.tbl.gz", comp_id))
  columns <- readr::cols_only(MarkerName = "c", Allele1 = "c", Allele2 = "c", Freq1 = "d", Effect = "d", StdErr = "d", `P-value` = "d", N = "d")
  standardise <- function(chunk, index) {
    d <- data.frame(
      Chromosome = suppressWarnings(as.numeric(sub("^chr([^:]+):.*$", "\\1", chunk$MarkerName))),  # chrX becomes NA
      Position = as.numeric(sub("^[^:]+:([^:]+):.*$", "\\1", chunk$MarkerName)),
      EffectAllele = toupper(chunk$Allele1), NonEffectAllele = toupper(chunk$Allele2),
      EAF = chunk$Freq1, Beta = chunk$Effect, SE = chunk$StdErr, Pval = chunk$`P-value`, N = chunk$N
    )
    d <- d[!is.na(d$Chromosome), , drop = FALSE]
    d[keep(d), , drop = FALSE]
  }
  as.data.frame(readr::read_tsv_chunked(path, readr::DataFrameCallback$new(standardise), chunk_size = 1e6, col_types = columns))
}

# Metabolon COMP_ID for each metabolite, with names standardised as in the GWAS preparation step.
metabolite_comp_ids <- function(input_dir) {
  supplement <- utils::read.delim(file.path(input_dir, "metabolite_sup_associations.txt"), skip = 1, check.names = FALSE)
  ids <- unique(data.frame(Metabolite = gsub(":", "_", supplement$Biochemical), COMP_ID = supplement$COMP_ID))
  stats::setNames(ids$COMP_ID, ids$Metabolite)
}

filtered_results_path <- function(root, design, outcome) {
  suffix <- if (identical(design, "liberal")) "_Liberal" else ""
  file.path(root, design, sprintf("Filtered_%s_Results%s.tsv", outcome, suffix))
}
