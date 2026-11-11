source("R/utils.R")

paths <- data_paths(c(
  "METABOLOME_MR_INPUT_DIR",
  "METABOLOME_MR_OUTPUT_DIR",
  "METABOLOME_MR_EUR_PANEL_DIR"
))
plink_bin <- Sys.getenv("METABOLOME_MR_PLINK")

as_numeric <- function(values, label) {
  numeric_values <- suppressWarnings(as.numeric(gsub(",", "", as.character(values), fixed = TRUE)))
  numeric_values
}

fill_chromosomes <- function(values) {
  filled <- zoo::na.locf(values, na.rm = FALSE)
  filled
}

standardise_t2dm <- function(path) {
  raw <- as.data.frame(readr::read_csv(path, skip = 2, name_repair = "minimal", show_col_types = FALSE))
  names(raw)[6:8] <- c("Residual", "Risk", "Other")
  european <- raw[, -c(9:18), drop = FALSE]
  european <- european[, -c(14:28), drop = FALSE]
  european <- european[-1, , drop = FALSE]
  names(european)[9:13] <- c("Effective_sample_size", "EAF", "Log_OR", "SE", "Pval")

  data <- data.frame(
    Chromosome = as_numeric(fill_chromosomes(european[[2]]), "Type 2 diabetes chromosome"),
    Position = as_numeric(european[[4]], "Type 2 diabetes position"),
    SNP = as.character(european[[3]]),
    EffectAllele = toupper(as.character(european[[7]])),
    NonEffectAllele = toupper(as.character(european[[8]])),
    EAF = as_numeric(european$EAF, "Type 2 diabetes EAF"),
    Beta = as_numeric(european$Log_OR, "Type 2 diabetes log odds ratio"),
    SE = as_numeric(european$SE, "Type 2 diabetes standard error"),
    Pval = as_numeric(european$Pval, "Type 2 diabetes p-value"),
    stringsAsFactors = FALSE
  )
  data <- data[!is.na(data$Pval), , drop = FALSE]
  data <- data[data$Pval < 5e-8 & data$EAF > 0.01 & data$EAF < 0.99, , drop = FALSE]
  data
}

standardise_glycaemic <- function(path, trait, label) {
  raw <- as.data.frame(readr::read_csv(path, name_repair = "minimal", show_col_types = FALSE))
  trait_rows <- raw[as.character(raw[[6]]) == trait, , drop = FALSE]
  data <- data.frame(
    Chromosome = as_numeric(trait_rows[[7]], sprintf("%s chromosome", label)),
    Position = as_numeric(trait_rows[[8]], sprintf("%s position", label)),
    SNP = as.character(trait_rows[[9]]),
    EffectAllele = toupper(as.character(trait_rows[[10]])),
    NonEffectAllele = toupper(as.character(trait_rows[[11]])),
    EAF = as_numeric(trait_rows[[15]], sprintf("%s EAF", label)),
    Beta = as_numeric(trait_rows[[16]], sprintf("%s effect", label)),
    SE = as_numeric(trait_rows[[17]], sprintf("%s standard error", label)),
    Pval = as_numeric(trait_rows[[18]], sprintf("%s p-value", label)),
    stringsAsFactors = FALSE
  )
  data <- data[!is.na(data$Pval), , drop = FALSE]
  data <- data[data$Pval < 5e-8 & data$EAF > 0.01 & data$EAF < 0.99, , drop = FALSE]
  data
}

clump_instruments <- function(data, label, eur_bfile) {
  if (nrow(data) == 1L) return(data)
  clump_input <- data[, c("SNP", "Pval"), drop = FALSE]
  names(clump_input) <- c("rsid", "pval")
  clumped_ids <- ieugwasr::ld_clump_local(
    clump_input,
    clump_kb = 10000,
    clump_r2 = 0.001,
    clump_p = 1,
    bfile = eur_bfile,
    plink_bin = plink_bin
  )
  clumped <- data[match(clumped_ids$rsid, data$SNP), , drop = FALSE]
  clumped
}

eur_bfile <- file.path(paths[["eur_panel_dir"]], "EUR")

input_dir <- paths[["input_dir"]]
output_dir <- file.path(paths[["output_dir"]], "06_reverse_mr", "outcome_instruments")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

selected <- list(
  T2DM = standardise_t2dm(file.path(input_dir, "t2dm_sig_snps_ancestry.csv")),
  FG = standardise_glycaemic(file.path(input_dir, "fg_hba1c_sig_snps_ancestry.csv"), "FG", "fasting glucose"),
  HBA1C = standardise_glycaemic(file.path(input_dir, "fg_hba1c_sig_snps_ancestry.csv"), "HbA1c", "HbA1c")
)

for (code in names(selected)) {
  selected_path <- file.path(output_dir, sprintf("%s_selected.tsv", code))
  readr::write_tsv(selected[[code]], selected_path)

  clumped <- clump_instruments(selected[[code]], code, eur_bfile)
  clumped_path <- file.path(output_dir, sprintf("%s_clumped.tsv", code))
  readr::write_tsv(clumped, clumped_path)
}
