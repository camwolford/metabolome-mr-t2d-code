source("config/environment.R")

paths <- archive_paths(c(
  "METABOLOME_MR_INPUT_DIR",
  "METABOLOME_MR_OUTPUT_DIR",
  "METABOLOME_MR_EUR_PANEL_DIR"
))
plink_bin <- require_executable("METABOLOME_MR_PLINK")

require_namespace <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("Required R package is unavailable: %s", package), call. = FALSE)
  }
}

require_file <- function(path, label) {
  if (!file.exists(path) || dir.exists(path)) {
    stop(sprintf("%s is missing: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

require_columns <- function(data, columns, label) {
  absent <- setdiff(columns, names(data))
  if (length(absent)) {
    stop(sprintf("%s is missing: %s", label, paste(absent, collapse = ", ")), call. = FALSE)
  }
  invisible(data)
}

require_nonempty <- function(data, label) {
  if (!nrow(data)) stop(sprintf("%s has zero rows.", label), call. = FALSE)
  invisible(data)
}

require_unique_keys <- function(data, keys, label) {
  duplicated_rows <- duplicated(data[keys])
  if (any(duplicated_rows)) {
    stop(sprintf("%s has duplicate keys: %s", label, paste(keys, collapse = ", ")), call. = FALSE)
  }
  invisible(data)
}

require_output <- function(path, label) {
  info <- file.info(path)
  if (!file.exists(path) || is.na(info$size) || info$size == 0L) {
    stop(sprintf("%s was not written: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

as_numeric <- function(values, label) {
  numeric_values <- suppressWarnings(as.numeric(gsub(",", "", as.character(values), fixed = TRUE)))
  if (all(is.na(numeric_values))) stop(sprintf("%s has no numeric values.", label), call. = FALSE)
  numeric_values
}

fill_chromosomes <- function(values) {
  filled <- zoo::na.locf(values, na.rm = FALSE)
  if (anyNA(filled)) stop("Type 2 diabetes source has an unfilled chromosome value.", call. = FALSE)
  filled
}

standardise_t2dm <- function(path) {
  require_file(path, "Type 2 diabetes outcome-instrument source")
  raw <- as.data.frame(readr::read_csv(path, skip = 2, name_repair = "minimal", show_col_types = FALSE))
  if (ncol(raw) != 38L) {
    stop(sprintf("Type 2 diabetes outcome-instrument source must contain 38 columns, found %d.", ncol(raw)), call. = FALSE)
  }
  names(raw)[6:8] <- c("Residual", "Risk", "Other")
  european <- raw[, -c(9:18), drop = FALSE]
  european <- european[, -c(14:28), drop = FALSE]
  european <- european[-1, , drop = FALSE]
  if (ncol(european) != 13L) {
    stop(sprintf("European type 2 diabetes source must contain 13 columns, found %d.", ncol(european)), call. = FALSE)
  }
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
  data <- data[data$Pval < 5e-8 & data$EAF > 0.01, , drop = FALSE]
  require_nonempty(data, "Selected type 2 diabetes outcome instruments")
  require_columns(data, c("Chromosome", "Position", "SNP", "EffectAllele", "NonEffectAllele", "EAF", "Beta", "SE", "Pval"), "Selected type 2 diabetes outcome instruments")
  if (any(!stats::complete.cases(data))) stop("Selected type 2 diabetes outcome instruments contain missing values.", call. = FALSE)
  require_unique_keys(data, "SNP", "Selected type 2 diabetes outcome instruments")
  data
}

standardise_glycaemic <- function(path, trait, label) {
  require_file(path, sprintf("%s outcome-instrument source", label))
  raw <- as.data.frame(readr::read_csv(path, name_repair = "minimal", show_col_types = FALSE))
  if (ncol(raw) < 18L) {
    stop(sprintf("%s outcome-instrument source must contain at least 18 columns, found %d.", label, ncol(raw)), call. = FALSE)
  }
  trait_rows <- raw[as.character(raw[[6]]) == trait, , drop = FALSE]
  require_nonempty(trait_rows, sprintf("%s rows", label))
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
  data <- data[data$Pval < 5e-8 & data$EAF > 0.01, , drop = FALSE]
  require_nonempty(data, sprintf("Selected %s outcome instruments", label))
  if (any(!stats::complete.cases(data))) stop(sprintf("Selected %s outcome instruments contain missing values.", label), call. = FALSE)
  require_unique_keys(data, "SNP", sprintf("Selected %s outcome instruments", label))
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
  require_nonempty(clumped_ids, sprintf("Clumped %s outcome instruments", label))
  require_columns(clumped_ids, "rsid", sprintf("Clumped %s outcome instruments", label))
  clumped <- data[match(clumped_ids$rsid, data$SNP), , drop = FALSE]
  if (anyNA(clumped$SNP)) stop(sprintf("Clumping returned an unknown %s variant.", label), call. = FALSE)
  require_unique_keys(clumped, "SNP", sprintf("Clumped %s outcome instruments", label))
  clumped
}

for (package in c("readr", "zoo", "ieugwasr")) require_namespace(package)

eur_bfile <- file.path(paths[["eur_panel_dir"]], "EUR")
if (!all(file.exists(paste0(eur_bfile, c(".bed", ".bim", ".fam"))))) {
  stop(sprintf("EUR PLINK panel prefix is incomplete: %s", eur_bfile), call. = FALSE)
}

input_dir <- paths[["input_dir"]]
output_dir <- file.path(paths[["output_dir"]], "06_reverse_mr", "outcome_instruments")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop(sprintf("Cannot create output directory: %s", output_dir), call. = FALSE)

selected <- list(
  T2DM = standardise_t2dm(file.path(input_dir, "t2dm_sig_snps_ancestry.csv")),
  FG = standardise_glycaemic(file.path(input_dir, "fg_hba1c_sig_snps_ancestry.csv"), "FG", "fasting glucose"),
  HBA1C = standardise_glycaemic(file.path(input_dir, "fg_hba1c_sig_snps_ancestry.csv"), "HbA1c", "HbA1c")
)

for (code in names(selected)) {
  selected_path <- file.path(output_dir, sprintf("%s_selected.tsv", code))
  readr::write_tsv(selected[[code]], selected_path)
  require_output(selected_path, sprintf("Selected %s outcome instruments", code))

  clumped <- clump_instruments(selected[[code]], code, eur_bfile)
  clumped_path <- file.path(output_dir, sprintf("%s_clumped.tsv", code))
  readr::write_tsv(clumped, clumped_path)
  require_output(clumped_path, sprintf("Clumped %s outcome instruments", code))
}
