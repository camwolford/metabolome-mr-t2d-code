source(file.path("config", "environment.R"))

paths <- archive_paths(c(
  "METABOLOME_MR_INPUT_DIR",
  "METABOLOME_MR_WORK_DIR",
  "METABOLOME_MR_OUTPUT_DIR"
))

require_file <- function(path, label) {
  if (!file.exists(path)) stop(sprintf("%s is missing: %s", label, path), call. = FALSE)
  invisible(path)
}

require_columns <- function(data, columns, label) {
  absent <- setdiff(columns, names(data))
  if (length(absent)) stop(sprintf("%s is missing: %s", label, paste(absent, collapse = ", ")), call. = FALSE)
  invisible(data)
}

require_nonempty <- function(data, label) {
  if (!nrow(data)) stop(sprintf("%s has zero rows.", label), call. = FALSE)
  invisible(data)
}

require_output <- function(path, label) {
  output_info <- file.info(path)
  if (!file.exists(path) || is.na(output_info$size) || output_info$size == 0L) {
    stop(sprintf("%s was not written: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

input_dir <- paths[["input_dir"]]
work_dir <- paths[["work_dir"]]

metabolite_input <- file.path(input_dir, "metabolite_sup_associations.txt")
fasting_glucose_input <- file.path(input_dir, "fasting_glucose_gwas_cleaned.csv")
hba1c_input <- file.path(input_dir, "hbA1c_gwas_cleaned.csv")

require_file(metabolite_input, "Metabolite association input")
require_file(fasting_glucose_input, "Fasting-glucose input")
require_file(hba1c_input, "HbA1c input")

metabolite_gwas <- utils::read.delim(metabolite_input, stringsAsFactors = FALSE)
if (nrow(metabolite_gwas) < 2L) {
  stop("Metabolite association input must contain a column-name row and at least one data row.", call. = FALSE)
}
if (ncol(metabolite_gwas) != 41L) {
  stop("Metabolite association input must contain exactly 41 columns.", call. = FALSE)
}

colnames(metabolite_gwas) <- metabolite_gwas[1, ]
metabolite_gwas <- metabolite_gwas[-1, ]
metabolite_gwas <- metabolite_gwas[, -c(
  1, 2, 3, 11, 12, 13, 14, 15, 17, 18, 19, 20, 22, 23, 24,
  25, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41
)]
if (ncol(metabolite_gwas) != 11L) {
  stop("Selected metabolite association fields must contain exactly 11 columns.", call. = FALSE)
}

colnames(metabolite_gwas)[4] <- "SNP"
colnames(metabolite_gwas)[5] <- "EffectAllele"
colnames(metabolite_gwas)[6] <- "NonEffectAllele"
colnames(metabolite_gwas)[7] <- "Metabolite"
colnames(metabolite_gwas)[8] <- "EAF"
colnames(metabolite_gwas)[9] <- "Beta"
colnames(metabolite_gwas)[10] <- "SE"
colnames(metabolite_gwas)[11] <- "Pval"
require_columns(
  metabolite_gwas,
  c("Chromosome", "Position", "SNP", "EffectAllele", "NonEffectAllele", "Metabolite", "EAF", "Beta", "SE", "Pval"),
  "Selected metabolite association fields"
)

metabolite_gwas$Chromosome <- as.numeric(metabolite_gwas$Chromosome)
metabolite_gwas$Position <- as.numeric(metabolite_gwas$Position)
metabolite_gwas$EAF <- as.numeric(metabolite_gwas$EAF)
metabolite_gwas$Beta <- as.numeric(metabolite_gwas$Beta)
metabolite_gwas$SE <- as.numeric(metabolite_gwas$SE)
metabolite_gwas$Pval <- as.numeric(metabolite_gwas$Pval)
metabolite_gwas$Pval_new <- 10^(-1 * metabolite_gwas$Pval)
metabolite_gwas$Pval <- metabolite_gwas$Pval_new
metabolite_gwas <- metabolite_gwas[, -12]
metabolite_gwas <- metabolite_gwas[stats::complete.cases(metabolite_gwas), ]
require_nonempty(metabolite_gwas, "Cleaned metabolite association table")
metabolite_gwas$Metabolite <- gsub(":", "_", metabolite_gwas$Metabolite)

metabolite_output <- file.path(work_dir, "metabolite_gwas_associations_cleaned.tsv")
utils::write.table(metabolite_gwas, metabolite_output, sep = "\t", row.names = FALSE)
require_output(metabolite_output, "Cleaned metabolite association table")

fasting_glucose_gwas <- utils::read.csv(fasting_glucose_input, stringsAsFactors = FALSE)
require_nonempty(fasting_glucose_gwas, "Fasting-glucose input")
if (ncol(fasting_glucose_gwas) != 9L) {
  stop("Fasting-glucose input must contain exactly 9 columns.", call. = FALSE)
}
colnames(fasting_glucose_gwas) <- c(
  "Chromosome", "Position", "EffectAllele", "NonEffectAllele", "EAF", "Beta", "SE", "Pval", "SampleSize"
)

fasting_glucose_output <- file.path(work_dir, "fasting_glucose_gwas_cleaned.tsv")
utils::write.table(fasting_glucose_gwas, fasting_glucose_output, sep = "\t", row.names = FALSE)
require_output(fasting_glucose_output, "Standardised fasting-glucose table")

hba1c_gwas <- utils::read.csv(hba1c_input, stringsAsFactors = FALSE)
require_nonempty(hba1c_gwas, "HbA1c input")
if (ncol(hba1c_gwas) != 9L) {
  stop("HbA1c input must contain exactly 9 columns.", call. = FALSE)
}
colnames(hba1c_gwas) <- c(
  "Chromosome", "Position", "EffectAllele", "NonEffectAllele", "EAF", "Beta", "SE", "Pval", "SampleSize"
)

hba1c_output <- file.path(work_dir, "hbA1c_gwas_cleaned.tsv")
utils::write.table(hba1c_gwas, hba1c_output, sep = "\t", row.names = FALSE)
require_output(hba1c_output, "Standardised HbA1c table")
