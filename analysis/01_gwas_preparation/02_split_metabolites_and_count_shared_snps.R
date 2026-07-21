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

work_dir <- paths[["work_dir"]]
metabolite_input <- file.path(work_dir, "metabolite_gwas_associations_cleaned.tsv")
require_file(metabolite_input, "Cleaned metabolite association table")

metabolite_gwas <- utils::read.delim(metabolite_input, stringsAsFactors = FALSE)
require_nonempty(metabolite_gwas, "Cleaned metabolite association table")
require_columns(metabolite_gwas, c("Metabolite", "SNP"), "Cleaned metabolite association table")

individual_metabolite_dir <- file.path(work_dir, "Individual_Metabolite_GWAS")
if (!dir.exists(individual_metabolite_dir)) dir.create(individual_metabolite_dir, recursive = TRUE)
if (!dir.exists(individual_metabolite_dir)) {
  stop(sprintf("Per-metabolite output directory was not created: %s", individual_metabolite_dir), call. = FALSE)
}

unique_metabolites <- unique(metabolite_gwas$Metabolite)
for (metabolite in unique_metabolites) {
  metabolite_gwas_metabolite <- metabolite_gwas[metabolite_gwas$Metabolite == metabolite, ]
  require_nonempty(metabolite_gwas_metabolite, sprintf("Metabolite table for %s", metabolite))

  metabolite_output <- file.path(individual_metabolite_dir, paste0(metabolite, "_GWAS.tsv"))
  utils::write.table(
    metabolite_gwas_metabolite,
    file = metabolite_output,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  require_output(metabolite_output, sprintf("Metabolite table for %s", metabolite))
}

snp_counts <- sort(table(metabolite_gwas$SNP), decreasing = TRUE)
repeated_snps <- data.frame(
  SNP = names(snp_counts),
  Count = as.numeric(snp_counts),
  stringsAsFactors = FALSE
)
require_nonempty(repeated_snps, "Repeated SNP count table")

repeated_snps_output <- file.path(work_dir, "Repeated_SNPS.tsv")
utils::write.table(
  repeated_snps,
  file = repeated_snps_output,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
require_output(repeated_snps_output, "Repeated SNP count table")
