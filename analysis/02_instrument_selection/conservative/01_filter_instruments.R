source("config/environment.R")

paths <- archive_paths(c(
  "METABOLOME_MR_INPUT_DIR", "METABOLOME_MR_WORK_DIR", "METABOLOME_MR_OUTPUT_DIR"
))

require_columns <- function(data, columns, label) {
  absent <- setdiff(columns, names(data))
  if (length(absent)) stop(sprintf("%s is missing: %s", label, paste(absent, collapse = ", ")), call. = FALSE)
  invisible(data)
}

require_nonempty <- function(data, label) {
  if (!nrow(data)) stop(sprintf("%s has zero rows.", label), call. = FALSE)
  invisible(data)
}

metabolite_dir <- file.path(paths[["work_dir"]], "Individual_Metabolite_GWAS")
count_file <- file.path(paths[["work_dir"]], "Repeated_SNPS.tsv")
output_dir <- file.path(paths[["work_dir"]], "02_instrument_selection", "conservative", "Filtered_IVs")
if (!dir.exists(metabolite_dir)) stop(sprintf("Metabolite directory is missing: %s", metabolite_dir), call. = FALSE)
if (!file.exists(count_file)) stop(sprintf("Shared-SNP count file is missing: %s", count_file), call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop(sprintf("Cannot create output directory: %s", output_dir), call. = FALSE)

repeated_snps <- readr::read_tsv(count_file, show_col_types = FALSE)
require_nonempty(repeated_snps, "Shared-SNP count table")
require_columns(repeated_snps, c("SNP", "Count"), "Shared-SNP count table")
shared_snps <- repeated_snps$SNP[repeated_snps$Count > 5]
metabolite_files <- list.files(metabolite_dir, pattern = "\\.tsv$", full.names = TRUE)
if (!length(metabolite_files)) stop(sprintf("No metabolite files found in: %s", metabolite_dir), call. = FALSE)

for (metabolite_file in metabolite_files) {
  metabolite_data <- readr::read_tsv(metabolite_file, show_col_types = FALSE)
  require_nonempty(metabolite_data, basename(metabolite_file))
  require_columns(metabolite_data, c("SNP", "Pval", "EAF"), basename(metabolite_file))
  filtered_snps <- metabolite_data[
    metabolite_data$Pval < 5e-8 & metabolite_data$EAF > 0.01 & metabolite_data$EAF < 0.99 & !(metabolite_data$SNP %in% shared_snps),
    , drop = FALSE
  ]
  if (nrow(filtered_snps)) {
    metabolite_name <- sub("_GWAS.*$", "", basename(metabolite_file))
    output_file <- file.path(output_dir, paste0(metabolite_name, "_IVs.tsv"))
    readr::write_tsv(filtered_snps, output_file)
    if (!file.exists(output_file)) stop(sprintf("Filtered output was not written: %s", output_file), call. = FALSE)
  }
}
