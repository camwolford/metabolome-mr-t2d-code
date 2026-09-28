source("R/utils.R")

paths <- data_paths(c(
  "METABOLOME_MR_INPUT_DIR", "METABOLOME_MR_WORK_DIR", "METABOLOME_MR_OUTPUT_DIR"
))

metabolite_dir <- file.path(paths[["work_dir"]], "Individual_Metabolite_GWAS")
count_file <- file.path(paths[["work_dir"]], "Repeated_SNPS.tsv")
output_dir <- file.path(paths[["work_dir"]], "02_instrument_selection", "liberal", "Filtered_IVs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

repeated_snps <- readr::read_tsv(count_file, show_col_types = FALSE)
shared_snps <- repeated_snps$SNP[repeated_snps$Count > 9]
metabolite_files <- list.files(metabolite_dir, pattern = "\\.tsv$", full.names = TRUE)

for (metabolite_file in metabolite_files) {
  metabolite_data <- readr::read_tsv(metabolite_file, show_col_types = FALSE)
  filtered_snps <- metabolite_data[
    metabolite_data$Pval < 5e-6 & metabolite_data$EAF > 0.01 & metabolite_data$EAF < 0.99 & !(metabolite_data$SNP %in% shared_snps),
    , drop = FALSE
  ]
  if (nrow(filtered_snps)) {
    metabolite_name <- sub("_GWAS.*$", "", basename(metabolite_file))
    output_file <- file.path(output_dir, paste0(metabolite_name, "_IVs_Liberal.tsv"))
    readr::write_tsv(filtered_snps, output_file)
  }
}
