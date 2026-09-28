source("R/utils.R")

paths <- data_paths(c(
  "METABOLOME_MR_INPUT_DIR", "METABOLOME_MR_WORK_DIR", "METABOLOME_MR_OUTPUT_DIR", "METABOLOME_MR_EUR_PANEL_DIR"
))
plink_bin <- Sys.getenv("METABOLOME_MR_PLINK")

input_dir <- file.path(paths[["work_dir"]], "02_instrument_selection", "conservative", "Filtered_IVs")
output_dir <- file.path(paths[["work_dir"]], "02_instrument_selection", "conservative", "Clumped_IVs")
eur_bfile <- file.path(paths[["eur_panel_dir"]], "EUR")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
metabolite_files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)

for (metabolite_file in metabolite_files) {
  metabolite_data <- as.data.frame(readr::read_tsv(metabolite_file, show_col_types = FALSE))
  metabolite_data <- metabolite_data[!(metabolite_data$EffectAllele %in% c("D", "I")), , drop = FALSE]
  if (!nrow(metabolite_data)) next
  if (nrow(metabolite_data) > 1) {
    clump_input <- metabolite_data
    names(clump_input)[names(clump_input) == "SNP"] <- "rsid"
    names(clump_input)[names(clump_input) == "Pval"] <- "pval"
    metabolite_data <- ieugwasr::ld_clump_local(
      clump_input, clump_kb = 10000, clump_r2 = 0.001, clump_p = 1,
      bfile = eur_bfile, plink_bin = plink_bin
    )
    names(metabolite_data)[names(metabolite_data) == "rsid"] <- "SNP"
    names(metabolite_data)[names(metabolite_data) == "pval"] <- "Pval"
  }
  if (!nrow(metabolite_data)) next
  metabolite_name <- sub("_IVs.*$", "", basename(metabolite_file))
  output_file <- file.path(output_dir, paste0(metabolite_name, "_Clumped_IVs.tsv"))
  utils::write.table(metabolite_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
}
