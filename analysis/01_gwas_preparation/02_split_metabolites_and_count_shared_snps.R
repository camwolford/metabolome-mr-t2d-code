source("R/utils.R")

paths <- data_paths(c(
  "METABOLOME_MR_INPUT_DIR",
  "METABOLOME_MR_WORK_DIR",
  "METABOLOME_MR_OUTPUT_DIR"
))

work_dir <- paths[["work_dir"]]
metabolite_input <- file.path(work_dir, "metabolite_gwas_associations_cleaned.tsv")

metabolite_gwas <- utils::read.delim(metabolite_input, stringsAsFactors = FALSE)

individual_metabolite_dir <- file.path(work_dir, "Individual_Metabolite_GWAS")
if (!dir.exists(individual_metabolite_dir)) dir.create(individual_metabolite_dir, recursive = TRUE)

unique_metabolites <- unique(metabolite_gwas$Metabolite)
for (metabolite in unique_metabolites) {
  metabolite_gwas_metabolite <- metabolite_gwas[metabolite_gwas$Metabolite == metabolite, ]

  metabolite_output <- file.path(individual_metabolite_dir, paste0(metabolite, "_GWAS.tsv"))
  utils::write.table(
    metabolite_gwas_metabolite,
    file = metabolite_output,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

snp_counts <- sort(table(metabolite_gwas$SNP), decreasing = TRUE)
repeated_snps <- data.frame(
  SNP = names(snp_counts),
  Count = as.numeric(snp_counts),
  stringsAsFactors = FALSE
)

repeated_snps_output <- file.path(work_dir, "Repeated_SNPS.tsv")
utils::write.table(
  repeated_snps,
  file = repeated_snps_output,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
