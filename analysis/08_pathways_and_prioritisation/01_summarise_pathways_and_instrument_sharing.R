source("R/utils.R")

normalise_metabolite_name <- function(values) {
  values <- gsub(":", "_", values, fixed = TRUE)
  gsub("/", "_", values, fixed = TRUE)
}

instrument_files <- function(directory, label) {
  files <- list.files(directory, full.names = TRUE)
  files <- files[!grepl("~", files, fixed = TRUE)]
  files <- files[!file.info(files)$isdir]
  files
}

instrument_snp_lists <- function(directory, outcome, metabolites, label) {
  snps_by_metabolite <- stats::setNames(rep(NA_character_, length(metabolites)), metabolites)
  for (instrument_file in instrument_files(directory, label)) {
    instrument_data <- as.data.frame(readr::read_tsv(instrument_file, show_col_types = FALSE))

    metabolite_name <- strsplit(basename(instrument_file), "_Harmonised_IVs", fixed = TRUE)[[1]][1]
    metabolite_name <- gsub(outcome, "", metabolite_name, fixed = TRUE)
    if (!metabolite_name %in% metabolites) next
    snps_by_metabolite[[metabolite_name]] <- paste(instrument_data$SNP, collapse = ", ")
  }
  snps_by_metabolite
}

summarise_snp_sharing <- function(data, snp_column) {
  expanded_snps <- data |>
    dplyr::filter(!is.na(.data[[snp_column]])) |>
    dplyr::mutate(!!snp_column := strsplit(as.character(.data[[snp_column]]), ",\\s*")) |>
    tidyr::unnest(cols = dplyr::all_of(snp_column))

  if (!nrow(expanded_snps)) {
    return(data.frame(setNames(list(character(), integer(), character(), character()), c(snp_column, "n", "superpathway", "subpathway"))))
  }

  snp_counts <- expanded_snps |>
    dplyr::group_by(.data[[snp_column]]) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$n))
  snp_counts$superpathway <- NA_character_
  snp_counts$subpathway <- NA_character_

  for (snp in snp_counts[[snp_column]]) {
    matching_rows <- !is.na(data[[snp_column]]) & stringr::str_detect(data[[snp_column]], snp)
    snp_counts$superpathway[snp_counts[[snp_column]] == snp] <- paste(unique(data$superpathway[matching_rows]), collapse = ", ")
    snp_counts$subpathway[snp_counts[[snp_column]] == snp] <- paste(unique(data$subpathway[matching_rows]), collapse = ", ")
  }
  snp_counts
}

write_summary <- function(data, path, label) {
  readr::write_tsv(data, path)
}

paths <- data_paths(c("METABOLOME_MR_INPUT_DIR", "METABOLOME_MR_OUTPUT_DIR"))
annotation_file <- file.path(paths[["input_dir"]], "annotations", "responses_combined.tsv")
significant_file <- file.path(
  paths[["output_dir"]], "04_significance_filtering", "combined", "Full_Significant_Results_Manuscript.tsv"
)
instrument_root <- file.path(paths[["output_dir"]], "02_instrument_selection", "conservative")
output_dir <- file.path(paths[["output_dir"]], "08_pathways_and_prioritisation")

responses <- as.data.frame(readr::read_tsv(annotation_file, show_col_types = FALSE))
sig_results <- as.data.frame(readr::read_tsv(significant_file, show_col_types = FALSE))

responses$name <- normalise_metabolite_name(as.character(responses$name))
sig_metabolites <- as.character(sig_results$Metabolite)
sig_metabolites_data <- dplyr::filter(responses, .data$name %in% sig_metabolites)
if (nrow(sig_metabolites_data) != length(sig_metabolites)) {
  missing <- setdiff(sig_metabolites, sig_metabolites_data$name)
  stop(sprintf("Metabolite annotation is missing significant metabolites: %s", paste(missing, collapse = ", ")), call. = FALSE)
}

superpathway_counts <- sig_metabolites_data |>
  dplyr::group_by(.data$superpathway) |>
  dplyr::summarise(n = dplyr::n_distinct(.data$name), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(.data$n))
subpathway_counts <- sig_metabolites_data |>
  dplyr::group_by(.data$subpathway) |>
  dplyr::summarise(n = dplyr::n_distinct(.data$name), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(.data$n))

instrument_lists <- data.frame(name = sig_metabolites, stringsAsFactors = FALSE)
instrument_lists$T2DM_SNPs <- instrument_snp_lists(
  file.path(instrument_root, "Harmonised_T2DM_IVs"), "T2DM", sig_metabolites, "Type 2 diabetes harmonised instruments"
)
instrument_lists$FG_SNPs <- instrument_snp_lists(
  file.path(instrument_root, "Harmonised_FG_IVs"), "FG", sig_metabolites, "Fasting-glucose harmonised instruments"
)
instrument_lists$HBA1C_SNPs <- instrument_snp_lists(
  file.path(instrument_root, "Harmonised_HBA1C_IVs"), "HBA1C", sig_metabolites, "HbA1c harmonised instruments"
)
sig_metabolites_data$T2DM_SNPs <- instrument_lists$T2DM_SNPs[match(sig_metabolites_data$name, instrument_lists$name)]
sig_metabolites_data$FG_SNPs <- instrument_lists$FG_SNPs[match(sig_metabolites_data$name, instrument_lists$name)]
sig_metabolites_data$HBA1C_SNPs <- instrument_lists$HBA1C_SNPs[match(sig_metabolites_data$name, instrument_lists$name)]

t2dm_snp_counts <- summarise_snp_sharing(sig_metabolites_data, "T2DM_SNPs")
fg_snp_counts <- summarise_snp_sharing(sig_metabolites_data, "FG_SNPs")
hba1c_snp_counts <- summarise_snp_sharing(sig_metabolites_data, "HBA1C_SNPs")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_summary(superpathway_counts, file.path(output_dir, "superpathway_counts.tsv"), "Superpathway summary")
write_summary(subpathway_counts, file.path(output_dir, "subpathway_counts.tsv"), "Subpathway summary")
write_summary(t2dm_snp_counts, file.path(output_dir, "T2DM_snp_counts.tsv"), "Type 2 diabetes SNP-sharing summary")
write_summary(fg_snp_counts, file.path(output_dir, "FG_snp_counts.tsv"), "Fasting-glucose SNP-sharing summary")
write_summary(hba1c_snp_counts, file.path(output_dir, "HBA1C_snp_counts.tsv"), "HbA1c SNP-sharing summary")
write_summary(sig_metabolites_data, file.path(output_dir, "significant_metabolites_data.tsv"), "Significant-metabolite annotation summary")
