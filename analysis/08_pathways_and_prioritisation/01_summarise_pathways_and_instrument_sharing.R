source(file.path("config", "environment.R"))

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

require_directory <- function(path, label) {
  if (!dir.exists(path)) {
    stop(sprintf("%s directory is missing: %s", label, path), call. = FALSE)
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

require_unique_values <- function(values, label) {
  if (anyNA(values) || any(!nzchar(as.character(values)))) {
    stop(sprintf("%s has missing keys.", label), call. = FALSE)
  }
  if (anyDuplicated(values)) stop(sprintf("%s has duplicate keys.", label), call. = FALSE)
  invisible(values)
}

require_output <- function(path, label) {
  info <- file.info(path)
  if (!file.exists(path) || is.na(info$size) || info$size == 0L) {
    stop(sprintf("%s was not written: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

normalise_metabolite_name <- function(values) {
  values <- gsub(":", "_", values, fixed = TRUE)
  gsub("/", "_", values, fixed = TRUE)
}

instrument_files <- function(directory, label) {
  require_directory(directory, label)
  files <- list.files(directory, full.names = TRUE)
  files <- files[!grepl("~", files, fixed = TRUE)]
  files <- files[!file.info(files)$isdir]
  if (!length(files)) stop(sprintf("%s has no instrument files: %s", label, directory), call. = FALSE)
  files
}

instrument_snp_lists <- function(directory, outcome, metabolites, label) {
  snps_by_metabolite <- stats::setNames(rep(NA_character_, length(metabolites)), metabolites)
  for (instrument_file in instrument_files(directory, label)) {
    require_file(instrument_file, basename(instrument_file))
    instrument_data <- as.data.frame(readr::read_tsv(instrument_file, show_col_types = FALSE))
    require_nonempty(instrument_data, basename(instrument_file))
    require_columns(instrument_data, "SNP", basename(instrument_file))
    if (anyNA(instrument_data$SNP) || any(!nzchar(as.character(instrument_data$SNP)))) {
      stop(sprintf("%s has missing SNP values.", basename(instrument_file)), call. = FALSE)
    }

    metabolite_name <- strsplit(basename(instrument_file), "_Harmonised_IVs", fixed = TRUE)[[1]][1]
    metabolite_name <- gsub(outcome, "", metabolite_name, fixed = TRUE)
    if (!metabolite_name %in% metabolites) next
    if (!is.na(snps_by_metabolite[[metabolite_name]])) {
      stop(sprintf("%s has duplicate instrument files for %s.", label, metabolite_name), call. = FALSE)
    }
    snps_by_metabolite[[metabolite_name]] <- paste(instrument_data$SNP, collapse = ", ")
  }
  snps_by_metabolite
}

summarise_snp_sharing <- function(data, snp_column) {
  require_columns(data, c(snp_column, "superpathway", "subpathway"), "Significant metabolite annotation")
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
  require_output(path, label)
}

run_stage_08 <- function() {
  for (package in c("readr", "dplyr", "tidyr", "stringr")) require_namespace(package)

  paths <- archive_paths(c("METABOLOME_MR_INPUT_DIR", "METABOLOME_MR_OUTPUT_DIR"))
  annotation_file <- file.path(paths[["input_dir"]], "annotations", "responses_combined.tsv")
  significant_file <- file.path(
    paths[["output_dir"]], "04_significance_filtering", "combined", "Full_Significant_Results_Manuscript.tsv"
  )
  instrument_root <- file.path(paths[["output_dir"]], "02_instrument_selection", "conservative")
  output_dir <- file.path(paths[["output_dir"]], "08_pathways_and_prioritisation")

  require_file(annotation_file, "Metabolite annotation")
  require_file(significant_file, "Combined significant-result summary")
  responses <- as.data.frame(readr::read_tsv(annotation_file, show_col_types = FALSE))
  sig_results <- as.data.frame(readr::read_tsv(significant_file, show_col_types = FALSE))
  require_nonempty(responses, "Metabolite annotation")
  require_nonempty(sig_results, "Combined significant-result summary")
  require_columns(responses, c("name", "superpathway", "subpathway"), "Metabolite annotation")
  require_columns(sig_results, "Metabolite", "Combined significant-result summary")

  responses$name <- normalise_metabolite_name(as.character(responses$name))
  require_unique_values(responses$name, "Metabolite annotation")
  sig_metabolites <- as.character(sig_results$Metabolite)
  require_unique_values(sig_metabolites, "Combined significant-result summary")
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
  if (!dir.exists(output_dir)) stop(sprintf("Cannot create output directory: %s", output_dir), call. = FALSE)
  write_summary(superpathway_counts, file.path(output_dir, "superpathway_counts.tsv"), "Superpathway summary")
  write_summary(subpathway_counts, file.path(output_dir, "subpathway_counts.tsv"), "Subpathway summary")
  write_summary(t2dm_snp_counts, file.path(output_dir, "T2DM_snp_counts.tsv"), "Type 2 diabetes SNP-sharing summary")
  write_summary(fg_snp_counts, file.path(output_dir, "FG_snp_counts.tsv"), "Fasting-glucose SNP-sharing summary")
  write_summary(hba1c_snp_counts, file.path(output_dir, "HBA1C_snp_counts.tsv"), "HbA1c SNP-sharing summary")
  write_summary(sig_metabolites_data, file.path(output_dir, "significant_metabolites_data.tsv"), "Significant-metabolite annotation summary")
}

if (sys.nframe() == 0L) run_stage_08()
