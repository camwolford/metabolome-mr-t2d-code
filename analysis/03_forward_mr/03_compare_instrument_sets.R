source("config/environment.R")

paths <- archive_paths("METABOLOME_MR_OUTPUT_DIR")

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

require_unique_metabolites <- function(data, label) {
  if (anyDuplicated(data$Metabolite)) stop(sprintf("%s has duplicate metabolites.", label), call. = FALSE)
  invisible(data)
}

require_output <- function(path, label) {
  output_info <- file.info(path)
  if (!file.exists(path) || is.na(output_info$size) || output_info$size == 0L) {
    stop(sprintf("%s was not written: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

compare_outcome <- function(conservative_file, liberal_file, output_file, label) {
  require_file(conservative_file, sprintf("Conservative %s MR result", label))
  require_file(liberal_file, sprintf("Liberal %s MR result", label))
  conservative_results <- readr::read_tsv(conservative_file, show_col_types = FALSE)
  liberal_results <- readr::read_tsv(liberal_file, show_col_types = FALSE)
  require_nonempty(conservative_results, sprintf("Conservative %s MR result", label))
  require_nonempty(liberal_results, sprintf("Liberal %s MR result", label))
  require_columns(conservative_results, c("Metabolite", "Number_of_IVs"), sprintf("Conservative %s MR result", label))
  require_columns(liberal_results, c("Metabolite", "Number_of_IVs"), sprintf("Liberal %s MR result", label))
  require_unique_metabolites(conservative_results, sprintf("Conservative %s MR result", label))
  require_unique_metabolites(liberal_results, sprintf("Liberal %s MR result", label))

  liberal_results <- dplyr::filter(liberal_results, Metabolite %in% conservative_results$Metabolite)
  comparison_df <- dplyr::inner_join(
    conservative_results,
    liberal_results,
    by = "Metabolite",
    suffix = c("_strict", "_liberal")
  )
  metabolites_with_fewer_IVs <- dplyr::filter(
    comparison_df,
    Number_of_IVs_strict < Number_of_IVs_liberal
  )
  utils::write.table(metabolites_with_fewer_IVs, file = output_file, sep = "\t", row.names = FALSE)
  require_output(output_file, sprintf("%s instrument-count comparison", label))
}

forward_mr_root <- file.path(paths[["output_dir"]], "03_forward_mr")
conservative_dir <- file.path(forward_mr_root, "conservative")
liberal_dir <- file.path(forward_mr_root, "liberal")
comparison_dir <- file.path(forward_mr_root, "instrument_comparison")
dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(comparison_dir)) stop(sprintf("Cannot create output directory: %s", comparison_dir), call. = FALSE)

compare_outcome(
  file.path(conservative_dir, "T2DM_MR_Results.tsv"),
  file.path(liberal_dir, "T2DM_MR_Results_Liberal.tsv"),
  file.path(comparison_dir, "metabolites_with_fewer_IVs_T2DM.tsv"),
  "T2DM"
)
compare_outcome(
  file.path(conservative_dir, "FG_MR_Results.tsv"),
  file.path(liberal_dir, "FG_MR_Results_Liberal.tsv"),
  file.path(comparison_dir, "metabolites_with_fewer_IVs_FG.tsv"),
  "fasting glucose"
)
compare_outcome(
  file.path(conservative_dir, "HBA1C_MR_Results.tsv"),
  file.path(liberal_dir, "HBA1C_MR_Results_Liberal.tsv"),
  file.path(comparison_dir, "metabolites_with_fewer_IVs_HBA1C.tsv"),
  "HbA1c"
)
