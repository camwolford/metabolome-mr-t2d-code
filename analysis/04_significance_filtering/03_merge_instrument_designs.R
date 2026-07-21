source("config/environment.R")

paths <- archive_paths("METABOLOME_MR_OUTPUT_DIR")

require_file <- function(path, label) {
  if (!file.exists(path)) stop(sprintf("%s is missing: %s", label, path), call. = FALSE)
  invisible(path)
}

require_nonempty <- function(data, label) {
  if (!nrow(data)) stop(sprintf("%s has zero rows.", label), call. = FALSE)
  invisible(data)
}

require_metabolites <- function(data, label) {
  if (!"Metabolite" %in% names(data)) stop(sprintf("%s is missing: Metabolite", label), call. = FALSE)
  if (anyNA(data$Metabolite) || any(!nzchar(data$Metabolite))) {
    stop(sprintf("%s has missing metabolite identifiers.", label), call. = FALSE)
  }
  if (anyDuplicated(data$Metabolite)) stop(sprintf("%s has duplicate metabolites.", label), call. = FALSE)
  invisible(data)
}

require_matching_schema <- function(left, right, left_label, right_label) {
  if (!identical(names(left), names(right))) {
    stop(sprintf("%s and %s do not have matching result schemas.", left_label, right_label), call. = FALSE)
  }
  invisible(left)
}

require_output <- function(path, label) {
  output_info <- file.info(path)
  if (!file.exists(path) || is.na(output_info$size) || output_info$size == 0L) {
    stop(sprintf("%s was not written: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

read_result <- function(path, label) {
  require_file(path, label)
  result <- readr::read_tsv(path, show_col_types = FALSE)
  require_nonempty(result, label)
  require_metabolites(result, label)
  result
}

combine_outcome_designs <- function(conservative_file, liberal_file, label) {
  conservative_results <- read_result(conservative_file, sprintf("Conservative %s significant-result summary", label))
  liberal_results <- read_result(liberal_file, sprintf("Liberal %s significant-result summary", label))
  require_matching_schema(
    conservative_results,
    liberal_results,
    sprintf("Conservative %s significant-result summary", label),
    sprintf("Liberal %s significant-result summary", label)
  )

  liberal_results <- liberal_results[!(liberal_results$Metabolite %in% conservative_results$Metabolite), , drop = FALSE]
  rbind(conservative_results, liberal_results)
}

require_matching_outcome_schemas <- function(t2dm_results, fg_results, hba1c_results) {
  t2dm_columns <- setdiff(names(t2dm_results), "Metabolite")
  if (!identical(t2dm_columns, setdiff(names(fg_results), "Metabolite")) ||
      !identical(t2dm_columns, setdiff(names(hba1c_results), "Metabolite"))) {
    stop("Combined outcome summaries do not have matching result schemas.", call. = FALSE)
  }
  invisible(t2dm_results)
}

suffix_outcome_columns <- function(data, outcome) {
  result_columns <- setdiff(names(data), "Metabolite")
  names(data)[match(result_columns, names(data))] <- paste0(result_columns, "_", outcome)
  data
}

merge_outcomes <- function(t2dm_results, fg_results, hba1c_results) {
  require_matching_outcome_schemas(t2dm_results, fg_results, hba1c_results)
  unique_metabolites <- unique(c(t2dm_results$Metabolite, fg_results$Metabolite, hba1c_results$Metabolite))
  full_significant_results <- data.frame(Metabolite = unique_metabolites)
  full_significant_results <- merge(
    full_significant_results,
    suffix_outcome_columns(t2dm_results, "T2DM"),
    by = "Metabolite",
    all.x = TRUE
  )
  full_significant_results <- merge(
    full_significant_results,
    suffix_outcome_columns(fg_results, "FG"),
    by = "Metabolite",
    all.x = TRUE
  )
  merge(
    full_significant_results,
    suffix_outcome_columns(hba1c_results, "HBA1C"),
    by = "Metabolite",
    all.x = TRUE
  )
}

stage_dir <- file.path(paths[["output_dir"]], "04_significance_filtering")
conservative_dir <- file.path(stage_dir, "conservative")
liberal_dir <- file.path(stage_dir, "liberal")
output_dir <- file.path(stage_dir, "combined")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop(sprintf("Cannot create output directory: %s", output_dir), call. = FALSE)

t2dm_results <- combine_outcome_designs(
  file.path(conservative_dir, "significant_T2DM_results.tsv"),
  file.path(liberal_dir, "significant_T2DM_results_liberal.tsv"),
  "type 2 diabetes"
)
fg_results <- combine_outcome_designs(
  file.path(conservative_dir, "significant_FG_results.tsv"),
  file.path(liberal_dir, "significant_FG_results_liberal.tsv"),
  "fasting glucose"
)
hba1c_results <- combine_outcome_designs(
  file.path(conservative_dir, "significant_HBA1C_results.tsv"),
  file.path(liberal_dir, "significant_HBA1C_results_liberal.tsv"),
  "HbA1c"
)

full_significant_results <- merge_outcomes(t2dm_results, fg_results, hba1c_results)
output_file <- file.path(output_dir, "Full_Significant_Results_Manuscript.tsv")
utils::write.table(full_significant_results, file = output_file, sep = "\t", row.names = FALSE)
require_output(output_file, "Combined significant-result summary")
