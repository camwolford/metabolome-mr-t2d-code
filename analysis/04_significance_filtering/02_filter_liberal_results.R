source("config/environment.R")

paths <- archive_paths("METABOLOME_MR_OUTPUT_DIR")

significance_threshold <- 2.826e-5

require_file <- function(path, label) {
  if (!file.exists(path)) stop(sprintf("%s is missing: %s", label, path), call. = FALSE)
  invisible(path)
}

require_nonempty <- function(data, label) {
  if (!nrow(data)) stop(sprintf("%s has zero rows.", label), call. = FALSE)
  invisible(data)
}

require_columns <- function(data, columns, label) {
  absent <- setdiff(columns, names(data))
  if (length(absent)) stop(sprintf("%s is missing: %s", label, paste(absent, collapse = ", ")), call. = FALSE)
  invisible(data)
}

require_numeric_columns <- function(data, columns, label) {
  non_numeric <- columns[!vapply(data[columns], is.numeric, logical(1))]
  if (length(non_numeric)) stop(sprintf("%s has non-numeric columns: %s", label, paste(non_numeric, collapse = ", ")), call. = FALSE)
  invisible(data)
}

require_metabolites <- function(data, label) {
  if (anyNA(data$Metabolite) || any(!nzchar(data$Metabolite))) {
    stop(sprintf("%s has missing metabolite identifiers.", label), call. = FALSE)
  }
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

filter_outcome <- function(input_file, output_file, label) {
  require_file(input_file, sprintf("Liberal %s forward-MR summary", label))
  results <- readr::read_tsv(input_file, show_col_types = FALSE)
  require_nonempty(results, sprintf("Liberal %s forward-MR summary", label))
  required_columns <- c(
    "Metabolite", "Number_of_IVs", "Fixed_IVW_Pval", "Random_IVW_Pval",
    "Egger_Pval", "Weighted_Median_Pval", "Weighted_Mode_Pval",
    "Fixed_IVW_Estimate", "Egger_Estimate", "Weighted_Median_Estimate",
    "Weighted_Mode_Estimate"
  )
  require_columns(results, required_columns, sprintf("Liberal %s forward-MR summary", label))
  require_numeric_columns(results, setdiff(required_columns, "Metabolite"), sprintf("Liberal %s forward-MR summary", label))
  require_metabolites(results, sprintf("Liberal %s forward-MR summary", label))

  significant_results <- results[which(
    results$Fixed_IVW_Pval < significance_threshold & results$Number_of_IVs <= 3
  ), , drop = FALSE]
  significant_results <- rbind(
    significant_results,
    results[which(results$Random_IVW_Pval < significance_threshold & results$Number_of_IVs > 3), , drop = FALSE]
  )

  holder <- significant_results[which(
    significant_results$Egger_Pval < 0.05 |
      significant_results$Weighted_Median_Pval < 0.05 |
      significant_results$Weighted_Mode_Pval < 0.05
  ), , drop = FALSE]
  holder <- holder[which(
    ifelse(holder$Egger_Pval < 0.05 & sign(holder$Egger_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE) &
      ifelse(holder$Weighted_Median_Pval < 0.05 & sign(holder$Weighted_Median_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE) &
      ifelse(holder$Weighted_Mode_Pval < 0.05 & sign(holder$Weighted_Mode_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE)
  ), , drop = FALSE]
  significant_results <- significant_results[which(is.na(significant_results$Egger_Pval)), , drop = FALSE]
  significant_results <- rbind(significant_results, holder)

  utils::write.table(significant_results, file = output_file, sep = "\t", row.names = FALSE)
  require_output(output_file, sprintf("Liberal %s significant-result summary", label))
}

forward_mr_dir <- file.path(paths[["output_dir"]], "03_forward_mr", "liberal")
output_dir <- file.path(paths[["output_dir"]], "04_significance_filtering", "liberal")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop(sprintf("Cannot create output directory: %s", output_dir), call. = FALSE)

filter_outcome(
  file.path(forward_mr_dir, "T2DM_MR_Results_Liberal.tsv"),
  file.path(output_dir, "significant_T2DM_results_liberal.tsv"),
  "type 2 diabetes"
)
filter_outcome(
  file.path(forward_mr_dir, "FG_MR_Results_Liberal.tsv"),
  file.path(output_dir, "significant_FG_results_liberal.tsv"),
  "fasting glucose"
)
filter_outcome(
  file.path(forward_mr_dir, "HBA1C_MR_Results_Liberal.tsv"),
  file.path(output_dir, "significant_HBA1C_results_liberal.tsv"),
  "HbA1c"
)
