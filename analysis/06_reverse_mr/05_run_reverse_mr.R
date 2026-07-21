source("config/environment.R")

paths <- archive_paths("METABOLOME_MR_OUTPUT_DIR")

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

require_unique_keys <- function(data, keys, label) {
  duplicated_rows <- duplicated(data[keys])
  if (any(duplicated_rows)) {
    stop(sprintf("%s has duplicate keys: %s", label, paste(keys, collapse = ", ")), call. = FALSE)
  }
  invisible(data)
}

require_output <- function(path, label) {
  info <- file.info(path)
  if (!file.exists(path) || is.na(info$size) || info$size == 0L) {
    stop(sprintf("%s was not written: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

read_tsv <- function(path, label, allow_empty = FALSE) {
  require_file(path, label)
  data <- as.data.frame(readr::read_tsv(path, show_col_types = FALSE))
  if (!allow_empty) require_nonempty(data, label)
  data
}

key_vector <- function(data, keys) {
  do.call(paste, c(data[keys], sep = "\r"))
}

require_same_keys <- function(left, right, keys, left_label, right_label) {
  if (!setequal(key_vector(left, keys), key_vector(right, keys))) {
    stop(sprintf("%s and %s have different %s memberships.", left_label, right_label, paste(keys, collapse = ", ")), call. = FALSE)
  }
  invisible(left)
}

result_field <- function(result, field, label) {
  if (!isS4(result)) stop(sprintf("%s is not an S4 result.", label), call. = FALSE)
  if (!(field %in% methods::slotNames(result))) {
    stop(sprintf("%s does not contain the %s slot.", label, field), call. = FALSE)
  }
  methods::slot(result, field)
}

extract_scalar <- function(result, field, label) {
  value <- result_field(result, field, label)
  if (is.null(value) || length(value) != 1L) {
    stop(sprintf("%s did not return %s.", label, field), call. = FALSE)
  }
  as.numeric(value)
}

extract_heterogeneity_statistic <- function(result, label) {
  value <- result_field(result, "Heter.Stat", label)
  if (is.null(value) || length(value) != 2L) {
    stop(sprintf("%s did not return a two-value Heter.Stat result.", label), call. = FALSE)
  }
  statistic <- suppressWarnings(as.numeric(value[1]))
  if (length(statistic) != 1L || !is.finite(statistic)) {
    stop(sprintf("%s did not return a valid Heter.Stat statistic.", label), call. = FALSE)
  }
  statistic
}

run_estimator <- function(expression, label) {
  tryCatch(
    expression(),
    error = function(error) stop(sprintf("%s failed: %s", label, conditionMessage(error)), call. = FALSE)
  )
}

run_reverse_mr <- function(data, metabolite, outcome) {
  required <- c("Metabolite", "Outcome", "SNP", "Beta", "SE", "out_Beta", "out_SE", "proxy")
  require_columns(data, required, sprintf("Harmonised %s instruments for %s", outcome, metabolite))
  require_unique_keys(data, "SNP", sprintf("Harmonised %s instruments for %s", outcome, metabolite))
  if (!all(as.character(data$Metabolite) == metabolite) || !all(as.character(data$Outcome) == outcome)) {
    stop(sprintf("Harmonised file has an unexpected metabolite or outcome: %s", metabolite), call. = FALSE)
  }
  if (any(!is.na(data$proxy))) {
    stop(sprintf("Direct-match-only reverse MR found a proxy for %s (%s).", metabolite, outcome), call. = FALSE)
  }
  numeric_columns <- c("Beta", "SE", "out_Beta", "out_SE")
  if (any(vapply(numeric_columns, function(column) !is.numeric(data[[column]]) || any(!is.finite(data[[column]])), logical(1)))) {
    stop(sprintf("Harmonised %s instruments for %s contain invalid effect estimates.", outcome, metabolite), call. = FALSE)
  }
  if (nrow(data) < 3L) {
    stop(sprintf("Raw reverse MR requires at least three instruments for %s (%s) because its source output includes weighted mode, weighted median and MR-Egger.", metabolite, outcome), call. = FALSE)
  }

  mr_input <- MendelianRandomization::mr_input(
    bx = data$Beta,
    bxse = data$SE,
    by = data$out_Beta,
    byse = data$out_SE
  )
  weighted_mode <- run_estimator(
    function() MendelianRandomization::mr_mbe(mr_input, weighting = "weighted"),
    sprintf("Weighted-mode reverse MR for %s (%s)", metabolite, outcome)
  )
  weighted_median <- run_estimator(
    function() MendelianRandomization::mr_median(mr_input, weighting = "weighted"),
    sprintf("Weighted-median reverse MR for %s (%s)", metabolite, outcome)
  )
  random_ivw <- run_estimator(
    function() MendelianRandomization::mr_ivw(mr_input, model = "random"),
    sprintf("Random-effects IVW reverse MR for %s (%s)", metabolite, outcome)
  )
  fixed_ivw <- run_estimator(
    function() MendelianRandomization::mr_ivw(mr_input, model = "fixed"),
    sprintf("Fixed-effects IVW reverse MR for %s (%s)", metabolite, outcome)
  )
  egger <- run_estimator(
    function() MendelianRandomization::mr_egger(mr_input),
    sprintf("MR-Egger reverse MR for %s (%s)", metabolite, outcome)
  )

  data.frame(
    Metabolite = metabolite,
    Outcome = outcome,
    Number_of_IVs = nrow(data),
    Number_of_Proxies = 0L,
    Weighted_Mode_Estimate = extract_scalar(weighted_mode, "Estimate", "Weighted-mode reverse MR"),
    Weighted_Mode_SE = extract_scalar(weighted_mode, "StdError", "Weighted-mode reverse MR"),
    Weighted_Mode_Pval = extract_scalar(weighted_mode, "Pvalue", "Weighted-mode reverse MR"),
    Weighted_Median_Estimate = extract_scalar(weighted_median, "Estimate", "Weighted-median reverse MR"),
    Weighted_Median_SE = extract_scalar(weighted_median, "StdError", "Weighted-median reverse MR"),
    Weighted_Median_Pval = extract_scalar(weighted_median, "Pvalue", "Weighted-median reverse MR"),
    Random_IVW_Estimate = extract_scalar(random_ivw, "Estimate", "Random-effects IVW reverse MR"),
    Random_IVW_SE = extract_scalar(random_ivw, "StdError", "Random-effects IVW reverse MR"),
    Random_IVW_Pval = extract_scalar(random_ivw, "Pvalue", "Random-effects IVW reverse MR"),
    Random_IVW_RSE = extract_scalar(random_ivw, "RSE", "Random-effects IVW reverse MR"),
    Random_IVW_HetStat = extract_heterogeneity_statistic(random_ivw, "Random-effects IVW reverse MR"),
    Random_IVW_FStat = extract_scalar(random_ivw, "Fstat", "Random-effects IVW reverse MR"),
    Fixed_IVW_Estimate = extract_scalar(fixed_ivw, "Estimate", "Fixed-effects IVW reverse MR"),
    Fixed_IVW_SE = extract_scalar(fixed_ivw, "StdError", "Fixed-effects IVW reverse MR"),
    Fixed_IVW_Pval = extract_scalar(fixed_ivw, "Pvalue", "Fixed-effects IVW reverse MR"),
    Fixed_IVW_RSE = extract_scalar(fixed_ivw, "RSE", "Fixed-effects IVW reverse MR"),
    Fixed_IVW_HetStat = extract_heterogeneity_statistic(fixed_ivw, "Fixed-effects IVW reverse MR"),
    Fixed_IVW_FStat = extract_scalar(fixed_ivw, "Fstat", "Fixed-effects IVW reverse MR"),
    Egger_Estimate = extract_scalar(egger, "Estimate", "MR-Egger reverse MR"),
    Egger_SE = extract_scalar(egger, "StdError.Est", "MR-Egger reverse MR"),
    Egger_Pval = extract_scalar(egger, "Pvalue.Est", "MR-Egger reverse MR"),
    Egger_Intercept = extract_scalar(egger, "Intercept", "MR-Egger reverse MR"),
    Egger_Intercept_SE = extract_scalar(egger, "StdError.Int", "MR-Egger reverse MR"),
    Egger_Intercept_Pval = extract_scalar(egger, "Pvalue.Int", "MR-Egger reverse MR"),
    Egger_RSE = extract_scalar(egger, "RSE", "MR-Egger reverse MR"),
    Egger_HetStat = extract_heterogeneity_statistic(egger, "MR-Egger reverse MR"),
    Egger_Isq = extract_scalar(egger, "I.sq", "MR-Egger reverse MR"),
    stringsAsFactors = FALSE
  )
}

require_namespace("readr")
require_namespace("MendelianRandomization")

outcomes <- c("T2DM", "FG", "HBA1C")
stage_root <- file.path(paths[["output_dir"]], "06_reverse_mr")
steiger_candidates <- read_tsv(file.path(stage_root, "steiger", "post_steiger_candidates.tsv"), "Post-Steiger candidates")
require_columns(steiger_candidates, c("Metabolite", "Retained_after_Steiger"), "Post-Steiger candidates")
require_unique_keys(steiger_candidates, "Metabolite", "Post-Steiger candidates")
if (nrow(steiger_candidates) != 54L) {
  stop(sprintf("Post-Steiger candidates must contain 54 metabolites, found %d.", nrow(steiger_candidates)), call. = FALSE)
}
if (anyNA(as.logical(steiger_candidates$Retained_after_Steiger))) stop("Post-Steiger status contains missing values.", call. = FALSE)
candidate_names <- as.character(steiger_candidates$Metabolite)
expected_grid <- expand.grid(Metabolite = candidate_names, Outcome = outcomes, stringsAsFactors = FALSE)
require_unique_keys(expected_grid, c("Metabolite", "Outcome"), "Expected reverse-MR grid")
if (nrow(expected_grid) != 162L) stop("Expected reverse-MR grid must contain 162 metabolite--outcome pairs.", call. = FALSE)

matching <- read_tsv(file.path(stage_root, "matching", "matching_summary.tsv"), "Reverse-MR matching summary")
require_columns(matching, c("Metabolite", "Outcome", "Match_Status"), "Reverse-MR matching summary")
require_unique_keys(matching, c("Metabolite", "Outcome"), "Reverse-MR matching summary")
require_same_keys(matching, expected_grid, c("Metabolite", "Outcome"), "Reverse-MR matching summary", "Expected reverse-MR grid")
if (any(matching$Match_Status != "matched")) stop("Reverse-MR matching is incomplete.", call. = FALSE)

harmonised_root <- file.path(stage_root, "harmonised")
harmonisation <- read_tsv(file.path(harmonised_root, "harmonisation_summary.tsv"), "Reverse-MR harmonisation summary")
require_columns(harmonisation, c("Metabolite", "Outcome", "Number_Harmonised_IVs", "Harmonisation_Status", "Harmonised_File"), "Reverse-MR harmonisation summary")
require_unique_keys(harmonisation, c("Metabolite", "Outcome"), "Reverse-MR harmonisation summary")
require_same_keys(harmonisation, expected_grid, c("Metabolite", "Outcome"), "Reverse-MR harmonisation summary", "Expected reverse-MR grid")
if (any(harmonisation$Harmonisation_Status != "harmonised") || anyNA(harmonisation$Harmonised_File)) {
  stop("Reverse-MR harmonisation is incomplete.", call. = FALSE)
}

results <- do.call(rbind, lapply(seq_len(nrow(harmonisation)), function(index) {
  row <- harmonisation[index, , drop = FALSE]
  path <- file.path(harmonised_root, row$Outcome, row$Harmonised_File)
  data <- read_tsv(path, sprintf("Harmonised %s instruments for %s", row$Outcome, row$Metabolite))
  if (nrow(data) != row$Number_Harmonised_IVs) {
    stop(sprintf("Harmonised instrument count disagrees with the summary for %s (%s).", row$Metabolite, row$Outcome), call. = FALSE)
  }
  run_reverse_mr(data, row$Metabolite, row$Outcome)
}))
require_nonempty(results, "Reverse-MR results")
require_unique_keys(results, c("Metabolite", "Outcome"), "Reverse-MR results")
require_same_keys(results, expected_grid, c("Metabolite", "Outcome"), "Reverse-MR results", "Expected reverse-MR grid")

results_root <- file.path(stage_root, "results")
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(results_root)) stop(sprintf("Cannot create output directory: %s", results_root), call. = FALSE)
for (code in outcomes) {
  outcome_results <- results[results$Outcome == code, , drop = FALSE]
  outcome_path <- file.path(results_root, sprintf("%s_reverse_mr_raw.tsv", tolower(code)))
  readr::write_tsv(outcome_results, outcome_path)
  require_output(outcome_path, sprintf("%s reverse-MR results", code))
}
combined_path <- file.path(results_root, "reverse_mr_raw.tsv")
readr::write_tsv(results, combined_path)
require_output(combined_path, "Combined reverse-MR results")
