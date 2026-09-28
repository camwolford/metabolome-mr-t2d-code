source("R/utils.R")

paths <- data_paths("METABOLOME_MR_OUTPUT_DIR")

result_field <- function(result, field, label) {
  methods::slot(result, field)
}

extract_scalar <- function(result, field, label) {
  value <- result_field(result, field, label)
  as.numeric(value)
}

# Heter.Stat holds Cochran's Q and its p-value; part = 2 returns the p-value.
extract_heterogeneity_statistic <- function(result, label, part = 1L) {
  value <- result_field(result, "Heter.Stat", label)
  statistic <- suppressWarnings(as.numeric(value[part]))
  statistic
}

run_estimator <- function(expression, label) {
  tryCatch(
    expression(),
    error = function(error) stop(sprintf("%s failed: %s", label, conditionMessage(error)), call. = FALSE)
  )
}

run_reverse_mr <- function(data, metabolite, outcome) {

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
    Random_IVW_HetStat_P = extract_heterogeneity_statistic(random_ivw, "Random-effects IVW reverse MR", part = 2L),
    Random_IVW_FStat = extract_scalar(random_ivw, "Fstat", "Random-effects IVW reverse MR"),
    Fixed_IVW_Estimate = extract_scalar(fixed_ivw, "Estimate", "Fixed-effects IVW reverse MR"),
    Fixed_IVW_SE = extract_scalar(fixed_ivw, "StdError", "Fixed-effects IVW reverse MR"),
    Fixed_IVW_Pval = extract_scalar(fixed_ivw, "Pvalue", "Fixed-effects IVW reverse MR"),
    Fixed_IVW_RSE = extract_scalar(fixed_ivw, "RSE", "Fixed-effects IVW reverse MR"),
    Fixed_IVW_HetStat = extract_heterogeneity_statistic(fixed_ivw, "Fixed-effects IVW reverse MR"),
    Fixed_IVW_HetStat_P = extract_heterogeneity_statistic(fixed_ivw, "Fixed-effects IVW reverse MR", part = 2L),
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

outcomes <- c("T2DM", "FG", "HBA1C")
reverse_mr_dir <- file.path(paths[["output_dir"]], "06_reverse_mr")


harmonised_root <- file.path(reverse_mr_dir, "harmonised")
harmonisation <- as.data.frame(readr::read_tsv(file.path(harmonised_root, "harmonisation_summary.tsv"), show_col_types = FALSE))

results <- do.call(rbind, lapply(seq_len(nrow(harmonisation)), function(index) {
  row <- harmonisation[index, , drop = FALSE]
  path <- file.path(harmonised_root, row$Outcome, row$Harmonised_File)
  data <- as.data.frame(readr::read_tsv(path, show_col_types = FALSE))
  run_reverse_mr(data, row$Metabolite, row$Outcome)
}))

results_root <- file.path(reverse_mr_dir, "results")
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
for (code in outcomes) {
  outcome_results <- results[results$Outcome == code, , drop = FALSE]
  outcome_path <- file.path(results_root, sprintf("%s_reverse_mr_raw.tsv", tolower(code)))
  readr::write_tsv(outcome_results, outcome_path)
}
combined_path <- file.path(results_root, "reverse_mr_raw.tsv")
readr::write_tsv(results, combined_path)
