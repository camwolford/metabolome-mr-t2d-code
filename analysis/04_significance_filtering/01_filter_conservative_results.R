source("R/utils.R")

paths <- data_paths("METABOLOME_MR_OUTPUT_DIR")


filter_outcome <- function(input_file, output_file, label) {
  results <- readr::read_tsv(input_file, show_col_types = FALSE)

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
}

forward_mr_dir <- file.path(paths[["output_dir"]], "03_forward_mr", "conservative")
output_dir <- file.path(paths[["output_dir"]], "04_significance_filtering", "conservative")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

filter_outcome(
  file.path(forward_mr_dir, "T2DM_MR_Results.tsv"),
  file.path(output_dir, "significant_T2DM_results.tsv"),
  "type 2 diabetes"
)
filter_outcome(
  file.path(forward_mr_dir, "FG_MR_Results.tsv"),
  file.path(output_dir, "significant_FG_results.tsv"),
  "fasting glucose"
)
filter_outcome(
  file.path(forward_mr_dir, "HBA1C_MR_Results.tsv"),
  file.path(output_dir, "significant_HBA1C_results.tsv"),
  "HbA1c"
)
