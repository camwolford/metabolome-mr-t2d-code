source("R/utils.R")

paths <- data_paths("METABOLOME_MR_OUTPUT_DIR")

read_result <- function(path, label) {
  result <- readr::read_tsv(path, show_col_types = FALSE)
  result
}

combine_outcome_designs <- function(conservative_file, liberal_file, label) {
  conservative_results <- read_result(conservative_file, sprintf("Conservative %s significant-result summary", label))
  liberal_results <- read_result(liberal_file, sprintf("Liberal %s significant-result summary", label))

  liberal_results <- liberal_results[!(liberal_results$Metabolite %in% conservative_results$Metabolite), , drop = FALSE]
  rbind(conservative_results, liberal_results)
}

suffix_outcome_columns <- function(data, outcome) {
  result_columns <- setdiff(names(data), "Metabolite")
  names(data)[match(result_columns, names(data))] <- paste0(result_columns, "_", outcome)
  data
}

merge_outcomes <- function(t2dm_results, fg_results, hba1c_results) {
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

significance_dir <- file.path(paths[["output_dir"]], "04_significance_filtering")
conservative_dir <- file.path(significance_dir, "conservative")
liberal_dir <- file.path(significance_dir, "liberal")
output_dir <- file.path(significance_dir, "combined")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
