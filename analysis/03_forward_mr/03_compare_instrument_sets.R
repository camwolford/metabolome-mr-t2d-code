source("R/utils.R")

paths <- data_paths("METABOLOME_MR_OUTPUT_DIR")

compare_outcome <- function(conservative_file, liberal_file, output_file, label) {
  conservative_results <- readr::read_tsv(conservative_file, show_col_types = FALSE)
  liberal_results <- readr::read_tsv(liberal_file, show_col_types = FALSE)

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
}

forward_mr_root <- file.path(paths[["output_dir"]], "03_forward_mr")
conservative_dir <- file.path(forward_mr_root, "conservative")
liberal_dir <- file.path(forward_mr_root, "liberal")
comparison_dir <- file.path(forward_mr_root, "instrument_comparison")
dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)

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
