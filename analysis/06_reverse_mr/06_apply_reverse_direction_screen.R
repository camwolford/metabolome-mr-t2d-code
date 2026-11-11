source("R/utils.R")

# Exclude a metabolite if Steiger filtering indicates the reverse direction for any tested outcome, or if
# reverse MR (IVW or MR-Egger) is nominally significant (p < 0.05) for any outcome.
paths <- data_paths("METABOLOME_MR_OUTPUT_DIR")
sensitivity_dir <- file.path(paths[["output_dir"]], "05_sensitivity_analysis")
reverse_mr_dir <- file.path(paths[["output_dir"]], "06_reverse_mr")

candidates <- readr::read_tsv(file.path(sensitivity_dir, "combined", "Full_Filtered_Results_Manuscript.tsv"), show_col_types = FALSE)$Metabolite
steiger <- readr::read_tsv(file.path(reverse_mr_dir, "steiger", "steiger_results.tsv"), show_col_types = FALSE)
reverse <- readr::read_tsv(file.path(reverse_mr_dir, "results", "reverse_mr_raw.tsv"), show_col_types = FALSE)

steiger_flagged <- steiger$Metabolite[which(!as.logical(steiger$Direction_Flag))]
reverse_ivw_p <- ifelse(reverse$Number_of_IVs > 3, reverse$Random_IVW_Pval, reverse$Fixed_IVW_Pval)
reverse_flagged <- reverse$Metabolite[which(reverse_ivw_p < 0.05 | reverse$Egger_Pval < 0.05)]

screen <- data.frame(
  Metabolite = candidates,
  Steiger_Flagged = candidates %in% steiger_flagged,
  Reverse_MR_Flagged = candidates %in% reverse_flagged
)
screen$Retained <- !screen$Steiger_Flagged & !screen$Reverse_MR_Flagged
readr::write_tsv(screen, file.path(reverse_mr_dir, "reverse_direction_screen.tsv"))
