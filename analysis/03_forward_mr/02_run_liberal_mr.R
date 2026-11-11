source("R/utils.R")

paths <- data_paths("METABOLOME_MR_OUTPUT_DIR")

new_results <- function(number_of_metabolites) {
  data.frame(
    Metabolite = rep(NA_character_, number_of_metabolites),
    Number_of_IVs = rep(NA_integer_, number_of_metabolites),
    Number_of_Proxies = rep(NA_integer_, number_of_metabolites),
    Weighted_Mode_Estimate = rep(NA_real_, number_of_metabolites),
    Weighted_Mode_SE = rep(NA_real_, number_of_metabolites),
    Weighted_Mode_Pval = rep(NA_real_, number_of_metabolites),
    Weighted_Median_Estimate = rep(NA_real_, number_of_metabolites),
    Weighted_Median_SE = rep(NA_real_, number_of_metabolites),
    Weighted_Median_Pval = rep(NA_real_, number_of_metabolites),
    Random_IVW_Estimate = rep(NA_real_, number_of_metabolites),
    Random_IVW_SE = rep(NA_real_, number_of_metabolites),
    Random_IVW_Pval = rep(NA_real_, number_of_metabolites),
    Random_IVW_RSE = rep(NA_real_, number_of_metabolites),
    Random_IVW_HetStat = rep(NA_real_, number_of_metabolites),
    Random_IVW_HetStat_P = rep(NA_real_, number_of_metabolites),
    Random_IVW_FStat = rep(NA_real_, number_of_metabolites),
    Fixed_IVW_Estimate = rep(NA_real_, number_of_metabolites),
    Fixed_IVW_SE = rep(NA_real_, number_of_metabolites),
    Fixed_IVW_Pval = rep(NA_real_, number_of_metabolites),
    Fixed_IVW_RSE = rep(NA_real_, number_of_metabolites),
    Fixed_IVW_HetStat = rep(NA_real_, number_of_metabolites),
    Fixed_IVW_HetStat_P = rep(NA_real_, number_of_metabolites),
    Fixed_IVW_FStat = rep(NA_real_, number_of_metabolites),
    Egger_Estimate = rep(NA_real_, number_of_metabolites),
    Egger_SE = rep(NA_real_, number_of_metabolites),
    Egger_Pval = rep(NA_real_, number_of_metabolites),
    Egger_Intercept = rep(NA_real_, number_of_metabolites),
    Egger_Intercept_SE = rep(NA_real_, number_of_metabolites),
    Egger_Intercept_Pval = rep(NA_real_, number_of_metabolites),
    Egger_RSE = rep(NA_real_, number_of_metabolites),
    Egger_HetStat = rep(NA_real_, number_of_metabolites),
    Egger_HetStat_P = rep(NA_real_, number_of_metabolites),
    Egger_Isq = rep(NA_real_, number_of_metabolites)
  )
}

run_outcome_mr <- function(input_dir, outcome_beta, outcome_se, output_file) {
  metabolite_files <- list.files(input_dir, full.names = TRUE)

  results_df <- new_results(length(metabolite_files))

  for (metabolite_counter in seq_along(metabolite_files)) {
    metabolite_file <- metabolite_files[[metabolite_counter]]
    metabolite_data <- readr::read_tsv(metabolite_file, show_col_types = FALSE)

    results_df$Metabolite[metabolite_counter] <- metabolite_data$Metabolite[1]
    results_df$Number_of_IVs[metabolite_counter] <- nrow(metabolite_data)
    results_df$Number_of_Proxies[metabolite_counter] <- sum(!is.na(metabolite_data$proxy))

    if (nrow(metabolite_data) > 2) {
      MRObject <- MendelianRandomization::mr_input(
        bx = metabolite_data$Beta,
        bxse = metabolite_data$SE,
        by = metabolite_data[[outcome_beta]],
        byse = metabolite_data[[outcome_se]]
      )
      MR_weighted_mode_out <- MendelianRandomization::mr_mbe(MRObject, weighting = "weighted")
      MR_weighted_median_out <- MendelianRandomization::mr_median(MRObject, weighting = "weighted")
      MR_random_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "random")
      MR_fixed_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "fixed")
      MR_egger_out <- MendelianRandomization::mr_egger(MRObject)

      results_df$Weighted_Mode_Estimate[metabolite_counter] <- MR_weighted_mode_out$Estimate
      results_df$Weighted_Mode_SE[metabolite_counter] <- MR_weighted_mode_out$StdError
      results_df$Weighted_Mode_Pval[metabolite_counter] <- MR_weighted_mode_out$Pvalue
      results_df$Weighted_Median_Estimate[metabolite_counter] <- MR_weighted_median_out$Estimate
      results_df$Weighted_Median_SE[metabolite_counter] <- MR_weighted_median_out$StdError
      results_df$Weighted_Median_Pval[metabolite_counter] <- MR_weighted_median_out$Pvalue
      results_df$Random_IVW_Estimate[metabolite_counter] <- MR_random_IVW_out$Estimate
      results_df$Random_IVW_SE[metabolite_counter] <- MR_random_IVW_out$StdError
      results_df$Random_IVW_Pval[metabolite_counter] <- MR_random_IVW_out$Pvalue
      results_df$Random_IVW_RSE[metabolite_counter] <- MR_random_IVW_out$RSE
      results_df$Random_IVW_HetStat[metabolite_counter] <- MR_random_IVW_out$Heter.Stat[1]
      results_df$Random_IVW_HetStat_P[metabolite_counter] <- MR_random_IVW_out$Heter.Stat[2]
      results_df$Random_IVW_FStat[metabolite_counter] <- MR_random_IVW_out$Fstat
      results_df$Fixed_IVW_Estimate[metabolite_counter] <- MR_fixed_IVW_out$Estimate
      results_df$Fixed_IVW_SE[metabolite_counter] <- MR_fixed_IVW_out$StdError
      results_df$Fixed_IVW_Pval[metabolite_counter] <- MR_fixed_IVW_out$Pvalue
      results_df$Fixed_IVW_RSE[metabolite_counter] <- MR_fixed_IVW_out$RSE
      results_df$Fixed_IVW_HetStat[metabolite_counter] <- MR_fixed_IVW_out$Heter.Stat[1]
      results_df$Fixed_IVW_HetStat_P[metabolite_counter] <- MR_fixed_IVW_out$Heter.Stat[2]
      results_df$Fixed_IVW_FStat[metabolite_counter] <- MR_fixed_IVW_out$Fstat
      results_df$Egger_Estimate[metabolite_counter] <- MR_egger_out$Estimate
      results_df$Egger_SE[metabolite_counter] <- MR_egger_out$StdError.Est
      results_df$Egger_Pval[metabolite_counter] <- MR_egger_out$Pvalue.Est
      results_df$Egger_Intercept[metabolite_counter] <- MR_egger_out$Intercept
      results_df$Egger_Intercept_SE[metabolite_counter] <- MR_egger_out$StdError.Int
      results_df$Egger_Intercept_Pval[metabolite_counter] <- MR_egger_out$Pvalue.Int
      results_df$Egger_RSE[metabolite_counter] <- MR_egger_out$RSE
      results_df$Egger_HetStat[metabolite_counter] <- MR_egger_out$Heter.Stat[1]
      results_df$Egger_HetStat_P[metabolite_counter] <- MR_egger_out$Heter.Stat[2]
      results_df$Egger_Isq[metabolite_counter] <- MR_egger_out$I.sq
    }

    if (nrow(metabolite_data) == 2) {
      MRObject <- MendelianRandomization::mr_input(
        bx = metabolite_data$Beta,
        bxse = metabolite_data$SE,
        by = metabolite_data[[outcome_beta]],
        byse = metabolite_data[[outcome_se]]
      )
      MR_random_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "random")
      MR_fixed_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "fixed")

      results_df$Random_IVW_Estimate[metabolite_counter] <- MR_random_IVW_out$Estimate
      results_df$Random_IVW_SE[metabolite_counter] <- MR_random_IVW_out$StdError
      results_df$Random_IVW_Pval[metabolite_counter] <- MR_random_IVW_out$Pvalue
      results_df$Random_IVW_RSE[metabolite_counter] <- MR_random_IVW_out$RSE
      results_df$Random_IVW_HetStat[metabolite_counter] <- MR_random_IVW_out$Heter.Stat[1]
      results_df$Random_IVW_HetStat_P[metabolite_counter] <- MR_random_IVW_out$Heter.Stat[2]
      results_df$Random_IVW_FStat[metabolite_counter] <- MR_random_IVW_out$Fstat
      results_df$Fixed_IVW_Estimate[metabolite_counter] <- MR_fixed_IVW_out$Estimate
      results_df$Fixed_IVW_SE[metabolite_counter] <- MR_fixed_IVW_out$StdError
      results_df$Fixed_IVW_Pval[metabolite_counter] <- MR_fixed_IVW_out$Pvalue
      results_df$Fixed_IVW_RSE[metabolite_counter] <- MR_fixed_IVW_out$RSE
      results_df$Fixed_IVW_HetStat[metabolite_counter] <- MR_fixed_IVW_out$Heter.Stat[1]
      results_df$Fixed_IVW_HetStat_P[metabolite_counter] <- MR_fixed_IVW_out$Heter.Stat[2]
      results_df$Fixed_IVW_FStat[metabolite_counter] <- MR_fixed_IVW_out$Fstat
    }

    if (nrow(metabolite_data) == 1) {
      MRObject <- MendelianRandomization::mr_input(
        bx = metabolite_data$Beta,
        bxse = metabolite_data$SE,
        by = metabolite_data[[outcome_beta]],
        byse = metabolite_data[[outcome_se]]
      )
      MR_WR_out <- MendelianRandomization::mr_ivw(MRObject)
      results_df$Fixed_IVW_Estimate[metabolite_counter] <- MR_WR_out$Estimate
      results_df$Fixed_IVW_SE[metabolite_counter] <- MR_WR_out$StdError
      results_df$Fixed_IVW_Pval[metabolite_counter] <- MR_WR_out$Pvalue
      results_df$Fixed_IVW_FStat[metabolite_counter] <- MR_WR_out$Fstat
    }
  }

  utils::write.table(results_df, file = output_file, sep = "\t", row.names = FALSE)
  results_df
}

merge_outcome_results <- function(t2dm_results, fg_results, hba1c_results, output_file) {
  unique_metabolites <- unique(c(t2dm_results$Metabolite, fg_results$Metabolite, hba1c_results$Metabolite))
  full_results_df <- data.frame(Metabolite = unique_metabolites)
  full_results_df <- merge(full_results_df, t2dm_results, by = "Metabolite", all.x = TRUE)
  full_results_df <- merge(full_results_df, fg_results, by = "Metabolite", all.x = TRUE)
  full_results_df <- merge(full_results_df, hba1c_results, by = "Metabolite", all.x = TRUE)
  colnames(full_results_df) <- gsub("\\.x", "_T2DM", colnames(full_results_df))
  colnames(full_results_df) <- gsub("\\.y", "_FG", colnames(full_results_df))
  colnames(full_results_df)[60:ncol(full_results_df)] <- paste(colnames(full_results_df)[60:ncol(full_results_df)], "_HBA1C", sep = "")
  utils::write.table(full_results_df, file = output_file, sep = "\t", row.names = FALSE)
}

set.seed(20260706)

harmonised_root <- file.path(paths[["output_dir"]], "02_instrument_selection", "liberal")
output_dir <- file.path(paths[["output_dir"]], "03_forward_mr", "liberal")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

t2dm_results <- run_outcome_mr(
  file.path(harmonised_root, "Harmonised_T2DM_IVs_Liberal"),
  "t2dm_Beta",
  "t2dm_SE",
  file.path(output_dir, "T2DM_MR_Results_Liberal.tsv")
)
fg_results <- run_outcome_mr(
  file.path(harmonised_root, "Harmonised_FG_IVs_Liberal"),
  "fg_Beta",
  "fg_SE",
  file.path(output_dir, "FG_MR_Results_Liberal.tsv")
)
hba1c_results <- run_outcome_mr(
  file.path(harmonised_root, "Harmonised_HBA1C_IVs_Liberal"),
  "hba1c_Beta",
  "hba1c_SE",
  file.path(output_dir, "HBA1C_MR_Results_Liberal.tsv")
)
merge_outcome_results(
  t2dm_results,
  fg_results,
  hba1c_results,
  file.path(output_dir, "Full_MR_Results_Liberal.tsv")
)
