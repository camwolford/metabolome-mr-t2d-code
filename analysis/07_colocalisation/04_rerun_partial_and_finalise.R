source("R/utils.R")

# Where only some instruments colocalise, re-run MR with the colocalising instruments. The association is kept
# if the re-run estimate has the same direction as the original and passes the Bonferroni threshold.
# A metabolite is prioritised if at least one of its associations is kept.
paths <- data_paths("METABOLOME_MR_OUTPUT_DIR")
sensitivity_dir <- file.path(paths[["output_dir"]], "05_sensitivity_analysis")
output_dir <- file.path(paths[["output_dir"]], "07_colocalisation")
prefixes <- c(T2DM = "t2dm", FG = "fg", HBA1C = "hba1c")

loci <- as.data.frame(readr::read_tsv(file.path(output_dir, "coloc_loci.tsv"), show_col_types = FALSE))
associations <- as.data.frame(readr::read_tsv(file.path(output_dir, "coloc_associations.tsv"), show_col_types = FALSE))

# Wald ratio (first-order SE) for one instrument; fixed-effect IVW for two or three; random-effects IVW for four or more.
run_mr <- function(bx, bxse, by, byse) {
  if (length(bx) == 1) {
    estimate <- by / bx
    se <- byse / abs(bx)
    return(c(Estimate = estimate, SE = se, Pval = 2 * stats::pnorm(-abs(estimate / se))))
  }
  fit <- MendelianRandomization::mr_ivw(MendelianRandomization::mr_input(bx = bx, bxse = bxse, by = by, byse = byse),
                                       model = if (length(bx) <= 3) "fixed" else "random")
  c(Estimate = fit$Estimate, SE = fit$StdError, Pval = fit$Pvalue)
}

original_estimate <- function(metabolite, outcome, design) {
  original <- readr::read_tsv(filtered_results_path(sensitivity_dir, design, outcome), show_col_types = FALSE)
  original <- original[original$Metabolite == metabolite, ]
  if (original$Number_of_IVs <= 3) original$Fixed_IVW_Estimate else original$Random_IVW_Estimate
}

reruns <- list()
associations$Retained <- associations$Colocalisation == "colocalised"
for (i in which(associations$Colocalisation == "partial")) {
  a <- associations[i, ]
  prefix <- prefixes[[a$Outcome]]
  colocalising <- loci$SNP[loci$Metabolite == a$Metabolite & loci$Outcome == a$Outcome & loci$Colocalised]
  instruments <- readr::read_tsv(instrument_files(harmonised_directory(paths[["output_dir"]], a$Design, a$Outcome))[[a$Metabolite]], show_col_types = FALSE)
  instruments <- instruments[instruments$SNP %in% colocalising, ]
  result <- run_mr(instruments$Beta, instruments$SE, instruments[[paste0(prefix, "_Beta")]], instruments[[paste0(prefix, "_SE")]])
  same_direction <- sign(result[["Estimate"]]) == sign(original_estimate(a$Metabolite, a$Outcome, a$Design))
  associations$Retained[i] <- same_direction && result[["Pval"]] <= significance_threshold
  reruns[[length(reruns) + 1]] <- data.frame(
    Metabolite = a$Metabolite, Outcome = a$Outcome, Design = a$Design, Number_of_IVs = nrow(instruments),
    Estimate = result[["Estimate"]], SE = result[["SE"]], Pval = result[["Pval"]],
    Same_Direction = same_direction, Retained = associations$Retained[i]
  )
}

readr::write_tsv(do.call(rbind, reruns), file.path(output_dir, "rerun_mr_results.tsv"))
readr::write_tsv(associations, file.path(output_dir, "final_associations.tsv"))
readr::write_tsv(data.frame(Metabolite = unique(associations$Metabolite[associations$Retained])),
                 file.path(output_dir, "prioritised_metabolites.tsv"))
