source("R/utils.R")

paths <- data_paths("METABOLOME_MR_OUTPUT_DIR")

harmonise_matched <- function(data, metabolite, outcome) {
  data <- data[!is.na(data$EAF), , drop = FALSE]
  if (!nrow(data)) return(list(data = data, removed = 0L))

  exposure <- data[, c("SNP", "Beta", "SE", "EffectAllele", "NonEffectAllele", "EAF"), drop = FALSE]
  outcome_data <- data[, c("SNP", "out_Beta", "out_SE", "out_EffectAllele", "out_NonEffectAllele", "out_EAF"), drop = FALSE]
  exposure$EffectAllele <- normalise_allele(exposure$EffectAllele)
  exposure$NonEffectAllele <- normalise_allele(exposure$NonEffectAllele)
  outcome_data$out_EffectAllele <- normalise_allele(outcome_data$out_EffectAllele)
  outcome_data$out_NonEffectAllele <- normalise_allele(outcome_data$out_NonEffectAllele)
  names(exposure) <- c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")
  names(outcome_data) <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome")
  exposure$exposure <- "exposure"
  exposure$id.exposure <- "exposure"
  outcome_data$outcome <- "outcome"
  outcome_data$id.outcome <- "outcome"

  harmonised <- TwoSampleMR::harmonise_data(exposure, outcome_data, action = 2)
  remove <- harmonised$SNP[!harmonised$mr_keep]  # incompatible or ambiguous alleles
  data <- data[!(data$SNP %in% remove), , drop = FALSE]
  harmonised <- harmonised[!(harmonised$SNP %in% remove), , drop = FALSE]
  if (!nrow(data)) return(list(data = data, removed = length(remove)))

  index <- match(data$SNP, harmonised$SNP)
  data$Beta <- harmonised$beta.exposure[index]
  data$SE <- harmonised$se.exposure[index]
  data$EffectAllele <- harmonised$effect_allele.exposure[index]
  data$NonEffectAllele <- harmonised$other_allele.exposure[index]
  data$EAF <- harmonised$eaf.exposure[index]
  data$out_Beta <- harmonised$beta.outcome[index]
  data$out_SE <- harmonised$se.outcome[index]
  data$out_EffectAllele <- harmonised$effect_allele.outcome[index]
  data$out_NonEffectAllele <- harmonised$other_allele.outcome[index]
  data$out_EAF <- harmonised$eaf.outcome[index]
  list(data = data, removed = length(remove))
}

outcomes <- c("T2DM", "FG", "HBA1C")
reverse_mr_dir <- file.path(paths[["output_dir"]], "06_reverse_mr")

matching_root <- file.path(reverse_mr_dir, "matching")
matching <- as.data.frame(readr::read_tsv(file.path(matching_root, "matching_summary.tsv"), show_col_types = FALSE))

matched_root <- file.path(reverse_mr_dir, "matched")
harmonised_root <- file.path(reverse_mr_dir, "harmonised")
for (code in outcomes) dir.create(file.path(harmonised_root, code), recursive = TRUE, showWarnings = FALSE)

summary_rows <- lapply(seq_len(nrow(matching)), function(index) {
  row <- matching[index, , drop = FALSE]
  matched_path <- file.path(matched_root, row$Outcome, row$Matched_File)
  data <- as.data.frame(readr::read_tsv(matched_path, show_col_types = FALSE))
  result <- harmonise_matched(data, row$Metabolite, row$Outcome)
  relative_file <- NA_character_
  status <- "no_harmonised_instruments"
  if (nrow(result$data)) {
    relative_file <- sub("_matched\\.tsv$", "_harmonised.tsv", row$Matched_File)
    output_path <- file.path(harmonised_root, row$Outcome, relative_file)
    readr::write_tsv(result$data, output_path)
    status <- "harmonised"
  }
  data.frame(
    Metabolite = as.character(row$Metabolite),
    Outcome = as.character(row$Outcome),
    Number_Matched_IVs = nrow(data),
    Number_Palindromic_Ambiguous_Removed = result$removed,
    Number_Harmonised_IVs = nrow(result$data),
    Harmonisation_Status = status,
    Harmonised_File = relative_file,
    stringsAsFactors = FALSE
  )
})
harmonisation_summary <- do.call(rbind, summary_rows)
summary_path <- file.path(harmonised_root, "harmonisation_summary.tsv")
readr::write_tsv(harmonisation_summary, summary_path)

