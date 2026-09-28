source("R/utils.R")

paths <- data_paths(c(
  "METABOLOME_MR_INPUT_DIR",
  "METABOLOME_MR_WORK_DIR",
  "METABOLOME_MR_OUTPUT_DIR"
))

marker_name <- function(chromosome, position, effect_allele, other_allele) {
  paste("chr", chromosome, position, normalise_allele(effect_allele), normalise_allele(other_allele), sep = "_")
}

safe_file_stem <- function(value) {
  gsub("[^[:alnum:]_.-]", "_", value)
}

read_outcome_instruments <- function(path, code) {
  data <- as.data.frame(readr::read_tsv(path, show_col_types = FALSE))
  required <- c("Chromosome", "Position", "SNP", "EffectAllele", "NonEffectAllele", "EAF", "Beta", "SE", "Pval")
  data$MarkerName <- marker_name(data$Chromosome, data$Position, data$EffectAllele, data$NonEffectAllele)
  data$reverseMarkerName <- marker_name(data$Chromosome, data$Position, data$NonEffectAllele, data$EffectAllele)
  data$Fstat <- (data$Beta / data$SE)^2
  data
}

match_instruments <- function(instruments, metabolite_gwas, metabolite, outcome) {
  direct <- match(instruments$MarkerName, metabolite_gwas$MarkerName)
  reversed <- match(instruments$reverseMarkerName, metabolite_gwas$MarkerName)
  matched_index <- ifelse(!is.na(direct), direct, reversed)
  matched <- !is.na(matched_index)

  unmatched <- data.frame(
    Metabolite = rep(metabolite, sum(!matched)),
    Outcome = rep(outcome, sum(!matched)),
    SNP = as.character(instruments$SNP[!matched]),
    Chromosome = instruments$Chromosome[!matched],
    Position = instruments$Position[!matched],
    EffectAllele = as.character(instruments$EffectAllele[!matched]),
    NonEffectAllele = as.character(instruments$NonEffectAllele[!matched]),
    Reason = "no_direct_or_reversed_allele_match",
    stringsAsFactors = FALSE
  )
  if (!any(matched)) return(list(matched = NULL, unmatched = unmatched))

  exposure <- instruments[matched, , drop = FALSE]
  outcome_rows <- metabolite_gwas[matched_index[matched], , drop = FALSE]
  matched_data <- data.frame(
    Metabolite = rep(metabolite, nrow(exposure)),
    Outcome = rep(outcome, nrow(exposure)),
    Chromosome = exposure$Chromosome,
    Position = exposure$Position,
    MarkerName = exposure$MarkerName,
    reverseMarkerName = exposure$reverseMarkerName,
    SNP = as.character(exposure$SNP),
    EffectAllele = normalise_allele(exposure$EffectAllele),
    NonEffectAllele = normalise_allele(exposure$NonEffectAllele),
    EAF = exposure$EAF,
    Beta = exposure$Beta,
    SE = exposure$SE,
    Pval = exposure$Pval,
    Fstat = exposure$Fstat,
    proxy = NA_character_,
    out_Chromosome = outcome_rows$Chromosome,
    out_Position = outcome_rows$Position,
    out_EffectAllele = normalise_allele(outcome_rows$EffectAllele),
    out_NonEffectAllele = normalise_allele(outcome_rows$NonEffectAllele),
    out_EAF = outcome_rows$EAF,
    out_Beta = outcome_rows$Beta,
    out_SE = outcome_rows$SE,
    out_Pval = outcome_rows$Pval,
    stringsAsFactors = FALSE
  )
  list(matched = matched_data, unmatched = unmatched)
}

outcomes <- c("T2DM", "FG", "HBA1C")
output_root <- file.path(paths[["output_dir"]], "06_reverse_mr")
steiger_path <- file.path(output_root, "steiger", "post_steiger_candidates.tsv")
candidates <- as.data.frame(readr::read_tsv(steiger_path, show_col_types = FALSE))
candidate_names <- as.character(candidates$Metabolite)


instrument_root <- file.path(output_root, "outcome_instruments")
outcome_instruments <- stats::setNames(lapply(outcomes, function(code) {
  read_outcome_instruments(file.path(instrument_root, sprintf("%s_clumped.tsv", code)), code)
}), outcomes)

comp_ids <- metabolite_comp_ids(paths[["input_dir"]])
instrument_positions <- unique(unlist(lapply(outcome_instruments, `[[`, "Position")))
matched_root <- file.path(output_root, "matched")
matching_root <- file.path(output_root, "matching")
for (code in outcomes) dir.create(file.path(matched_root, code), recursive = TRUE, showWarnings = FALSE)
dir.create(matching_root, recursive = TRUE, showWarnings = FALSE)

matching_rows <- list()
unmatched_rows <- list()
row_index <- 0L

for (metabolite in candidate_names) {
  # Only the rows at reverse-MR instrument positions are needed from the metabolite GWAS.
  metabolite_gwas <- read_metabolite_gwas(comp_ids[[metabolite]], function(d) d$Position %in% instrument_positions)
  metabolite_gwas$MarkerName <- marker_name(metabolite_gwas$Chromosome, metabolite_gwas$Position, metabolite_gwas$EffectAllele, metabolite_gwas$NonEffectAllele)
  for (code in outcomes) {
    result <- match_instruments(outcome_instruments[[code]], metabolite_gwas, metabolite, code)
    matched_count <- if (is.null(result$matched)) 0L else nrow(result$matched)
    unmatched_count <- nrow(result$unmatched)
    relative_file <- NA_character_
    if (matched_count) {
      relative_file <- paste0(safe_file_stem(metabolite), "_", code, "_matched.tsv")
      matched_path <- file.path(matched_root, code, relative_file)
      readr::write_tsv(result$matched, matched_path)
    }
    if (unmatched_count) unmatched_rows[[length(unmatched_rows) + 1L]] <- result$unmatched
    row_index <- row_index + 1L
    matching_rows[[row_index]] <- data.frame(
      Metabolite = metabolite,
      Outcome = code,
      Number_Requested_IVs = nrow(outcome_instruments[[code]]),
      Number_Direct_Matches = matched_count,
      Number_Unmatched_IVs = unmatched_count,
      Match_Status = if (unmatched_count) "unmatched_instruments" else "matched",
      Matched_File = relative_file,
      stringsAsFactors = FALSE
    )
  }
}

matching_summary <- do.call(rbind, matching_rows)
matching_path <- file.path(matching_root, "matching_summary.tsv")
readr::write_tsv(matching_summary, matching_path)

unmatched <- if (length(unmatched_rows)) {
  do.call(rbind, unmatched_rows)
} else {
  data.frame(
    Metabolite = character(), Outcome = character(), SNP = character(), Chromosome = numeric(), Position = numeric(),
    EffectAllele = character(), NonEffectAllele = character(), Reason = character(), stringsAsFactors = FALSE
  )
}
unmatched_path <- file.path(matching_root, "unmatched_instruments.tsv")
readr::write_tsv(unmatched, unmatched_path)

