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

normalise_allele <- function(values) {
  values <- as.character(values)
  values[values %in% c("TRUE", "True", "true")] <- "T"
  toupper(values)
}

harmonise_matched <- function(data, metabolite, outcome) {
  required <- c(
    "SNP", "Beta", "SE", "EffectAllele", "NonEffectAllele", "EAF",
    "out_Beta", "out_SE", "out_EffectAllele", "out_NonEffectAllele", "out_EAF"
  )
  require_columns(data, required, sprintf("Matched %s instruments for %s", outcome, metabolite))
  require_unique_keys(data, "SNP", sprintf("Matched %s instruments for %s", outcome, metabolite))
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
  required_harmonised <- c(
    "SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure",
    "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome",
    "palindromic", "ambiguous"
  )
  require_columns(harmonised, required_harmonised, sprintf("Harmonised %s instruments for %s", outcome, metabolite))
  remove <- harmonised$SNP[harmonised$palindromic & harmonised$ambiguous]
  data <- data[!(data$SNP %in% remove), , drop = FALSE]
  harmonised <- harmonised[!(harmonised$SNP %in% remove), , drop = FALSE]
  if (!nrow(data)) return(list(data = data, removed = length(remove)))

  index <- match(data$SNP, harmonised$SNP)
  if (anyNA(index)) {
    stop(sprintf("Harmonisation did not return every retained instrument for %s (%s).", metabolite, outcome), call. = FALSE)
  }
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

require_namespace("readr")
require_namespace("TwoSampleMR")

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

matching_root <- file.path(stage_root, "matching")
matching <- read_tsv(file.path(matching_root, "matching_summary.tsv"), "Reverse-MR matching summary")
require_columns(
  matching,
  c("Metabolite", "Outcome", "Number_Requested_IVs", "Number_Direct_Matches", "Number_Unmatched_IVs", "Match_Status", "Matched_File"),
  "Reverse-MR matching summary"
)
require_unique_keys(matching, c("Metabolite", "Outcome"), "Reverse-MR matching summary")
require_same_keys(matching, expected_grid, c("Metabolite", "Outcome"), "Reverse-MR matching summary", "Expected reverse-MR grid")
if (any(matching$Match_Status != "matched") || any(matching$Number_Unmatched_IVs != 0L) || anyNA(matching$Matched_File)) {
  stop("Reverse-MR matching is incomplete; inspect the external matching log before harmonising.", call. = FALSE)
}

matched_root <- file.path(stage_root, "matched")
harmonised_root <- file.path(stage_root, "harmonised")
for (code in outcomes) dir.create(file.path(harmonised_root, code), recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(harmonised_root)) stop(sprintf("Cannot create output directory: %s", harmonised_root), call. = FALSE)

summary_rows <- lapply(seq_len(nrow(matching)), function(index) {
  row <- matching[index, , drop = FALSE]
  matched_path <- file.path(matched_root, row$Outcome, row$Matched_File)
  data <- read_tsv(matched_path, sprintf("Matched %s instruments for %s", row$Outcome, row$Metabolite))
  result <- harmonise_matched(data, row$Metabolite, row$Outcome)
  relative_file <- NA_character_
  status <- "no_harmonised_instruments"
  if (nrow(result$data)) {
    relative_file <- sub("_matched\\.tsv$", "_harmonised.tsv", row$Matched_File)
    output_path <- file.path(harmonised_root, row$Outcome, relative_file)
    readr::write_tsv(result$data, output_path)
    require_output(output_path, sprintf("Harmonised %s instruments for %s", row$Outcome, row$Metabolite))
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
require_unique_keys(harmonisation_summary, c("Metabolite", "Outcome"), "Reverse-MR harmonisation summary")
require_same_keys(harmonisation_summary, expected_grid, c("Metabolite", "Outcome"), "Reverse-MR harmonisation summary", "Expected reverse-MR grid")
summary_path <- file.path(harmonised_root, "harmonisation_summary.tsv")
readr::write_tsv(harmonisation_summary, summary_path)
require_output(summary_path, "Reverse-MR harmonisation summary")

if (any(harmonisation_summary$Harmonisation_Status != "harmonised")) {
  stop("Reverse-MR harmonisation removed every instrument for at least one expected outcome--metabolite pair.", call. = FALSE)
}
