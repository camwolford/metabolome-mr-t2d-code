source("config/environment.R")

paths <- archive_paths(c(
  "METABOLOME_MR_FULL_METABOLITE_GWAS_DIR",
  "METABOLOME_MR_WORK_DIR",
  "METABOLOME_MR_OUTPUT_DIR"
))

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

marker_name <- function(chromosome, position, effect_allele, other_allele) {
  paste("chr", chromosome, position, normalise_allele(effect_allele), normalise_allele(other_allele), sep = "_")
}

safe_file_stem <- function(value) {
  gsub("[^[:alnum:]_.-]", "_", value)
}

read_outcome_instruments <- function(path, code) {
  data <- read_tsv(path, sprintf("Clumped %s outcome instruments", code))
  required <- c("Chromosome", "Position", "SNP", "EffectAllele", "NonEffectAllele", "EAF", "Beta", "SE", "Pval")
  require_columns(data, required, sprintf("Clumped %s outcome instruments", code))
  require_unique_keys(data, "SNP", sprintf("Clumped %s outcome instruments", code))
  if (any(!stats::complete.cases(data[required]))) {
    stop(sprintf("Clumped %s outcome instruments contain missing values.", code), call. = FALSE)
  }
  data$MarkerName <- marker_name(data$Chromosome, data$Position, data$EffectAllele, data$NonEffectAllele)
  data$reverseMarkerName <- marker_name(data$Chromosome, data$Position, data$NonEffectAllele, data$EffectAllele)
  require_unique_keys(data, "MarkerName", sprintf("Clumped %s outcome instruments", code))
  data$Fstat <- (data$Beta / data$SE)^2
  if (any(!is.finite(data$Fstat))) stop(sprintf("Clumped %s outcome instruments have invalid F-statistics.", code), call. = FALSE)
  data
}

read_metabolite_gwas <- function(path, metabolite) {
  data <- read_tsv(path, sprintf("Metabolite GWAS for %s", metabolite))
  required <- c("Metabolite", "SNP", "Chromosome", "Position", "EffectAllele", "NonEffectAllele", "EAF", "Beta", "SE", "Pval")
  require_columns(data, required, sprintf("Metabolite GWAS for %s", metabolite))
  if (!all(as.character(data$Metabolite) == metabolite)) {
    stop(sprintf("Metabolite GWAS contains the wrong metabolite: %s", path), call. = FALSE)
  }
  require_unique_keys(data, "SNP", sprintf("Metabolite GWAS for %s", metabolite))
  data$MarkerName <- marker_name(data$Chromosome, data$Position, data$EffectAllele, data$NonEffectAllele)
  require_unique_keys(data, "MarkerName", sprintf("Metabolite GWAS for %s", metabolite))
  data
}

match_instruments <- function(instruments, metabolite_gwas, metabolite, outcome) {
  direct <- match(instruments$MarkerName, metabolite_gwas$MarkerName)
  reversed <- match(instruments$reverseMarkerName, metabolite_gwas$MarkerName)
  if (any(!is.na(direct) & !is.na(reversed))) {
    stop(sprintf("Direct and reversed matches are both present for %s (%s).", metabolite, outcome), call. = FALSE)
  }
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
  require_unique_keys(matched_data, "SNP", sprintf("Matched %s instruments for %s", outcome, metabolite))
  list(matched = matched_data, unmatched = unmatched)
}

require_namespace("readr")

outcomes <- c("T2DM", "FG", "HBA1C")
output_root <- file.path(paths[["output_dir"]], "06_reverse_mr")
steiger_path <- file.path(output_root, "steiger", "post_steiger_candidates.tsv")
candidates <- read_tsv(steiger_path, "Post-Steiger candidates")
require_columns(candidates, c("Metabolite", "Retained_after_Steiger"), "Post-Steiger candidates")
require_unique_keys(candidates, "Metabolite", "Post-Steiger candidates")
if (nrow(candidates) != 54L) {
  stop(sprintf("Post-Steiger candidates must contain 54 metabolites, found %d.", nrow(candidates)), call. = FALSE)
}
if (anyNA(as.logical(candidates$Retained_after_Steiger))) stop("Post-Steiger status contains missing values.", call. = FALSE)
candidate_names <- as.character(candidates$Metabolite)
safe_stems <- safe_file_stem(candidate_names)
if (anyDuplicated(safe_stems)) stop("Candidate metabolite names collide after filename sanitisation.", call. = FALSE)

expected_grid <- expand.grid(
  Metabolite = candidate_names,
  Outcome = outcomes,
  stringsAsFactors = FALSE
)
require_unique_keys(expected_grid, c("Metabolite", "Outcome"), "Expected reverse-MR grid")
if (nrow(expected_grid) != 162L) stop("Expected reverse-MR grid must contain 162 metabolite--outcome pairs.", call. = FALSE)

instrument_root <- file.path(output_root, "outcome_instruments")
outcome_instruments <- stats::setNames(lapply(outcomes, function(code) {
  read_outcome_instruments(file.path(instrument_root, sprintf("%s_clumped.tsv", code)), code)
}), outcomes)

metabolite_root <- paths[["full_metabolite_gwas_dir"]]
if (!dir.exists(metabolite_root)) stop(sprintf("Per-metabolite GWAS directory is missing: %s", metabolite_root), call. = FALSE)
matched_root <- file.path(output_root, "matched")
matching_root <- file.path(output_root, "matching")
for (code in outcomes) dir.create(file.path(matched_root, code), recursive = TRUE, showWarnings = FALSE)
dir.create(matching_root, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(matched_root) || !dir.exists(matching_root)) stop("Cannot create reverse-MR matching output directories.", call. = FALSE)

matching_rows <- list()
unmatched_rows <- list()
row_index <- 0L

for (metabolite in candidate_names) {
  metabolite_path <- file.path(metabolite_root, paste0(metabolite, "_GWAS.tsv"))
  metabolite_gwas <- read_metabolite_gwas(metabolite_path, metabolite)
  for (code in outcomes) {
    result <- match_instruments(outcome_instruments[[code]], metabolite_gwas, metabolite, code)
    matched_count <- if (is.null(result$matched)) 0L else nrow(result$matched)
    unmatched_count <- nrow(result$unmatched)
    relative_file <- NA_character_
    if (matched_count) {
      relative_file <- paste0(safe_file_stem(metabolite), "_", code, "_matched.tsv")
      matched_path <- file.path(matched_root, code, relative_file)
      readr::write_tsv(result$matched, matched_path)
      require_output(matched_path, sprintf("Matched %s instruments for %s", code, metabolite))
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
require_unique_keys(matching_summary, c("Metabolite", "Outcome"), "Reverse-MR matching summary")
require_same_keys(matching_summary, expected_grid, c("Metabolite", "Outcome"), "Reverse-MR matching summary", "Expected reverse-MR grid")
matching_path <- file.path(matching_root, "matching_summary.tsv")
readr::write_tsv(matching_summary, matching_path)
require_output(matching_path, "Reverse-MR matching summary")

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
require_output(unmatched_path, "Unmatched reverse-MR instruments")

if (nrow(unmatched)) {
  stop(sprintf("%d reverse-MR instruments lacked a direct match; inspect %s before continuing.", nrow(unmatched), unmatched_path), call. = FALSE)
}
