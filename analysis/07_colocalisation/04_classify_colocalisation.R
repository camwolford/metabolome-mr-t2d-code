source("config/environment.R")

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

require_output <- function(path, label) {
  info <- file.info(path)
  if (!file.exists(path) || is.na(info$size) || info$size == 0L) {
    stop(sprintf("%s was not written: %s", label, path), call. = FALSE)
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
  if (any(duplicated(data[keys]))) {
    stop(sprintf("%s has duplicate keys: %s", label, paste(keys, collapse = ", ")), call. = FALSE)
  }
  invisible(data)
}

require_numeric <- function(data, columns, label) {
  for (column in columns) {
    values <- suppressWarnings(as.numeric(data[[column]]))
    if (any(!is.finite(values))) {
      stop(sprintf("%s has non-finite numeric values in %s.", label, column), call. = FALSE)
    }
  }
  invisible(data)
}

read_tsv <- function(path, label, allow_empty = FALSE) {
  require_file(path, label)
  data <- as.data.frame(readr::read_tsv(path, show_col_types = FALSE))
  if (!allow_empty) require_nonempty(data, label)
  data
}

as_flag <- function(values, label, allow_na = FALSE) {
  if (is.logical(values)) {
    if (!allow_na && anyNA(values)) stop(sprintf("%s has missing boolean values.", label), call. = FALSE)
    return(values)
  }
  mapped <- c("TRUE" = TRUE, "FALSE" = FALSE)[toupper(as.character(values))]
  if (!allow_na && anyNA(mapped)) stop(sprintf("%s has non-boolean values.", label), call. = FALSE)
  unname(as.logical(mapped))
}

require_relative_file <- function(work_stage, relative_path, label) {
  if (length(relative_path) != 1L || is.na(relative_path) || !nzchar(relative_path) || grepl("^/", relative_path) || grepl("(^|/)[.][.](/|$)", relative_path)) {
    stop(sprintf("%s must be a relative path below the Stage 07 work directory.", label), call. = FALSE)
  }
  path <- file.path(work_stage, relative_path)
  require_file(path, label)
  resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
  root <- normalizePath(work_stage, winslash = "/", mustWork = TRUE)
  if (!startsWith(resolved, paste0(root, "/"))) {
    stop(sprintf("%s must resolve below the Stage 07 work directory.", label), call. = FALSE)
  }
  resolved
}

raw_audit_reason <- function(nsnps, h4) {
  if (nsnps >= 100 && h4 >= 0.80) return("meets_nsnps_and_H4_rule")
  if (nsnps < 100 && h4 < 0.80) return("nsnps_below_100_and_H4_below_0.80")
  if (nsnps < 100) return("nsnps_below_100")
  "H4_below_0.80"
}

read_raw_coloc <- function(row, work_stage) {
  locus <- as.character(row$Locus_ID)
  path <- require_relative_file(work_stage, as.character(row$Relative_Coloc_Path), sprintf("PwCoCo .coloc output for %s", locus))
  data <- read_tsv(path, sprintf("PwCoCo .coloc output for %s", locus))
  posterior_columns <- paste0("H", 0:4)
  required <- c("Dataset1", "Dataset2", "SNP1", "SNP2", "nsnps", posterior_columns)
  require_columns(data, required, sprintf("PwCoCo .coloc output for %s", locus))
  require_numeric(data, c("nsnps", posterior_columns), sprintf("PwCoCo .coloc output for %s", locus))
  if (any(data$nsnps < 0 | data$nsnps != as.integer(data$nsnps)) || any(as.matrix(data[posterior_columns]) < 0 | as.matrix(data[posterior_columns]) > 1)) {
    stop(sprintf("PwCoCo .coloc output has invalid posterior values for %s.", locus), call. = FALSE)
  }
  raw_pass <- data$nsnps >= 100 & data$H4 >= 0.80
  data.frame(
    Candidate_ID = rep(as.character(row$Candidate_ID), nrow(data)),
    Association_ID = rep(as.character(row$Association_ID), nrow(data)),
    Metabolite = rep(as.character(row$Metabolite), nrow(data)),
    Outcome = rep(as.character(row$Outcome), nrow(data)),
    Instrument_Design = rep(as.character(row$Instrument_Design), nrow(data)),
    Locus_ID = rep(locus, nrow(data)),
    SNP = rep(as.character(row$SNP), nrow(data)),
    Chromosome = rep(as.character(row$Chromosome), nrow(data)),
    Position = rep(as.numeric(row$Position), nrow(data)),
    Raw_Row_ID = paste(locus, seq_len(nrow(data)), sep = "_"),
    Dataset1 = as.character(data$Dataset1),
    Dataset2 = as.character(data$Dataset2),
    SNP1 = as.character(data$SNP1),
    SNP2 = as.character(data$SNP2),
    nsnps = as.numeric(data$nsnps),
    H0 = as.numeric(data$H0),
    H1 = as.numeric(data$H1),
    H2 = as.numeric(data$H2),
    H3 = as.numeric(data$H3),
    H4 = as.numeric(data$H4),
    Raw_Row_Pass = raw_pass,
    Raw_Audit_Reason = vapply(seq_len(nrow(data)), function(index) raw_audit_reason(data$nsnps[index], data$H4[index]), character(1)),
    stringsAsFactors = FALSE
  )
}

summarise_locus <- function(data) {
  max_h4_index <- which.max(data$H4)[[1]]
  pass <- any(data$Raw_Row_Pass)
  reason <- if (pass) {
    "at_least_one_raw_row_meets_global_rule"
  } else if (!any(data$nsnps >= 100)) {
    "all_raw_rows_have_nsnps_below_100"
  } else if (!any(data$H4 >= 0.80)) {
    "all_raw_rows_have_H4_below_0.80"
  } else {
    "no_raw_row_meets_both_global_thresholds"
  }
  data.frame(
    Candidate_ID = data$Candidate_ID[[1]],
    Association_ID = data$Association_ID[[1]],
    Metabolite = data$Metabolite[[1]],
    Outcome = data$Outcome[[1]],
    Instrument_Design = data$Instrument_Design[[1]],
    Locus_ID = data$Locus_ID[[1]],
    SNP = data$SNP[[1]],
    Chromosome = data$Chromosome[[1]],
    Position = data$Position[[1]],
    Number_Raw_Results = nrow(data),
    Max_H4 = data$H4[[max_h4_index]],
    Nsnps_at_Max_H4 = data$nsnps[[max_h4_index]],
    Max_nsnps = max(data$nsnps),
    Locus_Coloc_Pass = pass,
    Locus_Audit_Reason = reason,
    stringsAsFactors = FALSE
  )
}

summarise_association <- function(association, eligibility, locus_summary) {
  rows <- eligibility[eligibility$Association_ID == association$Association_ID, , drop = FALSE]
  assessed <- rows[rows$Colocalisation_Assessed, , drop = FALSE]
  if (!nrow(assessed)) {
    reasons <- unique(as.character(rows$Eligibility_Reason))
    proxy_only <- length(reasons) == 1L && identical(reasons, "proxy_outcome_instrument")
    return(data.frame(
      Candidate_ID = association$Candidate_ID,
      Association_ID = association$Association_ID,
      Metabolite = association$Metabolite,
      Outcome = association$Outcome,
      Instrument_Design = association$Instrument_Design,
      Number_Assessed_Loci = 0L,
      Number_Passing_Loci = 0L,
      Coloc_Assessed = FALSE,
      Coloc_MR_Pass = NA,
      Coloc_Audit_Reason = if (proxy_only) "not_assessed_proxy_only" else paste(reasons, collapse = ";"),
      stringsAsFactors = FALSE
    ))
  }
  loci <- locus_summary[locus_summary$Association_ID == association$Association_ID, , drop = FALSE]
  if (nrow(loci) != nrow(assessed)) {
    stop(sprintf("Locus classification is incomplete for %s.", association$Association_ID), call. = FALSE)
  }
  pass <- all(loci$Locus_Coloc_Pass)
  data.frame(
    Candidate_ID = association$Candidate_ID,
    Association_ID = association$Association_ID,
    Metabolite = association$Metabolite,
    Outcome = association$Outcome,
    Instrument_Design = association$Instrument_Design,
    Number_Assessed_Loci = nrow(loci),
    Number_Passing_Loci = sum(loci$Locus_Coloc_Pass),
    Coloc_Assessed = TRUE,
    Coloc_MR_Pass = pass,
    Coloc_Audit_Reason = if (pass) "all_assessed_loci_pass_global_rule" else "one_or_more_loci_fail_global_rule",
    stringsAsFactors = FALSE
  )
}

summarise_candidate <- function(status, associations) {
  rows <- associations[associations$Candidate_ID == status$Candidate_ID, , drop = FALSE]
  assessed <- rows$Coloc_Assessed
  passes <- rows$Coloc_MR_Pass[!is.na(rows$Coloc_MR_Pass)]
  candidate_assessed <- any(assessed)
  candidate_pass <- if (candidate_assessed) any(passes) else NA
  reason <- if (!status$Pre_Coloc_Eligible) {
    status$Pre_Coloc_Reason
  } else if (!candidate_assessed) {
    "not_assessed_proxy_only"
  } else if (candidate_pass) {
    "one_or_more_associations_pass_global_rule"
  } else {
    "no_assessed_association_passes_global_rule"
  }
  data.frame(
    Candidate_ID = status$Candidate_ID,
    Metabolite = status$Metabolite,
    Pre_Coloc_Eligible = status$Pre_Coloc_Eligible,
    Number_Associations = nrow(rows),
    Number_Assessed_Associations = sum(assessed),
    Number_Assessed_Loci = sum(rows$Number_Assessed_Loci),
    Number_Passing_Loci = sum(rows$Number_Passing_Loci),
    Coloc_Assessed = candidate_assessed,
    Coloc_MR_Pass = candidate_pass,
    Coloc_Audit_Reason = reason,
    stringsAsFactors = FALSE
  )
}

write_and_check <- function(data, path, label, keys) {
  readr::write_tsv(data, path)
  require_output(path, label)
  written <- read_tsv(path, sprintf("Written %s", label), allow_empty = TRUE)
  if (nrow(written)) require_unique_keys(written, keys, sprintf("Written %s", label))
  invisible(path)
}

main <- function() {
  require_namespace("readr")
  paths <- archive_paths(c("METABOLOME_MR_WORK_DIR", "METABOLOME_MR_OUTPUT_DIR"))
  work_stage <- file.path(paths[["work_dir"]], "07_colocalisation")
  output_stage <- file.path(paths[["output_dir"]], "07_colocalisation")
  classification_dir <- file.path(output_stage, "classification")
  dir.create(classification_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(classification_dir)) stop(sprintf("Cannot create output directory: %s", classification_dir), call. = FALSE)

  eligibility <- read_tsv(file.path(work_stage, "regions", "locus_eligibility.tsv"), "Stage 07 locus eligibility")
  eligibility_columns <- c("Candidate_ID", "Association_ID", "Metabolite", "Outcome", "Instrument_Design", "Locus_ID", "SNP", "Chromosome", "Position", "Colocalisation_Assessed", "Eligibility_Reason")
  require_columns(eligibility, eligibility_columns, "Stage 07 locus eligibility")
  require_unique_keys(eligibility, "Locus_ID", "Stage 07 locus eligibility")
  eligibility$Colocalisation_Assessed <- as_flag(eligibility$Colocalisation_Assessed, "Stage 07 colocalisation eligibility")
  assessed <- eligibility[eligibility$Colocalisation_Assessed, , drop = FALSE]
  require_nonempty(assessed, "Colocalisation-assessed loci")

  run_manifest <- read_tsv(file.path(work_stage, "pwcoco", "pwcoco_run_manifest.tsv"), "PwCoCo run manifest")
  run_columns <- c("Candidate_ID", "Association_ID", "Metabolite", "Outcome", "Instrument_Design", "Locus_ID", "SNP", "Chromosome", "Position", "Relative_Coloc_Path")
  require_columns(run_manifest, run_columns, "PwCoCo run manifest")
  require_unique_keys(run_manifest, "Locus_ID", "PwCoCo run manifest")
  if (!setequal(run_manifest$Locus_ID, assessed$Locus_ID)) {
    stop("PwCoCo run manifest and assessed-locus eligibility disagree.", call. = FALSE)
  }
  raw <- do.call(rbind, lapply(seq_len(nrow(run_manifest)), function(index) read_raw_coloc(run_manifest[index, , drop = FALSE], work_stage)))
  require_unique_keys(raw, "Raw_Row_ID", "Raw PwCoCo audit")
  locus_summary <- do.call(rbind, lapply(split(raw, raw$Locus_ID), summarise_locus))
  require_unique_keys(locus_summary, "Locus_ID", "Locus classification")

  associations <- unique(eligibility[c("Candidate_ID", "Association_ID", "Metabolite", "Outcome", "Instrument_Design")])
  require_unique_keys(associations, "Association_ID", "Stage 07 association membership")
  association_summary <- do.call(rbind, lapply(seq_len(nrow(associations)), function(index) summarise_association(associations[index, , drop = FALSE], eligibility, locus_summary)))
  require_unique_keys(association_summary, "Association_ID", "Association colocalisation classification")

  status <- read_tsv(file.path(work_stage, "candidate_status_manifest.tsv"), "Stage 07 candidate-status manifest")
  require_columns(status, c("Candidate_ID", "Metabolite", "Pre_Coloc_Eligible", "Pre_Coloc_Reason"), "Stage 07 candidate-status manifest")
  require_unique_keys(status, "Candidate_ID", "Stage 07 candidate-status manifest")
  status$Pre_Coloc_Eligible <- as_flag(status$Pre_Coloc_Eligible, "Stage 07 pre-colocalisation eligibility")
  if (!setequal(associations$Candidate_ID, status$Candidate_ID)) {
    stop("Stage 07 association membership and candidate-status manifest disagree.", call. = FALSE)
  }
  candidate_summary <- do.call(rbind, lapply(seq_len(nrow(status)), function(index) summarise_candidate(status[index, , drop = FALSE], association_summary)))
  require_unique_keys(candidate_summary, "Candidate_ID", "Candidate colocalisation classification")

  write_and_check(raw, file.path(classification_dir, "raw_pwcoco_rows.tsv"), "Raw PwCoCo audit", "Raw_Row_ID")
  write_and_check(locus_summary, file.path(classification_dir, "locus_classification.tsv"), "Locus classification", "Locus_ID")
  write_and_check(association_summary, file.path(classification_dir, "association_colocalisation.tsv"), "Association colocalisation classification", "Association_ID")
  write_and_check(candidate_summary, file.path(classification_dir, "candidate_colocalisation.tsv"), "Candidate colocalisation classification", "Candidate_ID")
}

if (sys.nframe() == 0L) main()
