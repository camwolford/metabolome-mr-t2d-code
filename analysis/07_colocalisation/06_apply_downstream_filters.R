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

require_expected_candidate_count <- function(values, expected, label) {
  if (!is.logical(values) || anyNA(values)) {
    stop(sprintf("%s must be complete logical values.", label), call. = FALSE)
  }
  observed <- sum(values)
  if (observed != expected) {
    stop(sprintf("%s must contain %d candidates, found %d.", label, expected, observed), call. = FALSE)
  }
  invisible(values)
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
  dir.create(output_stage, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_stage)) stop(sprintf("Cannot create output directory: %s", output_stage), call. = FALSE)

  status <- read_tsv(file.path(work_stage, "candidate_status_manifest.tsv"), "Stage 07 candidate-status manifest")
  status_columns <- c("Candidate_ID", "Metabolite", "Steiger_Excluded", "Steiger_Reason", "Reverse_MR_Assessed", "Reverse_MR_Excluded", "Reverse_MR_Reason", "Pre_Coloc_Eligible", "Pre_Coloc_Reason")
  require_columns(status, status_columns, "Stage 07 candidate-status manifest")
  require_unique_keys(status, "Candidate_ID", "Stage 07 candidate-status manifest")
  if (nrow(status) != 54L) stop(sprintf("Stage 07 candidate-status manifest must contain 54 candidates, found %d.", nrow(status)), call. = FALSE)
  for (column in c("Steiger_Excluded", "Reverse_MR_Assessed", "Reverse_MR_Excluded", "Pre_Coloc_Eligible")) {
    status[[column]] <- as_flag(status[[column]], sprintf("Stage 07 %s", column))
  }
  require_expected_candidate_count(status$Pre_Coloc_Eligible, 26L, "Stage 07 pre-colocalisation eligibility")

  classification <- read_tsv(file.path(output_stage, "classification", "candidate_colocalisation.tsv"), "Candidate colocalisation classification")
  classification_columns <- c("Candidate_ID", "Metabolite", "Coloc_Assessed", "Coloc_MR_Pass", "Coloc_Audit_Reason")
  require_columns(classification, classification_columns, "Candidate colocalisation classification")
  require_unique_keys(classification, "Candidate_ID", "Candidate colocalisation classification")
  classification$Coloc_Assessed <- as_flag(classification$Coloc_Assessed, "Candidate colocalisation assessment")
  classification$Coloc_MR_Pass <- as_flag(classification$Coloc_MR_Pass, "Candidate colocalisation pass", allow_na = TRUE)

  rerun <- read_tsv(file.path(output_stage, "candidate_rerun_status.tsv"), "Candidate re-run status")
  rerun_columns <- c("Candidate_ID", "Metabolite", "Pre_Coloc_Eligible", "Coloc_Assessed", "Candidate_Coloc_MR_Pass", "Candidate_Rerun_Reason")
  require_columns(rerun, rerun_columns, "Candidate re-run status")
  require_unique_keys(rerun, "Candidate_ID", "Candidate re-run status")
  rerun$Pre_Coloc_Eligible <- as_flag(rerun$Pre_Coloc_Eligible, "Candidate re-run pre-colocalisation eligibility")
  rerun$Coloc_Assessed <- as_flag(rerun$Coloc_Assessed, "Candidate re-run assessment")
  rerun$Candidate_Coloc_MR_Pass <- as_flag(rerun$Candidate_Coloc_MR_Pass, "Candidate re-run pass", allow_na = TRUE)

  for (table in list(classification, rerun)) {
    if (!setequal(table$Candidate_ID, status$Candidate_ID)) {
      stop("Candidate-status, classification, and re-run memberships must agree.", call. = FALSE)
    }
  }
  classification <- classification[match(status$Candidate_ID, classification$Candidate_ID), , drop = FALSE]
  rerun <- rerun[match(status$Candidate_ID, rerun$Candidate_ID), , drop = FALSE]
  if (!identical(status$Metabolite, classification$Metabolite) || !identical(status$Metabolite, rerun$Metabolite)) {
    stop("Candidate identifiers map to inconsistent metabolite names.", call. = FALSE)
  }
  if (!all(status$Pre_Coloc_Eligible == rerun$Pre_Coloc_Eligible) || !all(classification$Coloc_Assessed == rerun$Coloc_Assessed)) {
    stop("Candidate hand-off status values disagree.", call. = FALSE)
  }
  require_expected_candidate_count(rerun$Coloc_Assessed, 24L, "Stage 07 colocalisation assessment")

  final_pass <- rerun$Candidate_Coloc_MR_Pass
  final_retained <- status$Pre_Coloc_Eligible & !is.na(final_pass) & final_pass
  require_expected_candidate_count(final_retained, 19L, "Stage 07 final retained set")
  final_reason <- ifelse(
    status$Steiger_Excluded,
    status$Steiger_Reason,
    ifelse(
      status$Reverse_MR_Excluded,
      status$Reverse_MR_Reason,
      ifelse(
        !rerun$Coloc_Assessed,
        "not_assessed_proxy_only",
        ifelse(final_retained, "retained_after_colocalisation", rerun$Candidate_Rerun_Reason)
      )
    )
  )
  downstream_status <- data.frame(
    Candidate_ID = status$Candidate_ID,
    Metabolite = status$Metabolite,
    Steiger_Excluded = status$Steiger_Excluded,
    Steiger_Reason = status$Steiger_Reason,
    Reverse_MR_Assessed = status$Reverse_MR_Assessed,
    Reverse_MR_Excluded = status$Reverse_MR_Excluded,
    Reverse_MR_Reason = status$Reverse_MR_Reason,
    Pre_Coloc_Eligible = status$Pre_Coloc_Eligible,
    Pre_Coloc_Reason = status$Pre_Coloc_Reason,
    Coloc_Assessed = rerun$Coloc_Assessed,
    Initial_Coloc_MR_Pass = classification$Coloc_MR_Pass,
    Coloc_MR_Pass = final_pass,
    Coloc_Reason = rerun$Candidate_Rerun_Reason,
    Final_Retained = final_retained,
    Final_Reason = final_reason,
    stringsAsFactors = FALSE
  )
  require_unique_keys(downstream_status, "Candidate_ID", "Downstream candidate status")
  write_and_check(downstream_status, file.path(output_stage, "candidate_downstream_status.tsv"), "Downstream candidate status", "Candidate_ID")

  retained <- downstream_status[downstream_status$Final_Retained, , drop = FALSE]
  write_and_check(retained, file.path(output_stage, "final_retained_candidates.tsv"), "Final retained candidates", "Candidate_ID")
}

if (sys.nframe() == 0L) main()
