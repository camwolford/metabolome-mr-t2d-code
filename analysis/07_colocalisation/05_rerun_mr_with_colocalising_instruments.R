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

outcomes <- data.frame(
  code = c("T2DM", "FG", "HBA1C"),
  prefix = c("t2dm", "fg", "hba1c"),
  stringsAsFactors = FALSE
)

stage05_arm_path <- function(root, design, outcome) {
  suffix <- if (identical(design, "liberal")) "_Liberal" else ""
  file.path(root, design, sprintf("Filtered_%s_Results%s.tsv", outcome, suffix))
}

harmonised_directory <- function(root, design, outcome) {
  suffix <- if (identical(design, "liberal")) "_Liberal" else ""
  file.path(root, "02_instrument_selection", design, paste0("Harmonised_", outcome, "_IVs", suffix))
}

select_rerun_estimator <- function(number_of_instruments) {
  if (length(number_of_instruments) != 1L || is.na(number_of_instruments) || number_of_instruments < 1L || number_of_instruments != as.integer(number_of_instruments)) {
    stop("The re-run instrument count must be one positive integer.", call. = FALSE)
  }
  if (number_of_instruments == 1L) return("wald_ratio")
  if (number_of_instruments <= 3L) return("fixed_ivw")
  "random_ivw"
}

estimate_wald_ratio <- function(bx, bxse, by, byse) {
  if (length(bx) != 1L || any(!is.finite(c(bx, bxse, by, byse))) || bx == 0 || bxse <= 0 || byse <= 0) {
    stop("Wald-ratio inputs must be finite with non-zero exposure effect and positive standard errors.", call. = FALSE)
  }
  estimate <- by / bx
  standard_error <- sqrt((byse^2 / bx^2) + (by^2 * bxse^2 / bx^4))
  pvalue <- 2 * stats::pnorm(abs(estimate / standard_error), lower.tail = FALSE)
  list(Estimate = estimate, SE = standard_error, Pval = pvalue)
}

run_rerun_mr <- function(data, outcome, metabolite) {
  prefix <- outcomes$prefix[match(outcome, outcomes$code)]
  if (is.na(prefix)) stop(sprintf("Unknown outcome: %s", outcome), call. = FALSE)
  required <- c("SNP", "Beta", "SE", paste0(prefix, "_Beta"), paste0(prefix, "_SE"))
  require_columns(data, required, sprintf("Re-run instruments for %s", metabolite))
  require_unique_keys(data, "SNP", sprintf("Re-run instruments for %s", metabolite))
  require_numeric(data, required[-1], sprintf("Re-run instruments for %s", metabolite))
  if (any(data$SE <= 0 | data[[paste0(prefix, "_SE")]] <= 0)) {
    stop(sprintf("Re-run instruments have non-positive standard errors for %s.", metabolite), call. = FALSE)
  }
  number_of_instruments <- nrow(data)
  if (!number_of_instruments) stop(sprintf("Re-run instrument set is empty for %s (%s).", metabolite, outcome), call. = FALSE)
  estimator <- select_rerun_estimator(number_of_instruments)
  if (identical(estimator, "wald_ratio")) {
    result <- estimate_wald_ratio(data$Beta[[1]], data$SE[[1]], data[[paste0(prefix, "_Beta")]][[1]], data[[paste0(prefix, "_SE")]][[1]])
    return(c(list(Number_of_Remaining_IVs = number_of_instruments, Selected_Estimator = estimator), result))
  }
  input <- MendelianRandomization::mr_input(
    bx = data$Beta,
    bxse = data$SE,
    by = data[[paste0(prefix, "_Beta")]],
    byse = data[[paste0(prefix, "_SE")]]
  )
  model <- if (identical(estimator, "fixed_ivw")) "fixed" else "random"
  result <- MendelianRandomization::mr_ivw(input, model = model)
  values <- c(Estimate = as.numeric(result$Estimate), SE = as.numeric(result$StdError), Pval = as.numeric(result$Pvalue))
  if (any(!is.finite(values)) || values[["SE"]] <= 0 || values[["Pval"]] < 0 || values[["Pval"]] > 1) {
    stop(sprintf("%s re-run MR returned invalid values for %s.", estimator, metabolite), call. = FALSE)
  }
  c(list(Number_of_Remaining_IVs = number_of_instruments, Selected_Estimator = estimator), as.list(values))
}

build_lookup <- function(directory, required_columns, label) {
  if (!dir.exists(directory)) stop(sprintf("%s directory is missing: %s", label, directory), call. = FALSE)
  files <- list.files(directory, pattern = "\\.tsv$", full.names = TRUE)
  files <- files[!grepl("~", basename(files), fixed = TRUE)]
  if (!length(files)) stop(sprintf("%s directory has no TSV files: %s", label, directory), call. = FALSE)
  lookup <- do.call(rbind, lapply(files, function(path) {
    data <- read_tsv(path, basename(path))
    require_columns(data, required_columns, basename(path))
    metabolites <- unique(as.character(data$Metabolite))
    if (length(metabolites) != 1L || !nzchar(metabolites)) {
      stop(sprintf("%s file must contain one metabolite: %s", label, path), call. = FALSE)
    }
    data.frame(Metabolite = metabolites, File = path, stringsAsFactors = FALSE)
  }))
  require_unique_keys(lookup, "Metabolite", label)
  lookup
}

lookup_file <- function(lookup, metabolite, label) {
  matches <- lookup$File[lookup$Metabolite == metabolite]
  if (length(matches) != 1L) stop(sprintf("%s must contain one file for %s.", label, metabolite), call. = FALSE)
  matches[[1]]
}

build_forward_results <- function(output_root) {
  stage05_root <- file.path(output_root, "05_sensitivity_analysis")
  rows <- do.call(rbind, lapply(c("conservative", "liberal"), function(design) {
    do.call(rbind, lapply(outcomes$code, function(outcome) {
      path <- stage05_arm_path(stage05_root, design, outcome)
      data <- read_tsv(path, sprintf("Stage 05 %s %s result", design, outcome), allow_empty = TRUE)
      required <- c("Metabolite", "Number_of_IVs", "Fixed_IVW_Estimate", "Random_IVW_Estimate")
      require_columns(data, required, sprintf("Stage 05 %s %s result", design, outcome))
      if (!nrow(data)) return(data.frame(Metabolite = character(), Outcome = character(), Instrument_Design = character(), Number_of_IVs = integer(), Fixed_IVW_Estimate = numeric(), Random_IVW_Estimate = numeric(), stringsAsFactors = FALSE))
      require_unique_keys(data, "Metabolite", sprintf("Stage 05 %s %s result", design, outcome))
      require_numeric(data, "Number_of_IVs", sprintf("Stage 05 %s %s result", design, outcome))
      data.frame(
        Metabolite = as.character(data$Metabolite),
        Outcome = outcome,
        Instrument_Design = design,
        Number_of_IVs = as.integer(data$Number_of_IVs),
        Fixed_IVW_Estimate = as.numeric(data$Fixed_IVW_Estimate),
        Random_IVW_Estimate = as.numeric(data$Random_IVW_Estimate),
        stringsAsFactors = FALSE
      )
    }))
  }))
  require_unique_keys(rows, c("Metabolite", "Outcome", "Instrument_Design"), "Stage 05 forward-result membership")
  rows
}

selected_original_effect <- function(row) {
  estimator <- select_rerun_estimator(as.integer(row$Number_of_IVs))
  value <- if (identical(estimator, "random_ivw")) row$Random_IVW_Estimate else row$Fixed_IVW_Estimate
  if (length(value) != 1L || !is.finite(value)) {
    stop(sprintf("Selected original %s effect is missing for %s (%s).", estimator, row$Metabolite, row$Outcome), call. = FALSE)
  }
  list(Estimator = estimator, Estimate = as.numeric(value))
}

rerun_reason <- function(threshold_pass, direction_retained, estimate) {
  if (!is.finite(estimate) || estimate == 0 || is.na(direction_retained)) return("rerun_direction_undefined")
  if (threshold_pass && direction_retained) return("retained_after_rerun")
  if (!threshold_pass && !direction_retained) return("rerun_pvalue_above_threshold_and_direction_changed")
  if (!threshold_pass) return("rerun_pvalue_above_threshold")
  "rerun_direction_changed"
}

run_association <- function(association, raw, eligibility, forward_results, instrument_lookups) {
  base <- data.frame(
    Candidate_ID = as.character(association$Candidate_ID),
    Association_ID = as.character(association$Association_ID),
    Metabolite = as.character(association$Metabolite),
    Outcome = as.character(association$Outcome),
    Instrument_Design = as.character(association$Instrument_Design),
    Coloc_Assessed = as.logical(association$Coloc_Assessed),
    Original_Coloc_MR_Pass = as.logical(association$Coloc_MR_Pass),
    Number_Loci_Retained_H4_ge_0.80 = 0L,
    Number_of_Remaining_IVs = NA_integer_,
    Selected_Estimator = NA_character_,
    Re_run_Estimate = NA_real_,
    Re_run_SE = NA_real_,
    Re_run_Pval = NA_real_,
    Original_Selected_Estimator = NA_character_,
    Original_Selected_Estimate = NA_real_,
    Direction_Retained = NA,
    Threshold_Pass = NA,
    Re_run_Pass = NA,
    Re_run_Reason = NA_character_,
    Final_Association_Coloc_MR_Pass = NA,
    stringsAsFactors = FALSE
  )
  if (isTRUE(association$Coloc_MR_Pass)) {
    base$Re_run_Reason <- "not_required_colocalisation_pass"
    base$Final_Association_Coloc_MR_Pass <- TRUE
    return(base)
  }
  if (is.na(association$Coloc_MR_Pass)) {
    base$Re_run_Reason <- as.character(association$Coloc_Audit_Reason)
    return(base)
  }

  high_h4_loci <- unique(raw$Locus_ID[raw$Association_ID == association$Association_ID & raw$H4 >= 0.80])
  high_h4_snps <- unique(eligibility$SNP[eligibility$Locus_ID %in% high_h4_loci])
  if (length(high_h4_loci) != length(high_h4_snps)) {
    stop(sprintf("High-H4 locus mapping is ambiguous for %s.", association$Association_ID), call. = FALSE)
  }
  lookup_key <- paste(association$Instrument_Design, association$Outcome, sep = "_")
  instrument_path <- lookup_file(instrument_lookups[[lookup_key]], association$Metabolite, sprintf("%s %s harmonised instruments", association$Instrument_Design, association$Outcome))
  instruments <- read_tsv(instrument_path, sprintf("Forward instruments for %s", association$Metabolite))
  prefix <- outcomes$prefix[match(association$Outcome, outcomes$code)]
  require_columns(instruments, c("Metabolite", "SNP", "Beta", "SE", "proxy", paste0(prefix, "_Beta"), paste0(prefix, "_SE")), sprintf("Forward instruments for %s", association$Metabolite))
  require_unique_keys(instruments, "SNP", sprintf("Forward instruments for %s", association$Metabolite))
  if (!all(as.character(instruments$Metabolite) == association$Metabolite)) {
    stop(sprintf("Forward instrument file has the wrong metabolite: %s", instrument_path), call. = FALSE)
  }
  remaining <- instruments[as.character(instruments$SNP) %in% high_h4_snps, , drop = FALSE]
  base$Number_Loci_Retained_H4_ge_0.80 <- length(high_h4_loci)
  if (!nrow(remaining)) {
    base$Number_of_Remaining_IVs <- 0L
    base$Re_run_Reason <- "no_colocalising_instruments"
    base$Final_Association_Coloc_MR_Pass <- FALSE
    return(base)
  }
  result <- run_rerun_mr(remaining, association$Outcome, association$Metabolite)
  original <- forward_results[
    forward_results$Metabolite == association$Metabolite & forward_results$Outcome == association$Outcome & forward_results$Instrument_Design == association$Instrument_Design,
    ,
    drop = FALSE
  ]
  if (nrow(original) != 1L) stop(sprintf("Original forward-result membership is ambiguous for %s.", association$Association_ID), call. = FALSE)
  original_effect <- selected_original_effect(original)
  direction_retained <- if (result$Estimate == 0 || original_effect$Estimate == 0) NA else sign(result$Estimate) == sign(original_effect$Estimate)
  threshold_pass <- result$Pval <= 2.826e-5
  pass <- isTRUE(threshold_pass) && isTRUE(direction_retained)
  base$Number_of_Remaining_IVs <- as.integer(result$Number_of_Remaining_IVs)
  base$Selected_Estimator <- result$Selected_Estimator
  base$Re_run_Estimate <- result$Estimate
  base$Re_run_SE <- result$SE
  base$Re_run_Pval <- result$Pval
  base$Original_Selected_Estimator <- original_effect$Estimator
  base$Original_Selected_Estimate <- original_effect$Estimate
  base$Direction_Retained <- direction_retained
  base$Threshold_Pass <- threshold_pass
  base$Re_run_Pass <- pass
  base$Re_run_Reason <- rerun_reason(threshold_pass, direction_retained, result$Estimate)
  base$Final_Association_Coloc_MR_Pass <- pass
  base
}

summarise_candidate <- function(status, reruns) {
  rows <- reruns[reruns$Candidate_ID == status$Candidate_ID, , drop = FALSE]
  if (!nrow(rows)) stop(sprintf("Re-run membership is absent for %s.", status$Candidate_ID), call. = FALSE)
  values <- rows$Final_Association_Coloc_MR_Pass
  assessed <- any(rows$Coloc_Assessed)
  final_pass <- if (any(values %in% TRUE, na.rm = TRUE)) TRUE else if (assessed) FALSE else NA
  reason <- if (!status$Pre_Coloc_Eligible) {
    status$Pre_Coloc_Reason
  } else if (!assessed) {
    "not_assessed_proxy_only"
  } else if (isTRUE(final_pass)) {
    "one_or_more_associations_retained"
  } else {
    "no_assessed_association_retained"
  }
  data.frame(
    Candidate_ID = status$Candidate_ID,
    Metabolite = status$Metabolite,
    Pre_Coloc_Eligible = status$Pre_Coloc_Eligible,
    Coloc_Assessed = assessed,
    Candidate_Coloc_MR_Pass = final_pass,
    Candidate_Rerun_Reason = reason,
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
  for (package in c("readr", "MendelianRandomization")) require_namespace(package)
  paths <- archive_paths(c("METABOLOME_MR_WORK_DIR", "METABOLOME_MR_OUTPUT_DIR"))
  work_stage <- file.path(paths[["work_dir"]], "07_colocalisation")
  output_stage <- file.path(paths[["output_dir"]], "07_colocalisation")
  dir.create(output_stage, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_stage)) stop(sprintf("Cannot create output directory: %s", output_stage), call. = FALSE)

  associations <- read_tsv(file.path(output_stage, "classification", "association_colocalisation.tsv"), "Association colocalisation classification")
  association_columns <- c("Candidate_ID", "Association_ID", "Metabolite", "Outcome", "Instrument_Design", "Coloc_Assessed", "Coloc_MR_Pass", "Coloc_Audit_Reason")
  require_columns(associations, association_columns, "Association colocalisation classification")
  require_unique_keys(associations, "Association_ID", "Association colocalisation classification")
  associations$Coloc_Assessed <- as_flag(associations$Coloc_Assessed, "Association colocalisation assessment")
  associations$Coloc_MR_Pass <- as_flag(associations$Coloc_MR_Pass, "Association colocalisation pass", allow_na = TRUE)

  raw <- read_tsv(file.path(output_stage, "classification", "raw_pwcoco_rows.tsv"), "Raw PwCoCo audit")
  require_columns(raw, c("Association_ID", "Locus_ID", "H4"), "Raw PwCoCo audit")
  require_numeric(raw, "H4", "Raw PwCoCo audit")
  eligibility <- read_tsv(file.path(work_stage, "regions", "locus_eligibility.tsv"), "Stage 07 locus eligibility")
  require_columns(eligibility, c("Locus_ID", "Association_ID", "SNP"), "Stage 07 locus eligibility")
  require_unique_keys(eligibility, "Locus_ID", "Stage 07 locus eligibility")

  forward_results <- build_forward_results(paths[["output_dir"]])
  instrument_lookups <- list()
  for (design in unique(associations$Instrument_Design)) {
    for (outcome in outcomes$code) {
      rows <- associations[associations$Instrument_Design == design & associations$Outcome == outcome, , drop = FALSE]
      if (!nrow(rows)) next
      prefix <- outcomes$prefix[match(outcome, outcomes$code)]
      key <- paste(design, outcome, sep = "_")
      instrument_lookups[[key]] <- build_lookup(
        harmonised_directory(paths[["output_dir"]], design, outcome),
        c("Metabolite", "SNP", "Beta", "SE", "proxy", paste0(prefix, "_Beta"), paste0(prefix, "_SE")),
        sprintf("%s %s harmonised instruments", design, outcome)
      )
    }
  }

  reruns <- do.call(rbind, lapply(seq_len(nrow(associations)), function(index) run_association(associations[index, , drop = FALSE], raw, eligibility, forward_results, instrument_lookups)))
  require_unique_keys(reruns, "Association_ID", "Colocalisation re-run results")

  status <- read_tsv(file.path(work_stage, "candidate_status_manifest.tsv"), "Stage 07 candidate-status manifest")
  require_columns(status, c("Candidate_ID", "Metabolite", "Pre_Coloc_Eligible", "Pre_Coloc_Reason"), "Stage 07 candidate-status manifest")
  require_unique_keys(status, "Candidate_ID", "Stage 07 candidate-status manifest")
  status$Pre_Coloc_Eligible <- as_flag(status$Pre_Coloc_Eligible, "Stage 07 pre-colocalisation eligibility")
  if (!setequal(status$Candidate_ID, reruns$Candidate_ID)) {
    stop("Candidate-status and re-run memberships disagree.", call. = FALSE)
  }
  candidate_reruns <- do.call(rbind, lapply(seq_len(nrow(status)), function(index) summarise_candidate(status[index, , drop = FALSE], reruns)))
  require_unique_keys(candidate_reruns, "Candidate_ID", "Candidate re-run status")

  write_and_check(reruns, file.path(output_stage, "rerun_mr_results.tsv"), "Combined colocalisation re-run results", "Association_ID")
  for (outcome in outcomes$code) {
    outcome_rows <- reruns[reruns$Outcome == outcome, , drop = FALSE]
    filename <- sprintf("%s_rerun_mr_results.tsv", tolower(outcome))
    write_and_check(outcome_rows, file.path(output_stage, filename), sprintf("%s colocalisation re-run results", outcome), "Association_ID")
  }
  write_and_check(candidate_reruns, file.path(output_stage, "candidate_rerun_status.tsv"), "Candidate re-run status", "Candidate_ID")
}

if (sys.nframe() == 0L) main()
