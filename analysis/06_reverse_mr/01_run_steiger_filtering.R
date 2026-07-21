source("config/environment.R")

paths <- archive_paths(c(
  "METABOLOME_MR_INPUT_DIR",
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

require_numeric <- function(data, columns, label) {
  for (column in columns) {
    if (!is.numeric(data[[column]]) || any(!is.finite(data[[column]]))) {
      stop(sprintf("%s has non-finite numeric values in %s.", label, column), call. = FALSE)
    }
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

outcomes <- data.frame(
  code = c("T2DM", "FG", "HBA1C"),
  label = c("Type 2 diabetes", "Fasting glucose", "HbA1c"),
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

read_stage05_arm <- function(path, outcome, design) {
  label <- sprintf("Stage 05 %s %s result", design, outcome)
  data <- read_tsv(path, label, allow_empty = TRUE)
  require_columns(data, c("Metabolite", "Number_of_IVs"), label)
  if (nrow(data)) {
    require_unique_keys(data, "Metabolite", label)
    require_numeric(data, "Number_of_IVs", label)
    if (any(data$Number_of_IVs < 1 | data$Number_of_IVs != as.integer(data$Number_of_IVs))) {
      stop(sprintf("%s has invalid instrument counts.", label), call. = FALSE)
    }
  }
  data.frame(
    Metabolite = as.character(data$Metabolite),
    Outcome = outcomes$label[outcomes$code == outcome],
    Instrument_Design = design,
    Number_of_IVs = as.integer(data$Number_of_IVs),
    stringsAsFactors = FALSE
  )
}

build_harmonised_lookup <- function(directory, label) {
  if (!dir.exists(directory)) stop(sprintf("%s directory is missing: %s", label, directory), call. = FALSE)
  files <- list.files(directory, pattern = "\\.tsv$", full.names = TRUE)
  files <- files[!grepl("~", basename(files), fixed = TRUE)]
  if (!length(files)) stop(sprintf("%s directory has no TSV files: %s", label, directory), call. = FALSE)

  lookup <- do.call(rbind, lapply(files, function(path) {
    data <- read_tsv(path, basename(path))
    require_columns(data, "Metabolite", basename(path))
    if (length(unique(data$Metabolite)) != 1L) {
      stop(sprintf("Harmonised instrument file must contain one metabolite: %s", path), call. = FALSE)
    }
    data.frame(Metabolite = as.character(data$Metabolite[[1]]), File = path, stringsAsFactors = FALSE)
  }))
  require_unique_keys(lookup, "Metabolite", label)
  lookup
}

find_harmonised_file <- function(lookup, metabolite, label) {
  matches <- lookup$File[lookup$Metabolite == metabolite]
  if (length(matches) != 1L) {
    stop(sprintf("%s must contain exactly one harmonised file for metabolite %s.", label, metabolite), call. = FALSE)
  }
  matches[[1]]
}

marker_name <- function(chromosome, position) {
  paste0("chr", chromosome, "_", position)
}

replace_zero_pvalues <- function(values, label) {
  if (any(!is.finite(values) | values < 0 | values > 1)) {
    stop(sprintf("%s contains invalid p-values.", label), call. = FALSE)
  }
  pmax(values, 1e-300)
}

read_outcome_gwas <- function(code) {
  if (identical(code, "T2DM")) {
    path <- file.path(paths[["input_dir"]], "t2dm_gwas_cleaned.tsv")
    required <- c("Chromosome", "Position", "Ncases", "Ncontrols", "Neff")
  } else if (identical(code, "FG")) {
    path <- file.path(paths[["work_dir"]], "fasting_glucose_gwas_cleaned.tsv")
    required <- c("Chromosome", "Position", "SampleSize")
  } else {
    path <- file.path(paths[["work_dir"]], "hbA1c_gwas_cleaned.tsv")
    required <- c("Chromosome", "Position", "SampleSize")
  }
  data <- read_tsv(path, sprintf("%s outcome GWAS", code))
  require_columns(data, required, sprintf("%s outcome GWAS", code))
  require_numeric(data, required, sprintf("%s outcome GWAS", code))
  data
}

calculate_steiger <- function(instruments, association, outcome, sample_size_map, outcome_gwas) {
  prefix <- outcome$prefix
  outcome_columns <- paste0(prefix, c("_Chromosome", "_Position", "_Beta", "_EAF", "_Pval"))
  require_columns(
    instruments,
    c("Metabolite", "SNP", "Pval", "Beta", "SE", outcome_columns),
    sprintf("Harmonised instruments for %s", association$Metabolite)
  )
  require_unique_keys(instruments, "SNP", sprintf("Harmonised instruments for %s", association$Metabolite))
  if (!all(as.character(instruments$Metabolite) == association$Metabolite)) {
    stop(sprintf("Harmonised instrument file contains the wrong metabolite: %s", association$Metabolite), call. = FALSE)
  }
  if (nrow(instruments) != association$Number_of_IVs) {
    stop(sprintf("Harmonised instrument count disagrees with Stage 05 for %s (%s).", association$Metabolite, association$Outcome), call. = FALSE)
  }
  require_numeric(instruments, c("Pval", paste0(prefix, c("_Beta", "_EAF", "_Pval"))), sprintf("Harmonised instruments for %s", association$Metabolite))

  sample_sizes <- sample_size_map$Samplesize[match(instruments$SNP, sample_size_map$SNP)]
  if (anyNA(sample_sizes) || any(!is.finite(sample_sizes) | sample_sizes <= 0)) {
    stop(sprintf("Metabolite sample sizes are incomplete for %s.", association$Metabolite), call. = FALSE)
  }

  source_marker <- marker_name(outcome_gwas$Chromosome, outcome_gwas$Position)
  instrument_marker <- marker_name(instruments[[paste0(prefix, "_Chromosome")]], instruments[[paste0(prefix, "_Position")]])
  outcome_rows <- match(instrument_marker, source_marker)
  if (anyNA(outcome_rows)) {
    stop(sprintf("Outcome sample-size metadata are incomplete for %s (%s).", association$Metabolite, association$Outcome), call. = FALSE)
  }

  pval_exposure <- replace_zero_pvalues(instruments$Pval, sprintf("Exposure p-values for %s", association$Metabolite))
  pval_outcome <- replace_zero_pvalues(instruments[[paste0(prefix, "_Pval")]], sprintf("Outcome p-values for %s", association$Metabolite))
  r_exposure <- TwoSampleMR::get_r_from_pn(pval_exposure, sample_sizes)

  if (identical(outcome$code, "T2DM")) {
    ncases <- outcome_gwas$Ncases[outcome_rows]
    ncontrols <- outcome_gwas$Ncontrols[outcome_rows]
    neff <- outcome_gwas$Neff[outcome_rows]
    if (any(!is.finite(ncases) | !is.finite(ncontrols) | !is.finite(neff) | ncases <= 0 | ncontrols <= 0 | neff <= 0)) {
      stop(sprintf("Type 2 diabetes sample-size metadata are invalid for %s.", association$Metabolite), call. = FALSE)
    }
    r_outcome <- TwoSampleMR::get_r_from_lor(
      instruments[[paste0(prefix, "_Beta")]],
      instruments[[paste0(prefix, "_EAF")]],
      ncases,
      ncontrols,
      prevalence = 0.063
    )
    sample_size_outcome <- neff
  } else {
    sample_size_outcome <- outcome_gwas$SampleSize[outcome_rows]
    if (any(!is.finite(sample_size_outcome) | sample_size_outcome <= 0)) {
      stop(sprintf("%s sample-size metadata are invalid for %s.", association$Outcome, association$Metabolite), call. = FALSE)
    }
    r_outcome <- TwoSampleMR::get_r_from_pn(pval_outcome, sample_size_outcome)
  }

  direction_data <- data.frame(
    SNP = instruments$SNP,
    r.exposure = r_exposure,
    r.outcome = r_outcome,
    id.exposure = "metabolite",
    id.outcome = outcome$code,
    pval.exposure = pval_exposure,
    pval.outcome = pval_outcome,
    samplesize.exposure = sample_sizes,
    samplesize.outcome = sample_size_outcome,
    exposure = "metabolite",
    outcome = outcome$code,
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(direction_data$r.exposure) | !is.finite(direction_data$r.outcome))) {
    stop(sprintf("Steiger variance estimates are invalid for %s (%s).", association$Metabolite, association$Outcome), call. = FALSE)
  }
  direction <- TwoSampleMR::directionality_test(direction_data)
  required <- c("snp_r2.exposure", "snp_r2.outcome", "correct_causal_direction", "steiger_pval")
  if (!all(required %in% names(direction))) {
    stop(sprintf("Steiger calculation returned incomplete output for %s (%s).", association$Metabolite, association$Outcome), call. = FALSE)
  }

  data.frame(
    Metabolite = association$Metabolite,
    Outcome = association$Outcome,
    Instrument_Design = association$Instrument_Design,
    R2_Exposure = as.numeric(direction$snp_r2.exposure[[1]]),
    R2_Outcome = as.numeric(direction$snp_r2.outcome[[1]]),
    Direction_Flag = as.logical(direction$correct_causal_direction[[1]]),
    Steiger_Pval = as.numeric(direction$steiger_pval[[1]]),
    stringsAsFactors = FALSE
  )
}

for (package in c("readr", "TwoSampleMR")) require_namespace(package)

output_root <- paths[["output_dir"]]
stage05_root <- file.path(output_root, "05_sensitivity_analysis")
stage06_root <- file.path(output_root, "06_reverse_mr")
steiger_dir <- file.path(stage06_root, "steiger")
dir.create(steiger_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(steiger_dir)) stop(sprintf("Cannot create output directory: %s", steiger_dir), call. = FALSE)

combined_path <- file.path(stage05_root, "combined", "Full_Filtered_Results_Manuscript.tsv")
combined <- read_tsv(combined_path, "Stage 05 post-sensitivity combined result")
count_columns <- paste0("Number_of_IVs_", outcomes$code)
require_columns(combined, c("Metabolite", count_columns), "Stage 05 post-sensitivity combined result")
require_unique_keys(combined, "Metabolite", "Stage 05 post-sensitivity combined result")
if (nrow(combined) != 54L) {
  stop(sprintf("Stage 05 post-sensitivity combined result must contain 54 metabolites, found %d.", nrow(combined)), call. = FALSE)
}

association_index <- do.call(rbind, lapply(outcomes$code, function(code) {
  conservative <- read_stage05_arm(stage05_arm_path(stage05_root, "conservative", code), code, "conservative")
  liberal <- read_stage05_arm(stage05_arm_path(stage05_root, "liberal", code), code, "liberal")
  liberal <- liberal[!(liberal$Metabolite %in% conservative$Metabolite), , drop = FALSE]
  retained <- rbind(conservative, liberal)
  require_unique_keys(retained, c("Metabolite", "Outcome"), sprintf("Stage 05 %s retained associations", code))
  expected <- data.frame(
    Metabolite = as.character(combined$Metabolite[!is.na(combined[[paste0("Number_of_IVs_", code)]])]),
    Outcome = outcomes$label[outcomes$code == code],
    stringsAsFactors = FALSE
  )
  require_same_keys(retained, expected, c("Metabolite", "Outcome"), sprintf("Stage 05 %s retained associations", code), "Stage 05 post-sensitivity combined result")
  retained
}))
require_nonempty(association_index, "Stage 05 retained associations")
require_unique_keys(association_index, c("Metabolite", "Outcome"), "Stage 05 retained associations")

sample_size_map <- read_tsv(file.path(paths[["input_dir"]], "metabolite_sample_sizes.tsv"), "External metabolite sample-size metadata")
require_columns(sample_size_map, c("SNP", "Samplesize"), "External metabolite sample-size metadata")
require_unique_keys(sample_size_map, "SNP", "External metabolite sample-size metadata")
if (any(is.na(sample_size_map$SNP) | !nzchar(as.character(sample_size_map$SNP)))) {
  stop("External metabolite sample-size metadata has an empty SNP key.", call. = FALSE)
}
require_numeric(sample_size_map, "Samplesize", "External metabolite sample-size metadata")
if (any(sample_size_map$Samplesize <= 0)) {
  stop("External metabolite sample-size metadata has a non-positive Samplesize value.", call. = FALSE)
}
outcome_gwas <- stats::setNames(lapply(outcomes$code, read_outcome_gwas), outcomes$code)

harmonised_lookups <- list()
for (design in unique(association_index$Instrument_Design)) {
  for (code in outcomes$code) {
    associations <- association_index[
      association_index$Instrument_Design == design & association_index$Outcome == outcomes$label[outcomes$code == code],
      ,
      drop = FALSE
    ]
    if (!nrow(associations)) next
    lookup_key <- paste(design, code, sep = "_")
    harmonised_lookups[[lookup_key]] <- build_harmonised_lookup(
      harmonised_directory(output_root, design, code),
      sprintf("%s %s harmonised instrument", design, code)
    )
  }
}

steiger_results <- do.call(rbind, lapply(seq_len(nrow(association_index)), function(index) {
  association <- association_index[index, , drop = FALSE]
  outcome <- outcomes[outcomes$label == association$Outcome, , drop = FALSE]
  lookup_key <- paste(association$Instrument_Design, outcome$code, sep = "_")
  instrument_file <- find_harmonised_file(
    harmonised_lookups[[lookup_key]],
    association$Metabolite,
    sprintf("%s %s harmonised instrument", association$Instrument_Design, outcome$code)
  )
  instruments <- read_tsv(instrument_file, basename(instrument_file))
  calculate_steiger(instruments, association, outcome, sample_size_map, outcome_gwas[[outcome$code]])
}))
require_unique_keys(steiger_results, c("Metabolite", "Outcome"), "Steiger results")
require_same_keys(steiger_results, association_index, c("Metabolite", "Outcome"), "Steiger results", "Stage 05 retained associations")

steiger_path <- file.path(steiger_dir, "steiger_results.tsv")
readr::write_tsv(steiger_results, steiger_path)
require_output(steiger_path, "Steiger results")

steiger_exclusions <- steiger_results[!is.na(steiger_results$Direction_Flag) & !steiger_results$Direction_Flag, c("Metabolite", "Outcome", "Instrument_Design"), drop = FALSE]
steiger_exclusions$Reason <- "incorrect_steiger_direction"
exclusions_path <- file.path(steiger_dir, "steiger_exclusions.tsv")
readr::write_tsv(steiger_exclusions, exclusions_path)
require_output(exclusions_path, "Steiger exclusions")

excluded_metabolites <- unique(steiger_exclusions$Metabolite)
post_steiger_candidates <- data.frame(
  Metabolite = as.character(combined$Metabolite),
  Retained_after_Steiger = !(combined$Metabolite %in% excluded_metabolites),
  stringsAsFactors = FALSE
)
require_unique_keys(post_steiger_candidates, "Metabolite", "Post-Steiger candidates")
candidate_path <- file.path(steiger_dir, "post_steiger_candidates.tsv")
readr::write_tsv(post_steiger_candidates, candidate_path)
require_output(candidate_path, "Post-Steiger candidates")
