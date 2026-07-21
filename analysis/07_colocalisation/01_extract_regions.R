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

require_numeric <- function(data, columns, label, positive = FALSE) {
  for (column in columns) {
    values <- suppressWarnings(as.numeric(data[[column]]))
    if (any(!is.finite(values))) {
      stop(sprintf("%s has non-finite numeric values in %s.", label, column), call. = FALSE)
    }
    if (positive && any(values <= 0)) {
      stop(sprintf("%s has non-positive values in %s.", label, column), call. = FALSE)
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

key_vector <- function(data, keys) {
  do.call(paste, c(data[keys], sep = "\r"))
}

require_same_keys <- function(left, right, keys, left_label, right_label) {
  if (!setequal(key_vector(left, keys), key_vector(right, keys))) {
    stop(sprintf("%s and %s have different %s memberships.", left_label, right_label, paste(keys, collapse = ", ")), call. = FALSE)
  }
  invisible(left)
}

stage07_token <- function(value) {
  if (length(value) != 1L || is.na(value) || !nzchar(as.character(value))) {
    stop("Identifiers must be one non-empty value.", call. = FALSE)
  }
  bytes <- as.integer(charToRaw(enc2utf8(as.character(value))))
  paste0("x", paste(sprintf("%02x", bytes), collapse = ""))
}

candidate_id <- function(metabolite) paste0("candidate_", stage07_token(metabolite))

association_id <- function(metabolite, outcome, design) {
  paste0("association_", stage07_token(paste(metabolite, outcome, design, sep = "\r")))
}

locus_id <- function(association, snp, chromosome, position) {
  paste0("locus_", stage07_token(paste(association, snp, chromosome, position, sep = "\r")))
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
    require_numeric(data, "Number_of_IVs", label, positive = TRUE)
    if (any(data$Number_of_IVs != as.integer(data$Number_of_IVs))) {
      stop(sprintf("%s has non-integer instrument counts.", label), call. = FALSE)
    }
  }
  data.frame(
    Metabolite = as.character(data$Metabolite),
    Outcome = outcome,
    Instrument_Design = design,
    Number_of_IVs = as.integer(data$Number_of_IVs),
    stringsAsFactors = FALSE
  )
}

build_association_index <- function(output_root) {
  stage05_root <- file.path(output_root, "05_sensitivity_analysis")
  combined <- read_tsv(
    file.path(stage05_root, "combined", "Full_Filtered_Results_Manuscript.tsv"),
    "Stage 05 post-sensitivity combined result"
  )
  require_columns(combined, c("Metabolite", paste0("Number_of_IVs_", outcomes$code)), "Stage 05 post-sensitivity combined result")
  require_unique_keys(combined, "Metabolite", "Stage 05 post-sensitivity combined result")
  if (nrow(combined) != 54L) {
    stop(sprintf("Stage 05 post-sensitivity combined result must contain 54 metabolites, found %d.", nrow(combined)), call. = FALSE)
  }

  associations <- do.call(rbind, lapply(outcomes$code, function(outcome) {
    conservative <- read_stage05_arm(stage05_arm_path(stage05_root, "conservative", outcome), outcome, "conservative")
    liberal <- read_stage05_arm(stage05_arm_path(stage05_root, "liberal", outcome), outcome, "liberal")
    liberal <- liberal[!(liberal$Metabolite %in% conservative$Metabolite), , drop = FALSE]
    retained <- rbind(conservative, liberal)
    require_unique_keys(retained, c("Metabolite", "Outcome"), sprintf("Stage 05 %s retained associations", outcome))
    retained
  }))
  require_nonempty(associations, "Stage 05 retained associations")
  require_unique_keys(associations, c("Metabolite", "Outcome"), "Stage 05 retained associations")
  if (!setequal(unique(associations$Metabolite), as.character(combined$Metabolite))) {
    stop("Stage 05 retained associations and the 54-member combined set disagree.", call. = FALSE)
  }
  associations$Candidate_ID <- vapply(associations$Metabolite, candidate_id, character(1))
  associations$Association_ID <- mapply(
    association_id,
    associations$Metabolite,
    associations$Outcome,
    associations$Instrument_Design,
    USE.NAMES = FALSE
  )
  require_unique_keys(associations, "Association_ID", "Stage 05 association identifiers")
  list(combined = combined, associations = associations)
}

as_flag <- function(values, label) {
  if (is.logical(values)) {
    if (anyNA(values)) stop(sprintf("%s has missing boolean values.", label), call. = FALSE)
    return(values)
  }
  mapped <- c("TRUE" = TRUE, "FALSE" = FALSE)[toupper(as.character(values))]
  if (anyNA(mapped)) stop(sprintf("%s has non-boolean values.", label), call. = FALSE)
  unname(as.logical(mapped))
}

build_candidate_status <- function(combined, associations, output_root) {
  stage06_root <- file.path(output_root, "06_reverse_mr")
  steiger <- read_tsv(file.path(stage06_root, "steiger", "steiger_results.tsv"), "Stage 06 Steiger results")
  require_columns(steiger, c("Metabolite", "Outcome", "Direction_Flag"), "Stage 06 Steiger results")
  require_unique_keys(steiger, c("Metabolite", "Outcome"), "Stage 06 Steiger results")
  expected_steiger <- data.frame(
    Metabolite = associations$Metabolite,
    Outcome = outcomes$label[match(associations$Outcome, outcomes$code)],
    stringsAsFactors = FALSE
  )
  require_same_keys(steiger, expected_steiger, c("Metabolite", "Outcome"), "Stage 06 Steiger results", "Stage 05 retained associations")
  steiger$Direction_Flag <- as_flag(steiger$Direction_Flag, "Stage 06 Steiger Direction_Flag")
  steiger_excluded <- tapply(!steiger$Direction_Flag, steiger$Metabolite, any)

  post_steiger <- read_tsv(file.path(stage06_root, "steiger", "post_steiger_candidates.tsv"), "Stage 06 post-Steiger candidates")
  require_columns(post_steiger, c("Metabolite", "Retained_after_Steiger"), "Stage 06 post-Steiger candidates")
  require_unique_keys(post_steiger, "Metabolite", "Stage 06 post-Steiger candidates")
  require_same_keys(post_steiger, combined, "Metabolite", "Stage 06 post-Steiger candidates", "Stage 05 combined result")
  post_steiger$Retained_after_Steiger <- as_flag(post_steiger$Retained_after_Steiger, "Stage 06 post-Steiger status")

  status <- data.frame(
    Candidate_ID = vapply(as.character(combined$Metabolite), candidate_id, character(1)),
    Metabolite = as.character(combined$Metabolite),
    Steiger_Excluded = unname(steiger_excluded[as.character(combined$Metabolite)]),
    stringsAsFactors = FALSE
  )
  if (anyNA(status$Steiger_Excluded)) stop("Stage 06 Steiger evidence is incomplete.", call. = FALSE)
  if (!all((!status$Steiger_Excluded) == post_steiger$Retained_after_Steiger[match(status$Metabolite, post_steiger$Metabolite)])) {
    stop("Stage 06 Steiger results disagree with post-Steiger status.", call. = FALSE)
  }
  status$Steiger_Reason <- ifelse(status$Steiger_Excluded, "incorrect_steiger_direction", "steiger_direction_retained")

  reverse <- read_tsv(file.path(stage06_root, "results", "reverse_mr_raw.tsv"), "Stage 06 reverse-MR results")
  reverse_columns <- c("Metabolite", "Outcome", "Number_of_IVs", "Random_IVW_Pval", "Fixed_IVW_Pval", "Egger_Pval")
  require_columns(reverse, reverse_columns, "Stage 06 reverse-MR results")
  require_unique_keys(reverse, c("Metabolite", "Outcome"), "Stage 06 reverse-MR results")
  expected_reverse <- expand.grid(Metabolite = status$Metabolite, Outcome = outcomes$code, stringsAsFactors = FALSE)
  if (nrow(expected_reverse) != 162L) stop("Expected reverse-MR grid must contain 162 metabolite--outcome pairs.", call. = FALSE)
  require_same_keys(reverse, expected_reverse, c("Metabolite", "Outcome"), "Stage 06 reverse-MR results", "Expected reverse-MR grid")
  require_numeric(reverse, c("Number_of_IVs", "Random_IVW_Pval", "Fixed_IVW_Pval", "Egger_Pval"), "Stage 06 reverse-MR results")
  if (any(reverse$Number_of_IVs < 1 | reverse$Number_of_IVs != as.integer(reverse$Number_of_IVs))) {
    stop("Stage 06 reverse-MR results have invalid instrument counts.", call. = FALSE)
  }
  reverse$Selected_IVW_Pval <- ifelse(reverse$Number_of_IVs > 3L, reverse$Random_IVW_Pval, reverse$Fixed_IVW_Pval)
  reverse$Reverse_MR_Flag <- reverse$Selected_IVW_Pval < 0.05 | reverse$Egger_Pval < 0.05
  reverse_excluded <- tapply(reverse$Reverse_MR_Flag, reverse$Metabolite, any)

  status$Reverse_MR_Assessed <- TRUE
  status$Reverse_MR_Excluded <- unname(reverse_excluded[status$Metabolite])
  if (anyNA(status$Reverse_MR_Excluded)) stop("Stage 06 reverse-MR evidence is incomplete.", call. = FALSE)
  status$Reverse_MR_Reason <- ifelse(status$Reverse_MR_Excluded, "reverse_MR_signal", "reverse_MR_not_flagged")
  status$Pre_Coloc_Eligible <- !status$Steiger_Excluded & !status$Reverse_MR_Excluded
  status$Pre_Coloc_Reason <- ifelse(
    status$Steiger_Excluded,
    "incorrect_steiger_direction",
    ifelse(status$Reverse_MR_Excluded, "reverse_MR_signal", "eligible_for_colocalisation")
  )
  require_unique_keys(status, "Candidate_ID", "Stage 07 candidate-status manifest")
  status
}

build_file_lookup <- function(directory, required_columns, key_column, label) {
  if (!dir.exists(directory)) stop(sprintf("%s directory is missing: %s", label, directory), call. = FALSE)
  files <- list.files(directory, pattern = "\\.tsv$", full.names = TRUE)
  files <- files[!grepl("~", basename(files), fixed = TRUE)]
  if (!length(files)) stop(sprintf("%s directory has no TSV files: %s", label, directory), call. = FALSE)
  lookup <- do.call(rbind, lapply(files, function(path) {
    data <- read_tsv(path, basename(path))
    require_columns(data, required_columns, basename(path))
    values <- unique(as.character(data[[key_column]]))
    if (length(values) != 1L || !nzchar(values)) {
      stop(sprintf("%s file must contain one %s: %s", label, key_column, path), call. = FALSE)
    }
    data.frame(Key = values, File = path, stringsAsFactors = FALSE)
  }))
  require_unique_keys(lookup, "Key", label)
  lookup
}

lookup_file <- function(lookup, key, label) {
  matches <- lookup$File[lookup$Key == key]
  if (length(matches) != 1L) {
    stop(sprintf("%s must contain exactly one file for %s.", label, key), call. = FALSE)
  }
  matches[[1]]
}

standardise_metabolite_region <- function(path, metabolite, sample_sizes) {
  data <- read_tsv(path, sprintf("Metabolite regional source for %s", metabolite))
  required <- c("Metabolite", "SNP", "Chromosome", "Position", "EffectAllele", "NonEffectAllele", "EAF", "Beta", "SE", "Pval")
  require_columns(data, required, sprintf("Metabolite regional source for %s", metabolite))
  if (!all(as.character(data$Metabolite) == metabolite)) {
    stop(sprintf("Metabolite regional source contains the wrong metabolite: %s", path), call. = FALSE)
  }
  require_unique_keys(data, "SNP", sprintf("Metabolite regional source for %s", metabolite))
  require_numeric(data, c("Chromosome", "Position", "EAF", "Beta", "SE", "Pval"), sprintf("Metabolite regional source for %s", metabolite))
  index <- match(as.character(data$SNP), sample_sizes$SNP)
  if (anyNA(index)) stop(sprintf("Metabolite sample sizes are incomplete for %s.", metabolite), call. = FALSE)
  result <- data.frame(
    SNP = paste(data$Chromosome, data$Position, sep = "_"),
    Source_SNP = as.character(data$SNP),
    Chromosome = as.character(data$Chromosome),
    Position = as.numeric(data$Position),
    EffectAllele = toupper(as.character(data$EffectAllele)),
    NonEffectAllele = toupper(as.character(data$NonEffectAllele)),
    EAF = as.numeric(data$EAF),
    Beta = as.numeric(data$Beta),
    SE = as.numeric(data$SE),
    Pval = as.numeric(data$Pval),
    N = as.numeric(sample_sizes$Samplesize[index]),
    stringsAsFactors = FALSE
  )
  require_unique_keys(result, "SNP", sprintf("Metabolite regional coordinates for %s", metabolite))
  require_numeric(result, c("Position", "EAF", "Beta", "SE", "Pval", "N"), sprintf("Metabolite regional coordinates for %s", metabolite), positive = FALSE)
  if (any(result$N <= 0)) stop(sprintf("Metabolite sample sizes are non-positive for %s.", metabolite), call. = FALSE)
  result
}

standardise_outcome_region <- function(path, outcome) {
  label <- sprintf("%s regional outcome source", outcome)
  data <- read_tsv(path, label)
  required <- c("Chromosome", "Position", "EffectAllele", "NonEffectAllele", "EAF", "Beta", "SE", "Pval")
  if (identical(outcome, "T2DM")) required <- c(required, "Ncases", "Ncontrols") else required <- c(required, "SampleSize")
  require_columns(data, required, label)
  data <- data[!is.na(data$EAF), , drop = FALSE]
  require_nonempty(data, sprintf("%s rows with non-missing EAF", label))
  require_numeric(data, required[vapply(required, function(column) column %in% c("Chromosome", "Position", "EAF", "Beta", "SE", "Pval", "Ncases", "Ncontrols", "SampleSize"), logical(1))], label)
  n <- if (identical(outcome, "T2DM")) as.numeric(data$Ncases) + as.numeric(data$Ncontrols) else as.numeric(data$SampleSize)
  ncase <- if (identical(outcome, "T2DM")) as.numeric(data$Ncases) else rep(NA_real_, nrow(data))
  if (any(!is.finite(n) | n <= 0)) stop(sprintf("%s has invalid total sample sizes.", label), call. = FALSE)
  result <- data.frame(
    SNP = paste(data$Chromosome, data$Position, sep = "_"),
    Source_SNP = if ("SNP" %in% names(data)) as.character(data$SNP) else paste(data$Chromosome, data$Position, sep = "_"),
    Chromosome = as.character(data$Chromosome),
    Position = as.numeric(data$Position),
    EffectAllele = toupper(as.character(data$EffectAllele)),
    NonEffectAllele = toupper(as.character(data$NonEffectAllele)),
    EAF = as.numeric(data$EAF),
    Beta = as.numeric(data$Beta),
    SE = as.numeric(data$SE),
    Pval = as.numeric(data$Pval),
    N = n,
    Ncase = ncase,
    stringsAsFactors = FALSE
  )
  require_unique_keys(result, "SNP", sprintf("%s regional outcome coordinates", outcome))
  result
}

main <- function() {
  for (package in c("readr")) require_namespace(package)
  paths <- archive_paths(c(
    "METABOLOME_MR_INPUT_DIR",
    "METABOLOME_MR_FULL_METABOLITE_GWAS_DIR",
    "METABOLOME_MR_WORK_DIR",
    "METABOLOME_MR_OUTPUT_DIR"
  ))
  input_dir <- paths[["input_dir"]]
  full_metabolite_gwas_dir <- paths[["full_metabolite_gwas_dir"]]
  work_dir <- paths[["work_dir"]]
  output_dir <- paths[["output_dir"]]
  work_stage <- file.path(work_dir, "07_colocalisation")
  region_dir <- file.path(work_stage, "regions")
  dir.create(region_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(region_dir)) stop(sprintf("Cannot create work directory: %s", region_dir), call. = FALSE)

  index <- build_association_index(output_dir)
  status <- build_candidate_status(index$combined, index$associations, output_dir)
  status_path <- file.path(work_stage, "candidate_status_manifest.tsv")
  readr::write_tsv(status, status_path)
  require_output(status_path, "Stage 07 candidate-status manifest")

  sample_sizes <- read_tsv(file.path(input_dir, "metabolite_sample_sizes.tsv"), "Metabolite sample-size metadata")
  require_columns(sample_sizes, c("SNP", "Samplesize"), "Metabolite sample-size metadata")
  require_unique_keys(sample_sizes, "SNP", "Metabolite sample-size metadata")
  require_numeric(sample_sizes, "Samplesize", "Metabolite sample-size metadata", positive = TRUE)
  sample_sizes$SNP <- as.character(sample_sizes$SNP)

  metabolite_lookup <- build_file_lookup(
    full_metabolite_gwas_dir,
    c("Metabolite", "SNP", "Chromosome", "Position", "EffectAllele", "NonEffectAllele", "EAF", "Beta", "SE", "Pval"),
    "Metabolite",
    "Per-metabolite regional source"
  )
  outcome_paths <- c(
    T2DM = file.path(input_dir, "t2dm_gwas_cleaned.tsv"),
    FG = file.path(work_dir, "fasting_glucose_gwas_cleaned.tsv"),
    HBA1C = file.path(work_dir, "hbA1c_gwas_cleaned.tsv")
  )
  outcome_regions <- lapply(names(outcome_paths), function(outcome) standardise_outcome_region(outcome_paths[[outcome]], outcome))
  names(outcome_regions) <- names(outcome_paths)

  harmonised_lookups <- list()
  for (design in unique(index$associations$Instrument_Design)) {
    for (outcome in outcomes$code) {
      rows <- index$associations[index$associations$Instrument_Design == design & index$associations$Outcome == outcome, , drop = FALSE]
      if (!nrow(rows)) next
      key <- paste(design, outcome, sep = "_")
      harmonised_lookups[[key]] <- build_file_lookup(
        harmonised_directory(output_dir, design, outcome),
        c("Metabolite", "SNP", "proxy", "Chromosome", "Position"),
        "Metabolite",
        sprintf("%s %s harmonised instrument", design, outcome)
      )
    }
  }

  locus_rows <- lapply(seq_len(nrow(index$associations)), function(row_index) {
    association <- index$associations[row_index, , drop = FALSE]
    status_row <- status[status$Candidate_ID == association$Candidate_ID, , drop = FALSE]
    if (nrow(status_row) != 1L) stop("Candidate-status membership is ambiguous.", call. = FALSE)
    lookup_key <- paste(association$Instrument_Design, association$Outcome, sep = "_")
    instrument_path <- lookup_file(
      harmonised_lookups[[lookup_key]],
      association$Metabolite,
      sprintf("%s %s harmonised instrument", association$Instrument_Design, association$Outcome)
    )
    instruments <- read_tsv(instrument_path, sprintf("Harmonised instruments for %s", association$Metabolite))
    require_columns(instruments, c("Metabolite", "SNP", "proxy", "Chromosome", "Position"), sprintf("Harmonised instruments for %s", association$Metabolite))
    require_unique_keys(instruments, "SNP", sprintf("Harmonised instruments for %s", association$Metabolite))
    if (!all(as.character(instruments$Metabolite) == association$Metabolite)) {
      stop(sprintf("Harmonised instrument file contains the wrong metabolite: %s", instrument_path), call. = FALSE)
    }
    require_numeric(instruments, c("Chromosome", "Position"), sprintf("Harmonised instruments for %s", association$Metabolite))
    if (nrow(instruments) != association$Number_of_IVs) {
      stop(sprintf("Harmonised instrument count disagrees with Stage 05 for %s (%s).", association$Metabolite, association$Outcome), call. = FALSE)
    }

    metabolite_path <- lookup_file(metabolite_lookup, association$Metabolite, "Per-metabolite regional source")
    metabolite_region <- NULL
    rows <- lapply(seq_len(nrow(instruments)), function(instrument_index) {
      instrument <- instruments[instrument_index, , drop = FALSE]
      locus <- locus_id(association$Association_ID, as.character(instrument$SNP), instrument$Chromosome, instrument$Position)
      stem <- stage07_token(locus)
      base <- data.frame(
        Candidate_ID = association$Candidate_ID,
        Association_ID = association$Association_ID,
        Metabolite = association$Metabolite,
        Outcome = association$Outcome,
        Instrument_Design = association$Instrument_Design,
        Locus_ID = locus,
        Locus_File_Stem = stem,
        SNP = as.character(instrument$SNP),
        Chromosome = as.character(instrument$Chromosome),
        Position = as.numeric(instrument$Position),
        Proxy_Status = ifelse(is.na(instrument$proxy), "non_proxy", "proxy"),
        Colocalisation_Assessed = FALSE,
        Eligibility_Reason = NA_character_,
        Metabolite_Region_File = NA_character_,
        Outcome_Region_File = NA_character_,
        stringsAsFactors = FALSE
      )
      if (!status_row$Pre_Coloc_Eligible) {
        base$Eligibility_Reason <- status_row$Pre_Coloc_Reason
        return(base)
      }
      if (!is.na(instrument$proxy)) {
        base$Eligibility_Reason <- "proxy_outcome_instrument"
        return(base)
      }
      if (is.null(metabolite_region)) {
        metabolite_region <<- standardise_metabolite_region(metabolite_path, association$Metabolite, sample_sizes)
      }
      lower <- max(0, as.numeric(instrument$Position) - 500000)
      upper <- as.numeric(instrument$Position) + 500000
      outcome_region <- outcome_regions[[association$Outcome]]
      metabolite_window <- metabolite_region[
        metabolite_region$Chromosome == as.character(instrument$Chromosome) & metabolite_region$Position >= lower & metabolite_region$Position <= upper,
        ,
        drop = FALSE
      ]
      outcome_window <- outcome_region[
        outcome_region$Chromosome == as.character(instrument$Chromosome) & outcome_region$Position >= lower & outcome_region$Position <= upper,
        ,
        drop = FALSE
      ]
      require_nonempty(metabolite_window, sprintf("Metabolite regional window for %s", locus))
      require_nonempty(outcome_window, sprintf("Outcome regional window for %s", locus))
      metabolite_file <- file.path(region_dir, paste0(stem, "_metabolite_region.tsv"))
      outcome_file <- file.path(region_dir, paste0(stem, "_outcome_region.tsv"))
      readr::write_tsv(metabolite_window, metabolite_file)
      readr::write_tsv(outcome_window, outcome_file)
      require_output(metabolite_file, sprintf("Metabolite regional window for %s", locus))
      require_output(outcome_file, sprintf("Outcome regional window for %s", locus))
      base$Colocalisation_Assessed <- TRUE
      base$Eligibility_Reason <- "eligible_non_proxy_instrument"
      base$Metabolite_Region_File <- basename(metabolite_file)
      base$Outcome_Region_File <- basename(outcome_file)
      base
    })
    do.call(rbind, rows)
  })
  locus_eligibility <- do.call(rbind, locus_rows)
  require_unique_keys(locus_eligibility, "Locus_ID", "Stage 07 locus eligibility")
  eligibility_path <- file.path(region_dir, "locus_eligibility.tsv")
  readr::write_tsv(locus_eligibility, eligibility_path)
  require_output(eligibility_path, "Stage 07 locus eligibility")
}

if (sys.nframe() == 0L) main()
