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

as_flag <- function(values, label) {
  if (is.logical(values)) {
    if (anyNA(values)) stop(sprintf("%s has missing boolean values.", label), call. = FALSE)
    return(values)
  }
  mapped <- c("TRUE" = TRUE, "FALSE" = FALSE)[toupper(as.character(values))]
  if (anyNA(mapped)) stop(sprintf("%s has non-boolean values.", label), call. = FALSE)
  unname(as.logical(mapped))
}

require_filename <- function(value, label) {
  if (length(value) != 1L || is.na(value) || !nzchar(value) || basename(value) != value) {
    stop(sprintf("%s must be one plain filename.", label), call. = FALSE)
  }
  value
}

read_region <- function(path, kind, locus) {
  data <- read_tsv(path, sprintf("%s regional file for %s", kind, locus))
  required <- c("SNP", "Source_SNP", "Chromosome", "Position", "EffectAllele", "NonEffectAllele", "EAF", "Beta", "SE", "Pval", "N")
  if (identical(kind, "Outcome")) required <- c(required, "Ncase")
  require_columns(data, required, sprintf("%s regional file for %s", kind, locus))
  require_unique_keys(data, "SNP", sprintf("%s regional file for %s", kind, locus))
  require_numeric(data, c("Position", "EAF", "Beta", "SE", "Pval", "N"), sprintf("%s regional file for %s", kind, locus))
  if (any(data$N <= 0 | data$EAF < 0 | data$EAF > 1 | data$SE <= 0 | data$Pval < 0 | data$Pval > 1)) {
    stop(sprintf("%s regional file has invalid numeric values for %s.", kind, locus), call. = FALSE)
  }
  if (identical(kind, "Outcome") && any(!is.na(data$Ncase) & (!is.finite(data$Ncase) | data$Ncase <= 0))) {
    stop(sprintf("Outcome regional file has invalid case counts for %s.", locus), call. = FALSE)
  }
  data$EffectAllele <- toupper(as.character(data$EffectAllele))
  data$NonEffectAllele <- toupper(as.character(data$NonEffectAllele))
  data
}

keep_snp_alleles <- function(data, locus, kind) {
  keep <- grepl("^[ACGT]$", data$EffectAllele) & grepl("^[ACGT]$", data$NonEffectAllele) & data$EffectAllele != data$NonEffectAllele
  data <- data[keep, , drop = FALSE]
  require_nonempty(data, sprintf("%s SNP alleles for %s", kind, locus))
  data
}

format_pwcoco <- function(data, a1, a2, eaf, beta, se, pval, sample_size, label, ncase = NULL) {
  output <- data.frame(
    SNP = as.character(data$SNP),
    A1 = toupper(as.character(data[[a1]])),
    A2 = toupper(as.character(data[[a2]])),
    A1_freq = as.numeric(data[[eaf]]),
    beta = as.numeric(data[[beta]]),
    se = as.numeric(data[[se]]),
    p = as.numeric(data[[pval]]),
    n = as.numeric(sample_size),
    stringsAsFactors = FALSE
  )
  if (any(!grepl("^[ACGT]$", output$A1) | !grepl("^[ACGT]$", output$A2) | output$A1 == output$A2)) {
    stop(sprintf("%s has non-SNP alleles after harmonisation.", label), call. = FALSE)
  }
  require_numeric(output, c("A1_freq", "beta", "se", "p", "n"), label)
  if (any(output$A1_freq < 0 | output$A1_freq > 1 | output$se <= 0 | output$p < 0 | output$p > 1 | output$n <= 0)) {
    stop(sprintf("%s has invalid PwCoCo values.", label), call. = FALSE)
  }
  flip <- output$A1 > output$A2
  if (any(flip)) {
    old_a1 <- output$A1[flip]
    output$A1[flip] <- output$A2[flip]
    output$A2[flip] <- old_a1
    output$A1_freq[flip] <- 1 - output$A1_freq[flip]
    output$beta[flip] <- -output$beta[flip]
  }
  output$SNP <- paste(output$SNP, output$A1, output$A2, sep = "_")
  require_unique_keys(output, "SNP", label)
  if (!is.null(ncase)) {
    output$ncase <- as.numeric(ncase)
    if (any(!is.finite(output$ncase) | output$ncase <= 0)) {
      stop(sprintf("%s has invalid case counts.", label), call. = FALSE)
    }
  }
  output
}

prepare_locus <- function(row, region_dir, ready_dir) {
  locus <- as.character(row$Locus_ID)
  metabolite_file <- file.path(region_dir, require_filename(row$Metabolite_Region_File, sprintf("Metabolite regional file for %s", locus)))
  outcome_file <- file.path(region_dir, require_filename(row$Outcome_Region_File, sprintf("Outcome regional file for %s", locus)))
  metabolite <- keep_snp_alleles(read_region(metabolite_file, "Metabolite", locus), locus, "Metabolite")
  outcome <- keep_snp_alleles(read_region(outcome_file, "Outcome", locus), locus, "Outcome")

  exposure_data <- metabolite[, c("SNP", "Beta", "SE", "EffectAllele", "NonEffectAllele", "EAF", "Pval"), drop = FALSE]
  names(exposure_data) <- c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "pval.exposure")
  exposure_data$exposure <- row$Candidate_ID
  exposure_data$id.exposure <- row$Candidate_ID
  outcome_data <- outcome[, c("SNP", "Beta", "SE", "EffectAllele", "NonEffectAllele", "EAF", "Pval"), drop = FALSE]
  names(outcome_data) <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "pval.outcome")
  outcome_data$outcome <- row$Outcome
  outcome_data$id.outcome <- row$Outcome
  harmonised <- TwoSampleMR::harmonise_data(exposure_data, outcome_data, action = 2)
  require_columns(
    harmonised,
    c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "pval.exposure", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "pval.outcome", "remove"),
    sprintf("Harmonised PwCoCo data for %s", locus)
  )
  if (anyNA(harmonised$remove)) stop(sprintf("Harmonisation returned missing removal flags for %s.", locus), call. = FALSE)
  harmonised <- harmonised[!harmonised$remove, , drop = FALSE]
  require_nonempty(harmonised, sprintf("Harmonised PwCoCo data for %s", locus))
  require_unique_keys(harmonised, "SNP", sprintf("Harmonised PwCoCo data for %s", locus))

  metabolite_n <- metabolite$N[match(harmonised$SNP, metabolite$SNP)]
  outcome_n <- outcome$N[match(harmonised$SNP, outcome$SNP)]
  if (anyNA(metabolite_n) || anyNA(outcome_n)) stop(sprintf("Sample-size membership is incomplete for %s.", locus), call. = FALSE)
  metabolite_ready <- format_pwcoco(
    harmonised,
    "effect_allele.exposure",
    "other_allele.exposure",
    "eaf.exposure",
    "beta.exposure",
    "se.exposure",
    "pval.exposure",
    metabolite_n,
    sprintf("Metabolite PwCoCo input for %s", locus)
  )
  outcome_ncase <- if (identical(row$Outcome, "T2DM")) outcome$Ncase[match(harmonised$SNP, outcome$SNP)] else NULL
  outcome_ready <- format_pwcoco(
    harmonised,
    "effect_allele.outcome",
    "other_allele.outcome",
    "eaf.outcome",
    "beta.outcome",
    "se.outcome",
    "pval.outcome",
    outcome_n,
    sprintf("Outcome PwCoCo input for %s", locus),
    ncase = outcome_ncase
  )
  stem <- as.character(row$Locus_File_Stem)
  metabolite_ready_file <- file.path(ready_dir, paste0(stem, "_metabolite_pwcoco.tsv"))
  outcome_ready_file <- file.path(ready_dir, paste0(stem, "_outcome_pwcoco.tsv"))
  readr::write_tsv(metabolite_ready, metabolite_ready_file)
  readr::write_tsv(outcome_ready, outcome_ready_file)
  require_output(metabolite_ready_file, sprintf("Metabolite PwCoCo input for %s", locus))
  require_output(outcome_ready_file, sprintf("Outcome PwCoCo input for %s", locus))
  data.frame(
    Candidate_ID = as.character(row$Candidate_ID),
    Association_ID = as.character(row$Association_ID),
    Metabolite = as.character(row$Metabolite),
    Outcome = as.character(row$Outcome),
    Instrument_Design = as.character(row$Instrument_Design),
    Locus_ID = locus,
    Locus_File_Stem = stem,
    SNP = as.character(row$SNP),
    Chromosome = as.character(row$Chromosome),
    Position = as.numeric(row$Position),
    Metabolite_Input_File = file.path("pwcoco", "ready", basename(metabolite_ready_file)),
    Outcome_Input_File = file.path("pwcoco", "ready", basename(outcome_ready_file)),
    stringsAsFactors = FALSE
  )
}

main <- function() {
  for (package in c("readr", "TwoSampleMR")) require_namespace(package)
  paths <- archive_paths(c("METABOLOME_MR_WORK_DIR"))
  work_stage <- file.path(paths[["work_dir"]], "07_colocalisation")
  region_dir <- file.path(work_stage, "regions")
  ready_dir <- file.path(work_stage, "pwcoco", "ready")
  dir.create(ready_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(ready_dir)) stop(sprintf("Cannot create work directory: %s", ready_dir), call. = FALSE)
  eligibility <- read_tsv(file.path(region_dir, "locus_eligibility.tsv"), "Stage 07 locus eligibility")
  required <- c("Candidate_ID", "Association_ID", "Metabolite", "Outcome", "Instrument_Design", "Locus_ID", "Locus_File_Stem", "SNP", "Chromosome", "Position", "Colocalisation_Assessed", "Metabolite_Region_File", "Outcome_Region_File")
  require_columns(eligibility, required, "Stage 07 locus eligibility")
  require_unique_keys(eligibility, "Locus_ID", "Stage 07 locus eligibility")
  eligibility$Colocalisation_Assessed <- as_flag(eligibility$Colocalisation_Assessed, "Stage 07 colocalisation eligibility")
  assessed <- eligibility[eligibility$Colocalisation_Assessed, , drop = FALSE]
  require_nonempty(assessed, "Colocalisation-assessed loci")
  ready <- do.call(rbind, lapply(seq_len(nrow(assessed)), function(index) prepare_locus(assessed[index, , drop = FALSE], region_dir, ready_dir)))
  require_unique_keys(ready, "Locus_ID", "PwCoCo-ready locus manifest")
  manifest_path <- file.path(work_stage, "pwcoco", "pwcoco_ready_manifest.tsv")
  readr::write_tsv(ready, manifest_path)
  require_output(manifest_path, "PwCoCo-ready locus manifest")
}

if (sys.nframe() == 0L) main()
