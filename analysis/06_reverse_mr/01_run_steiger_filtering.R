source("R/utils.R")

paths <- data_paths(c(
  "METABOLOME_MR_INPUT_DIR",
  "METABOLOME_MR_WORK_DIR",
  "METABOLOME_MR_OUTPUT_DIR"
))

outcomes <- data.frame(
  code = c("T2DM", "FG", "HBA1C"),
  label = c("Type 2 diabetes", "Fasting glucose", "HbA1c"),
  prefix = c("t2dm", "fg", "hba1c"),
  stringsAsFactors = FALSE
)

read_filtered_results <- function(path, outcome, design) {
  data <- as.data.frame(readr::read_tsv(path, show_col_types = FALSE))
  data.frame(
    Metabolite = as.character(data$Metabolite),
    Outcome = outcomes$label[outcomes$code == outcome],
    Instrument_Design = design,
    Number_of_IVs = as.integer(data$Number_of_IVs),
    stringsAsFactors = FALSE
  )
}

build_harmonised_lookup <- function(directory, label) {
  files <- list.files(directory, pattern = "\\.tsv$", full.names = TRUE)
  files <- files[!grepl("~", basename(files), fixed = TRUE)]

  lookup <- do.call(rbind, lapply(files, function(path) {
    data <- as.data.frame(readr::read_tsv(path, show_col_types = FALSE))
    data.frame(Metabolite = as.character(data$Metabolite[[1]]), File = path, stringsAsFactors = FALSE)
  }))
  lookup
}

find_harmonised_file <- function(lookup, metabolite, label) {
  matches <- lookup$File[lookup$Metabolite == metabolite]
  matches[[1]]
}

marker_name <- function(chromosome, position) {
  paste0("chr", chromosome, "_", position)
}

replace_zero_pvalues <- function(values, label) {
  pmax(values, 1e-300)
}

# Stream the outcome GWAS in chunks, keeping only the sample-size columns at instrument positions,
# so the multi-gigabyte files are never held in memory.
read_outcome_gwas <- function(code, positions) {
  if (identical(code, "T2DM")) {
    path <- file.path(paths[["input_dir"]], "t2dm_gwas_cleaned.tsv")
    columns <- readr::cols_only(Chromosome = "c", Position = "d", Ncases = "d", Ncontrols = "d", Neff = "d")
  } else {
    path <- file.path(paths[["work_dir"]], if (identical(code, "FG")) "fasting_glucose_gwas_cleaned.tsv" else "hbA1c_gwas_cleaned.tsv")
    columns <- readr::cols_only(Chromosome = "c", Position = "d", SampleSize = "d")
  }
  keep <- function(chunk, index) chunk[chunk$Position %in% positions, ]
  data <- as.data.frame(readr::read_tsv_chunked(path, readr::DataFrameCallback$new(keep), chunk_size = 1e6, col_types = columns))
  data <- data[!is.na(suppressWarnings(as.numeric(data$Chromosome))), , drop = FALSE]  # instruments are autosomal
  data$Marker <- marker_name(data$Chromosome, data$Position)
  data
}

calculate_steiger <- function(instruments, association, outcome, sample_size_map, outcome_gwas) {
  prefix <- outcome$prefix

  sample_sizes <- sample_size_map$Samplesize[match(instruments$SNP, sample_size_map$SNP)]

  instrument_marker <- marker_name(instruments[[paste0(prefix, "_Chromosome")]], instruments[[paste0(prefix, "_Position")]])
  outcome_rows <- match(instrument_marker, outcome_gwas$Marker)

  pval_exposure <- replace_zero_pvalues(instruments$Pval, sprintf("Exposure p-values for %s", association$Metabolite))
  pval_outcome <- replace_zero_pvalues(instruments[[paste0(prefix, "_Pval")]], sprintf("Outcome p-values for %s", association$Metabolite))
  r_exposure <- TwoSampleMR::get_r_from_pn(pval_exposure, sample_sizes)

  if (identical(outcome$code, "T2DM")) {
    ncases <- outcome_gwas$Ncases[outcome_rows]
    ncontrols <- outcome_gwas$Ncontrols[outcome_rows]
    neff <- outcome_gwas$Neff[outcome_rows]
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
  direction <- TwoSampleMR::directionality_test(direction_data)
  required <- c("snp_r2.exposure", "snp_r2.outcome", "correct_causal_direction", "steiger_pval")

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

output_root <- paths[["output_dir"]]
sensitivity_dir <- file.path(output_root, "05_sensitivity_analysis")
reverse_mr_dir <- file.path(output_root, "06_reverse_mr")
steiger_dir <- file.path(reverse_mr_dir, "steiger")
dir.create(steiger_dir, recursive = TRUE, showWarnings = FALSE)

combined_path <- file.path(sensitivity_dir, "combined", "Full_Filtered_Results_Manuscript.tsv")
combined <- as.data.frame(readr::read_tsv(combined_path, show_col_types = FALSE))

association_index <- do.call(rbind, lapply(outcomes$code, function(code) {
  conservative <- read_filtered_results(filtered_results_path(sensitivity_dir, "conservative", code), code, "conservative")
  liberal <- read_filtered_results(filtered_results_path(sensitivity_dir, "liberal", code), code, "liberal")
  liberal <- liberal[!(liberal$Metabolite %in% conservative$Metabolite), , drop = FALSE]
  retained <- rbind(conservative, liberal)
  retained
}))

sample_size_map <- as.data.frame(readr::read_tsv(file.path(paths[["input_dir"]], "metabolite_sample_sizes.tsv"), show_col_types = FALSE))

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

instruments <- lapply(seq_len(nrow(association_index)), function(index) {
  association <- association_index[index, , drop = FALSE]
  code <- outcomes$code[outcomes$label == association$Outcome]
  lookup_key <- paste(association$Instrument_Design, code, sep = "_")
  as.data.frame(readr::read_tsv(find_harmonised_file(harmonised_lookups[[lookup_key]], association$Metabolite, lookup_key), show_col_types = FALSE))
})
association_codes <- outcomes$code[match(association_index$Outcome, outcomes$label)]
outcome_gwas <- list()
for (code in unique(association_codes)) {
  prefix <- outcomes$prefix[outcomes$code == code]
  positions <- unique(unlist(lapply(instruments[association_codes == code], function(x) x[[paste0(prefix, "_Position")]])))
  outcome_gwas[[code]] <- read_outcome_gwas(code, positions)
}

steiger_results <- do.call(rbind, lapply(seq_len(nrow(association_index)), function(index) {
  outcome <- outcomes[outcomes$code == association_codes[index], , drop = FALSE]
  calculate_steiger(instruments[[index]], association_index[index, , drop = FALSE], outcome, sample_size_map, outcome_gwas[[outcome$code]])
}))

steiger_path <- file.path(steiger_dir, "steiger_results.tsv")
readr::write_tsv(steiger_results, steiger_path)

steiger_exclusions <- steiger_results[!is.na(steiger_results$Direction_Flag) & !steiger_results$Direction_Flag, c("Metabolite", "Outcome", "Instrument_Design"), drop = FALSE]
steiger_exclusions$Reason <- "incorrect_steiger_direction"
exclusions_path <- file.path(steiger_dir, "steiger_exclusions.tsv")
readr::write_tsv(steiger_exclusions, exclusions_path)

excluded_metabolites <- unique(steiger_exclusions$Metabolite)
post_steiger_candidates <- data.frame(
  Metabolite = as.character(combined$Metabolite),
  Retained_after_Steiger = !(combined$Metabolite %in% excluded_metabolites),
  stringsAsFactors = FALSE
)
candidate_path <- file.path(steiger_dir, "post_steiger_candidates.tsv")
readr::write_tsv(post_steiger_candidates, candidate_path)
