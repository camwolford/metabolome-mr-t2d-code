source("config/environment.R")

paths <- archive_paths("METABOLOME_MR_OUTPUT_DIR")

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

require_unique_keys <- function(data, columns, label) {
  keys <- do.call(paste, c(data[columns], sep = "\r"))
  if (anyDuplicated(keys)) stop(sprintf("%s has duplicate keys: %s", label, paste(columns, collapse = ", ")), call. = FALSE)
  invisible(data)
}

coerce_numeric_columns <- function(data, columns, label) {
  for (column in columns) {
    if (!is.numeric(data[[column]]) && all(is.na(data[[column]]))) {
      data[[column]] <- as.double(data[[column]])
    }
    if (!is.numeric(data[[column]])) {
      stop(sprintf("%s has a non-numeric %s column.", label, column), call. = FALSE)
    }
  }
  data
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

key_vector <- function(data, columns) {
  do.call(paste, c(data[columns], sep = "\r"))
}

require_same_keys <- function(left, right, columns, left_label, right_label) {
  if (!setequal(key_vector(left, columns), key_vector(right, columns))) {
    stop(sprintf("%s and %s have different %s memberships.", left_label, right_label, paste(columns, collapse = ", ")), call. = FALSE)
  }
  invisible(left)
}

outcomes <- data.frame(
  code = c("T2DM", "FG", "HBA1C"),
  label = c("Type 2 diabetes", "Fasting glucose", "HbA1c"),
  stringsAsFactors = FALSE
)

stage04_summary_path <- function(stage04_root, design, outcome) {
  suffix <- if (identical(design, "liberal")) "_liberal" else ""
  file.path(stage04_root, design, sprintf("significant_%s_results%s.tsv", outcome, suffix))
}

read_stage04_arm <- function(path, outcome, design) {
  label <- sprintf("Stage 04 %s %s result", design, outcome)
  data <- read_tsv(path, label, allow_empty = TRUE)
  required <- c("Metabolite", "Number_of_IVs", "Fixed_IVW_FStat", "Fixed_IVW_HetStat_P", "Egger_Intercept_Pval")
  require_columns(data, required, label)
  data <- coerce_numeric_columns(data, required[-1], label)
  if (nrow(data)) {
    require_unique_keys(data, "Metabolite", label)
    if (any(is.na(data$Metabolite) | !nzchar(data$Metabolite))) {
      stop(sprintf("%s has an empty metabolite name.", label), call. = FALSE)
    }
  }
  data
}

association_index <- function(data, outcome, design) {
  data.frame(
    Metabolite = as.character(data$Metabolite),
    Outcome = rep(outcomes$label[outcomes$code == outcome], nrow(data)),
    Instrument_Design = rep(design, nrow(data)),
    Number_of_IVs = as.integer(data$Number_of_IVs),
    stringsAsFactors = FALSE
  )
}

filter_sensitivity <- function(data, outcome, design) {
  passes_f_statistic <- !is.na(data$Fixed_IVW_FStat) & data$Fixed_IVW_FStat > 10
  fails_heterogeneity <- !is.na(data$Fixed_IVW_HetStat_P) & data$Fixed_IVW_HetStat_P <= 0.05
  fails_pleiotropy <- !is.na(data$Egger_Intercept_Pval) & data$Egger_Intercept_Pval <= 0.05
  retained <- data[passes_f_statistic & !fails_heterogeneity & !fails_pleiotropy, , drop = FALSE]

  reasons <- rbind(
    data.frame(Metabolite = data$Metabolite[!passes_f_statistic], Outcome = rep(outcomes$label[outcomes$code == outcome], sum(!passes_f_statistic)), Instrument_Design = rep(design, sum(!passes_f_statistic)), Reason = rep("F_statistic_not_greater_than_10", sum(!passes_f_statistic)), stringsAsFactors = FALSE),
    data.frame(Metabolite = data$Metabolite[fails_heterogeneity], Outcome = rep(outcomes$label[outcomes$code == outcome], sum(fails_heterogeneity)), Instrument_Design = rep(design, sum(fails_heterogeneity)), Reason = rep("Cochran_Q_pvalue_less_than_or_equal_to_0.05", sum(fails_heterogeneity)), stringsAsFactors = FALSE),
    data.frame(Metabolite = data$Metabolite[fails_pleiotropy], Outcome = rep(outcomes$label[outcomes$code == outcome], sum(fails_pleiotropy)), Instrument_Design = rep(design, sum(fails_pleiotropy)), Reason = rep("MR_Egger_intercept_pvalue_less_than_or_equal_to_0.05", sum(fails_pleiotropy)), stringsAsFactors = FALSE)
  )
  list(retained = retained, reasons = reasons)
}

make_precedence_reasons <- function(liberal, conservative, outcome) {
  duplicated_liberal <- liberal[liberal$Metabolite %in% conservative$Metabolite, , drop = FALSE]
  data.frame(
    Metabolite = duplicated_liberal$Metabolite,
    Outcome = rep(outcomes$label[outcomes$code == outcome], nrow(duplicated_liberal)),
    Instrument_Design = rep("liberal", nrow(duplicated_liberal)),
    Reason = rep("conservative_significant_membership_precedence", nrow(duplicated_liberal)),
    stringsAsFactors = FALSE
  )
}

rename_outcome_columns <- function(data, outcome) {
  columns <- setdiff(names(data), "Metabolite")
  names(data)[match(columns, names(data))] <- paste0(columns, "_", outcome)
  data
}

combine_outcomes <- function(results) {
  all_metabolites <- unique(unlist(lapply(results, function(data) data$Metabolite), use.names = FALSE))
  combined <- data.frame(Metabolite = all_metabolites, stringsAsFactors = FALSE)
  for (outcome in outcomes$code) {
    combined <- merge(combined, rename_outcome_columns(results[[outcome]], outcome), by = "Metabolite", all.x = TRUE, sort = FALSE)
  }
  combined
}

write_and_check <- function(data, path, label, required_columns, keys = NULL) {
  readr::write_tsv(data, path)
  require_output(path, label)
  written <- read_tsv(path, sprintf("Written %s", label), allow_empty = TRUE)
  require_columns(written, required_columns, sprintf("Written %s", label))
  if (!is.null(keys) && nrow(written)) require_unique_keys(written, keys, sprintf("Written %s", label))
  invisible(path)
}

if (!requireNamespace("readr", quietly = TRUE)) {
  stop("Required R package is unavailable: readr", call. = FALSE)
}

output_root <- paths[["output_dir"]]
stage04_root <- file.path(output_root, "04_significance_filtering")
stage05_root <- file.path(output_root, "05_sensitivity_analysis")

stage04_combined_path <- file.path(stage04_root, "combined", "Full_Significant_Results_Manuscript.tsv")
stage04_combined <- read_tsv(stage04_combined_path, "Stage 04 combined significant result")
require_columns(stage04_combined, c("Metabolite", paste0("Number_of_IVs_", outcomes$code)), "Stage 04 combined significant result")
require_unique_keys(stage04_combined, "Metabolite", "Stage 04 combined significant result")
if (nrow(stage04_combined) != 94L) {
  stop(sprintf("Stage 04 combined significant result must contain 94 metabolites, found %d.", nrow(stage04_combined)), call. = FALSE)
}

filtered_by_outcome <- list()
conservative_by_outcome <- list()
liberal_by_outcome <- list()
all_reasons <- list()
all_arm_associations <- list()

for (outcome in outcomes$code) {
  conservative <- read_stage04_arm(stage04_summary_path(stage04_root, "conservative", outcome), outcome, "conservative")
  liberal <- read_stage04_arm(stage04_summary_path(stage04_root, "liberal", outcome), outcome, "liberal")
  conservative_index <- association_index(conservative, outcome, "conservative")
  liberal_index <- association_index(liberal, outcome, "liberal")
  all_arms <- rbind(conservative_index, liberal_index)
  require_unique_keys(all_arms, c("Metabolite", "Outcome", "Instrument_Design"), sprintf("Stage 04 %s associations", outcome))
  all_arm_associations[[outcome]] <- all_arms

  precedence_reasons <- make_precedence_reasons(liberal, conservative, outcome)
  liberal <- liberal[!(liberal$Metabolite %in% conservative$Metabolite), , drop = FALSE]
  conservative_filtered <- filter_sensitivity(conservative, outcome, "conservative")
  liberal_filtered <- filter_sensitivity(liberal, outcome, "liberal")

  conservative_by_outcome[[outcome]] <- conservative_filtered$retained
  liberal_by_outcome[[outcome]] <- liberal_filtered$retained
  filtered_by_outcome[[outcome]] <- rbind(conservative_filtered$retained, liberal_filtered$retained)
  require_unique_keys(filtered_by_outcome[[outcome]], "Metabolite", sprintf("Stage 05 %s retained result", outcome))
  all_reasons[[outcome]] <- rbind(conservative_filtered$reasons, liberal_filtered$reasons, precedence_reasons)
}

expected_outcome_counts <- c(T2DM = 45L, FG = 20L, HBA1C = 23L)
observed_outcome_counts <- vapply(outcomes$code, function(outcome) nrow(filtered_by_outcome[[outcome]]), integer(1))
if (!identical(observed_outcome_counts, expected_outcome_counts)) {
  stop(
    sprintf(
      paste0(
        "Stage 05 retained association counts must be T2DM=45, FG=20, HBA1C=23; ",
        "found T2DM=%d, FG=%d, HBA1C=%d."
      ),
      observed_outcome_counts[["T2DM"]],
      observed_outcome_counts[["FG"]],
      observed_outcome_counts[["HBA1C"]]
    ),
    call. = FALSE
  )
}

all_arm_associations <- do.call(rbind, all_arm_associations)
stage04_pair_membership <- unique(all_arm_associations[c("Metabolite", "Outcome")])
combined_pair_membership <- do.call(rbind, lapply(outcomes$code, function(outcome) {
  data.frame(
    Metabolite = as.character(stage04_combined$Metabolite[!is.na(stage04_combined[[paste0("Number_of_IVs_", outcome)]])]),
    Outcome = outcomes$label[outcomes$code == outcome],
    stringsAsFactors = FALSE
  )
}))
require_same_keys(stage04_pair_membership, combined_pair_membership, c("Metabolite", "Outcome"), "Stage 04 arm associations", "Stage 04 combined significant result")

combined_results <- combine_outcomes(filtered_by_outcome)
require_nonempty(combined_results, "Post-sensitivity combined result")
require_unique_keys(combined_results, "Metabolite", "Post-sensitivity combined result")
if (nrow(combined_results) != 54L) {
  stop(sprintf("Post-sensitivity combined result must contain 54 metabolites, found %d.", nrow(combined_results)), call. = FALSE)
}

exclusions <- do.call(rbind, all_reasons)
require_columns(exclusions, c("Metabolite", "Outcome", "Instrument_Design", "Reason"), "Sensitivity-filter exclusions")

for (directory in c("conservative", "liberal", "combined")) {
  path <- file.path(stage05_root, directory)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop(sprintf("Cannot create output directory: %s", path), call. = FALSE)
}

for (outcome in outcomes$code) {
  conservative_output <- file.path(stage05_root, "conservative", sprintf("Filtered_%s_Results.tsv", outcome))
  liberal_output <- file.path(stage05_root, "liberal", sprintf("Filtered_%s_Results_Liberal.tsv", outcome))
  write_and_check(conservative_by_outcome[[outcome]], conservative_output, sprintf("conservative %s filtered result", outcome), names(conservative_by_outcome[[outcome]]), "Metabolite")
  write_and_check(liberal_by_outcome[[outcome]], liberal_output, sprintf("liberal %s filtered result", outcome), names(liberal_by_outcome[[outcome]]), "Metabolite")
}

combined_path <- file.path(stage05_root, "combined", "Full_Filtered_Results_Manuscript.tsv")
write_and_check(combined_results, combined_path, "post-sensitivity combined result", c("Metabolite", paste0("Number_of_IVs_", outcomes$code)), "Metabolite")

exclusions_path <- file.path(stage05_root, "combined", "sensitivity_filter_exclusions.tsv")
write_and_check(exclusions, exclusions_path, "sensitivity-filter exclusions", c("Metabolite", "Outcome", "Instrument_Design", "Reason"))
