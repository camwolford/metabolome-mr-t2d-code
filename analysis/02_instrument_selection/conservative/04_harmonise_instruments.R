source("config/environment.R")

paths <- archive_paths(c(
  "METABOLOME_MR_INPUT_DIR", "METABOLOME_MR_WORK_DIR", "METABOLOME_MR_OUTPUT_DIR"
))

require_columns <- function(data, columns, label) {
  absent <- setdiff(columns, names(data))
  if (length(absent)) stop(sprintf("%s is missing: %s", label, paste(absent, collapse = ", ")), call. = FALSE)
  invisible(data)
}

require_nonempty <- function(data, label) {
  if (!nrow(data)) stop(sprintf("%s has zero rows.", label), call. = FALSE)
  invisible(data)
}

normalise_true <- function(values) ifelse(values == TRUE, "T", values)

harmonise_matched <- function(data, outcome_prefix) {
  outcome_columns <- paste0(outcome_prefix, c("_Beta", "_SE", "_EffectAllele", "_NonEffectAllele", "_EAF"))
  require_columns(
    data,
    c("SNP", "Beta", "SE", "EffectAllele", "NonEffectAllele", "EAF", "proxy", "proxy_EAF", "proxy_EffectAllele", "proxy_NonEffectAllele", outcome_columns),
    "Matched instruments"
  )
  for (column in c("EffectAllele", "NonEffectAllele", paste0(outcome_prefix, c("_EffectAllele", "_NonEffectAllele")), "proxy_EffectAllele", "proxy_NonEffectAllele")) {
    data[[column]] <- normalise_true(data[[column]])
  }
  data <- data[!is.na(data$EAF), , drop = FALSE]
  if (!nrow(data)) return(data)
  exposure <- data[, c("SNP", "Beta", "SE", "EffectAllele", "NonEffectAllele", "EAF"), drop = FALSE]
  outcome <- data[, c("SNP", paste0(outcome_prefix, c("_Beta", "_SE", "_EffectAllele", "_NonEffectAllele", "_EAF"))), drop = FALSE]
  tolerance <- 0.08
  for (i in seq_len(nrow(data))) {
    if (is.na(data$proxy[i])) next
    exposure$EAF[i] <- data$proxy_EAF[i]
    exposure$EffectAllele[i] <- data$proxy_EffectAllele[i]
    exposure$NonEffectAllele[i] <- data$proxy_NonEffectAllele[i]
    lower <- data$EAF[i] - tolerance
    upper <- data$EAF[i] + tolerance
    if (data$proxy_EAF[i] < lower || data$proxy_EAF[i] > upper) {
      inverse <- 1 - data$EAF[i]
      if (data$proxy_EAF[i] > inverse - tolerance && data$proxy_EAF[i] < inverse + tolerance) {
        exposure$EffectAllele[i] <- data$proxy_NonEffectAllele[i]
        exposure$NonEffectAllele[i] <- data$proxy_EffectAllele[i]
        exposure$EAF[i] <- data$EAF[i]
      } else {
        data$SNP[i] <- NA
        exposure$SNP[i] <- NA
        outcome$SNP[i] <- NA
      }
    }
  }
  keep <- !is.na(data$SNP)
  data <- data[keep, , drop = FALSE]
  exposure <- exposure[keep, , drop = FALSE]
  outcome <- outcome[keep, , drop = FALSE]
  if (!nrow(data)) return(data)
  names(exposure) <- c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")
  names(outcome) <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome")
  exposure$exposure <- "exposure"
  exposure$id.exposure <- "exposure"
  outcome$outcome <- "outcome"
  outcome$id.outcome <- "outcome"
  harmonised <- TwoSampleMR::harmonise_data(exposure, outcome, action = 2)
  remove <- harmonised$SNP[harmonised$palindromic & harmonised$ambiguous]
  data <- data[!(data$SNP %in% remove), , drop = FALSE]
  harmonised <- harmonised[!(harmonised$SNP %in% remove), , drop = FALSE]
  for (i in seq_len(nrow(harmonised))) {
    rows <- data$SNP == harmonised$SNP[i]
    if (!any(rows)) next
    if (is.na(data$proxy[which(rows)[1]])) {
      data$Beta[rows] <- harmonised$beta.exposure[i]
      data$SE[rows] <- harmonised$se.exposure[i]
      data$EffectAllele[rows] <- harmonised$effect_allele.exposure[i]
      data$NonEffectAllele[rows] <- harmonised$other_allele.exposure[i]
      data$EAF[rows] <- harmonised$eaf.exposure[i]
    } else {
      data$Beta[rows] <- harmonised$beta.exposure[i]
      data$SE[rows] <- harmonised$se.exposure[i]
      data$proxy_EffectAllele[rows] <- harmonised$effect_allele.exposure[i]
      data$proxy_NonEffectAllele[rows] <- harmonised$other_allele.exposure[i]
      data$proxy_EAF[rows] <- harmonised$eaf.exposure[i]
    }
    data[[paste0(outcome_prefix, "_Beta")]][rows] <- harmonised$beta.outcome[i]
    data[[paste0(outcome_prefix, "_SE")]][rows] <- harmonised$se.outcome[i]
    data[[paste0(outcome_prefix, "_EffectAllele")]][rows] <- harmonised$effect_allele.outcome[i]
    data[[paste0(outcome_prefix, "_NonEffectAllele")]][rows] <- harmonised$other_allele.outcome[i]
    data[[paste0(outcome_prefix, "_EAF")]][rows] <- harmonised$eaf.outcome[i]
  }
  data
}

input_dir <- file.path(paths[["work_dir"]], "02_instrument_selection", "conservative", "Matched_IVs")
output_root <- file.path(paths[["output_dir"]], "02_instrument_selection", "conservative")
if (!dir.exists(input_dir)) stop(sprintf("Matched instrument directory is missing: %s", input_dir), call. = FALSE)
for (outcome in c("T2DM", "FG", "HBA1C")) dir.create(file.path(output_root, paste0("Harmonised_", outcome, "_IVs")), recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_root)) stop(sprintf("Cannot create output directory: %s", output_root), call. = FALSE)

outcomes <- c(T2DM = "t2dm", FG = "fg", HBA1C = "hba1c")
for (outcome in names(outcomes)) {
  metabolite_files <- list.files(input_dir, pattern = paste0(outcome, "_Matched_IVs\\.tsv$"), full.names = TRUE)
  for (metabolite_file in metabolite_files) {
    data <- as.data.frame(readr::read_tsv(metabolite_file, show_col_types = FALSE))
    require_nonempty(data, basename(metabolite_file))
    harmonised <- harmonise_matched(data, outcomes[[outcome]])
    if (!nrow(harmonised)) next
    metabolite_name <- sub(paste0(outcome, "_Matched_IVs.*$"), "", basename(metabolite_file))
    output_file <- file.path(output_root, paste0("Harmonised_", outcome, "_IVs"), paste0(metabolite_name, outcome, "_Harmonised_IVs.tsv"))
    utils::write.table(harmonised, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    if (!file.exists(output_file)) stop(sprintf("Harmonised output was not written: %s", output_file), call. = FALSE)
  }
}
