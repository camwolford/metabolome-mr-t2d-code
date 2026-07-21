source("config/environment.R")

paths <- archive_paths(c(
  "METABOLOME_MR_INPUT_DIR", "METABOLOME_MR_WORK_DIR", "METABOLOME_MR_OUTPUT_DIR", "METABOLOME_MR_EUR_PANEL_DIR"
))
plink_bin <- require_executable("METABOLOME_MR_PLINK")

require_columns <- function(data, columns, label) {
  absent <- setdiff(columns, names(data))
  if (length(absent)) stop(sprintf("%s is missing: %s", label, paste(absent, collapse = ", ")), call. = FALSE)
  invisible(data)
}

require_nonempty <- function(data, label) {
  if (!nrow(data)) stop(sprintf("%s has zero rows.", label), call. = FALSE)
  invisible(data)
}

marker_name <- function(data, effect, other) {
  paste("chr", data$Chromosome, "_", data$Position, "_", data[[effect]], "_", data[[other]], sep = "")
}

copy_outcome <- function(data, index, outcome, outcome_prefix) {
  data[[paste0(outcome_prefix, "_Chromosome")]][index] <- outcome$Chromosome[1]
  data[[paste0(outcome_prefix, "_Position")]][index] <- outcome$Position[1]
  data[[paste0(outcome_prefix, "_EffectAllele")]][index] <- outcome$EffectAllele[1]
  data[[paste0(outcome_prefix, "_NonEffectAllele")]][index] <- outcome$NonEffectAllele[1]
  data[[paste0(outcome_prefix, "_Beta")]][index] <- outcome$Beta[1]
  data[[paste0(outcome_prefix, "_SE")]][index] <- outcome$SE[1]
  data[[paste0(outcome_prefix, "_EAF")]][index] <- outcome$EAF[1]
  data[[paste0(outcome_prefix, "_Pval")]][index] <- outcome$Pval[1]
  data
}

match_proxy <- function(data, index, outcome_data, outcome_prefix, eur_bfile) {
  proxies <- gwasvcf::get_ld_proxies(
    rsid = data$SNP[index], bfile = eur_bfile, searchspace = NULL,
    tag_kb = 5000, tag_nsnp = 5000, tag_r2 = 0.8, threads = 1, out = tempfile()
  )
  if (!nrow(proxies)) return(list(data = data, found = FALSE))
  require_columns(proxies, c("SNP_B", "CHR_B", "BP_B", "A1", "A2", "B1", "B2", "MAF_B", "R"), "LD proxies")
  proxies <- proxies[proxies$R >= 0.8, , drop = FALSE]
  for (j in seq_len(nrow(proxies))) {
    for (orientation in list(c("B1", "B2"), c("B2", "B1"))) {
      outcome <- outcome_data[
        outcome_data$Chromosome == proxies$CHR_B[j] & outcome_data$Position == proxies$BP_B[j] &
          outcome_data$EffectAllele == proxies[[orientation[1]]][j] &
          outcome_data$NonEffectAllele == proxies[[orientation[2]]][j], , drop = FALSE
      ]
      if (!nrow(outcome)) next
      if (data$EffectAllele[index] == proxies$A1[j] && data$NonEffectAllele[index] == proxies$A2[j]) {
        data$proxy[index] <- proxies$SNP_B[j]
        data$proxy_EAF[index] <- proxies$MAF_B[j]
        data$proxy_R2[index] <- proxies$R[j]
        data$proxy_EffectAllele[index] <- proxies$B1[j]
        data$proxy_NonEffectAllele[index] <- proxies$B2[j]
        return(list(data = copy_outcome(data, index, outcome, outcome_prefix), found = TRUE))
      }
      if (data$EffectAllele[index] == proxies$A2[j] && data$NonEffectAllele[index] == proxies$A1[j]) {
        data$proxy[index] <- proxies$SNP_B[j]
        data$proxy_EAF[index] <- 1 - proxies$MAF_B[j]
        data$proxy_R2[index] <- proxies$R[j]
        data$proxy_EffectAllele[index] <- proxies$B2[j]
        data$proxy_NonEffectAllele[index] <- proxies$B1[j]
        return(list(data = copy_outcome(data, index, outcome, outcome_prefix), found = TRUE))
      }
    }
  }
  list(data = data, found = FALSE)
}

match_outcome <- function(metabolite_data, outcome_data, outcome_prefix, eur_bfile) {
  require_columns(metabolite_data, c("SNP", "Chromosome", "Position", "EffectAllele", "NonEffectAllele", "Beta", "SE", "EAF"), "Clumped instruments")
  metabolite_data$Fstat <- (metabolite_data$Beta / metabolite_data$SE)^2
  metabolite_data$EffectAllele <- ifelse(metabolite_data$EffectAllele == TRUE, "T", metabolite_data$EffectAllele)
  metabolite_data$NonEffectAllele <- ifelse(metabolite_data$NonEffectAllele == TRUE, "T", metabolite_data$NonEffectAllele)
  metabolite_data$tempMarkerName <- marker_name(metabolite_data, "EffectAllele", "NonEffectAllele")
  metabolite_data$reverseMarkerName <- marker_name(metabolite_data, "NonEffectAllele", "EffectAllele")
  for (column in c("proxy", "proxy_EAF", "proxy_R2", "proxy_EffectAllele", "proxy_NonEffectAllele")) metabolite_data[[column]] <- NA
  for (column in c("Chromosome", "Position", "EffectAllele", "NonEffectAllele", "Beta", "SE", "EAF", "Pval")) metabolite_data[[paste0(outcome_prefix, "_", column)]] <- NA
  for (i in seq_len(nrow(metabolite_data))) {
    match_found <- FALSE
    direct <- outcome_data[outcome_data$MarkerName == metabolite_data$tempMarkerName[i], , drop = FALSE]
    if (nrow(direct)) {
      metabolite_data <- copy_outcome(metabolite_data, i, direct, outcome_prefix)
      match_found <- TRUE
    }
    reverse <- outcome_data[outcome_data$MarkerName == metabolite_data$reverseMarkerName[i], , drop = FALSE]
    if (nrow(reverse)) {
      metabolite_data <- copy_outcome(metabolite_data, i, reverse, outcome_prefix)
      match_found <- TRUE
    }
    if (!match_found) {
      proxy_result <- match_proxy(metabolite_data, i, outcome_data, outcome_prefix, eur_bfile)
      metabolite_data <- proxy_result$data
    }
  }
  metabolite_data <- metabolite_data[!is.na(metabolite_data[[paste0(outcome_prefix, "_Beta")]]), , drop = FALSE]
  metabolite_data$tempMarkerName <- NULL
  metabolite_data$reverseMarkerName <- NULL
  metabolite_data
}

input_dir <- file.path(paths[["work_dir"]], "02_instrument_selection", "liberal", "Clumped_IVs")
output_dir <- file.path(paths[["work_dir"]], "02_instrument_selection", "liberal", "Matched_IVs")
eur_bfile <- file.path(paths[["eur_panel_dir"]], "EUR")
if (!dir.exists(input_dir)) stop(sprintf("Clumped instrument directory is missing: %s", input_dir), call. = FALSE)
if (!all(file.exists(paste0(eur_bfile, c(".bed", ".bim", ".fam"))))) stop(sprintf("EUR PLINK panel prefix is incomplete: %s", eur_bfile), call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop(sprintf("Cannot create output directory: %s", output_dir), call. = FALSE)
gwasvcf::set_plink(plink_bin)

outcome_files <- c(
  t2dm = file.path(paths[["input_dir"]], "t2dm_gwas_cleaned.tsv"),
  fg = file.path(paths[["work_dir"]], "fasting_glucose_gwas_cleaned.tsv"),
  hba1c = file.path(paths[["work_dir"]], "hbA1c_gwas_cleaned.tsv")
)
for (path in outcome_files) if (!file.exists(path)) stop(sprintf("Outcome GWAS is missing: %s", path), call. = FALSE)
outcomes <- lapply(outcome_files, function(path) {
  data <- as.data.frame(readr::read_tsv(path, show_col_types = FALSE))
  require_nonempty(data, basename(path))
  require_columns(data, c("Chromosome", "Position", "EffectAllele", "NonEffectAllele", "Beta", "SE", "EAF", "Pval"), basename(path))
  data
})
outcomes$fg <- outcomes$fg[!is.na(outcomes$fg$EAF), , drop = FALSE]
outcomes$hba1c <- outcomes$hba1c[!is.na(outcomes$hba1c$EAF), , drop = FALSE]
for (name in names(outcomes)) outcomes[[name]]$MarkerName <- marker_name(outcomes[[name]], "EffectAllele", "NonEffectAllele")
metabolite_files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)
if (!length(metabolite_files)) stop(sprintf("No clumped instruments found in: %s", input_dir), call. = FALSE)

for (metabolite_file in metabolite_files) {
  metabolite_data <- as.data.frame(readr::read_tsv(metabolite_file, show_col_types = FALSE))
  require_nonempty(metabolite_data, basename(metabolite_file))
  metabolite_name <- sub("_Clumped_IVs_Liberal.*$", "", basename(metabolite_file))
  for (outcome_name in names(outcomes)) {
    matched <- match_outcome(metabolite_data, outcomes[[outcome_name]], outcome_name, eur_bfile)
    if (!nrow(matched)) next
    suffix <- c(t2dm = "T2DM_Matched_IVs_Liberal", fg = "FG_Matched_IVs_Liberal", hba1c = "HBA1C_Matched_IVs_Liberal")[[outcome_name]]
    output_file <- file.path(output_dir, paste0(metabolite_name, suffix, ".tsv"))
    utils::write.table(matched, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    if (!file.exists(output_file)) stop(sprintf("Matched output was not written: %s", output_file), call. = FALSE)
  }
}
