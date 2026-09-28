source("R/utils.R")

paths <- data_paths("METABOLOME_MR_OUTPUT_DIR")
significance_dir <- file.path(paths[["output_dir"]], "04_significance_filtering")
sensitivity_dir <- file.path(paths[["output_dir"]], "05_sensitivity_analysis")
outcomes <- c("T2DM", "FG", "HBA1C")

read_significant_results <- function(design, outcome) {
  suffix <- if (identical(design, "liberal")) "_liberal" else ""
  path <- file.path(significance_dir, design, sprintf("significant_%s_results%s.tsv", outcome, suffix))
  data <- as.data.frame(readr::read_tsv(path, show_col_types = FALSE))
  # A column that is entirely NA is read as logical; make the filter columns numeric.
  for (column in c("Fixed_IVW_FStat", "Fixed_IVW_HetStat_P", "Egger_Intercept_Pval")) {
    data[[column]] <- as.numeric(data[[column]])
  }
  data
}

# Keep associations with F > 10 and no evidence of heterogeneity (Cochran's Q p <= 0.05)
# or directional pleiotropy (MR-Egger intercept p <= 0.05).
filter_sensitivity <- function(data) {
  passes_f <- !is.na(data$Fixed_IVW_FStat) & data$Fixed_IVW_FStat > 10
  heterogeneous <- !is.na(data$Fixed_IVW_HetStat_P) & data$Fixed_IVW_HetStat_P <= 0.05
  pleiotropic <- !is.na(data$Egger_Intercept_Pval) & data$Egger_Intercept_Pval <= 0.05
  data[passes_f & !heterogeneous & !pleiotropic, , drop = FALSE]
}

for (directory in c("conservative", "liberal", "combined")) {
  dir.create(file.path(sensitivity_dir, directory), recursive = TRUE, showWarnings = FALSE)
}

filtered <- list()
for (outcome in outcomes) {
  conservative <- read_significant_results("conservative", outcome)
  liberal <- read_significant_results("liberal", outcome)
  # The conservative result takes precedence: a liberal result is used only where the conservative one was not significant.
  liberal <- liberal[!(liberal$Metabolite %in% conservative$Metabolite), , drop = FALSE]
  conservative <- filter_sensitivity(conservative)
  liberal <- filter_sensitivity(liberal)
  readr::write_tsv(conservative, file.path(sensitivity_dir, "conservative", sprintf("Filtered_%s_Results.tsv", outcome)))
  readr::write_tsv(liberal, file.path(sensitivity_dir, "liberal", sprintf("Filtered_%s_Results_Liberal.tsv", outcome)))
  filtered[[outcome]] <- rbind(conservative, liberal)
}

# One row per metabolite; each outcome's columns carry the outcome as a suffix.
combined <- data.frame(Metabolite = unique(unlist(lapply(filtered, `[[`, "Metabolite"), use.names = FALSE)))
for (outcome in outcomes) {
  data <- filtered[[outcome]]
  columns <- setdiff(names(data), "Metabolite")
  names(data)[match(columns, names(data))] <- paste0(columns, "_", outcome)
  combined <- merge(combined, data, by = "Metabolite", all.x = TRUE, sort = FALSE)
}
readr::write_tsv(combined, file.path(sensitivity_dir, "combined", "Full_Filtered_Results_Manuscript.tsv"))
