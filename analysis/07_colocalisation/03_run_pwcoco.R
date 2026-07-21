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

read_tsv <- function(path, label) {
  require_file(path, label)
  data <- as.data.frame(readr::read_tsv(path, show_col_types = FALSE))
  require_nonempty(data, label)
  data
}

require_relative_file <- function(work_stage, relative_path, label) {
  if (length(relative_path) != 1L || is.na(relative_path) || !nzchar(relative_path) || grepl("^/", relative_path) || grepl("(^|/)[.][.](/|$)", relative_path)) {
    stop(sprintf("%s must be a relative path below the Stage 07 work directory.", label), call. = FALSE)
  }
  path <- file.path(work_stage, relative_path)
  require_file(path, label)
  resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
  root <- normalizePath(work_stage, winslash = "/", mustWork = TRUE)
  if (!startsWith(resolved, paste0(root, "/"))) {
    stop(sprintf("%s must resolve below the Stage 07 work directory.", label), call. = FALSE)
  }
  resolved
}

pwcoco_arguments <- function(bfile, metabolite_input, outcome_input, output_prefix, chromosome) {
  c(
    "--bfile", shQuote(bfile),
    "--sum_stats1", shQuote(metabolite_input),
    "--sum_stats2", shQuote(outcome_input),
    "--out", shQuote(output_prefix),
    "--chr", as.character(chromosome),
    "--maf", "0.01"
  )
}

main <- function() {
  require_namespace("readr")
  paths <- archive_paths("METABOLOME_MR_WORK_DIR")
  pwcoco <- require_executable("METABOLOME_MR_PWCOCO")
  work_stage <- file.path(paths[["work_dir"]], "07_colocalisation")
  pwcoco_root <- file.path(work_stage, "pwcoco")
  run_dir <- file.path(pwcoco_root, "runs")
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(run_dir)) stop(sprintf("Cannot create work directory: %s", run_dir), call. = FALSE)

  ready <- read_tsv(file.path(pwcoco_root, "pwcoco_ready_manifest.tsv"), "PwCoCo-ready locus manifest")
  required <- c("Candidate_ID", "Association_ID", "Metabolite", "Outcome", "Instrument_Design", "Locus_ID", "Locus_File_Stem", "SNP", "Chromosome", "Position", "Metabolite_Input_File", "Outcome_Input_File")
  require_columns(ready, required, "PwCoCo-ready locus manifest")
  require_unique_keys(ready, "Locus_ID", "PwCoCo-ready locus manifest")
  if (any(!grepl("^x[0-9a-f]+$", as.character(ready$Locus_File_Stem)))) {
    stop("PwCoCo-ready locus manifest has an unsafe output stem.", call. = FALSE)
  }

  runs <- lapply(seq_len(nrow(ready)), function(index) {
    row <- ready[index, , drop = FALSE]
    metabolite_input <- require_relative_file(work_stage, as.character(row$Metabolite_Input_File), sprintf("Metabolite PwCoCo input for %s", row$Locus_ID))
    outcome_input <- require_relative_file(work_stage, as.character(row$Outcome_Input_File), sprintf("Outcome PwCoCo input for %s", row$Locus_ID))
    bfile <- require_bfile_prefix("METABOLOME_MR_UKB_EUR_BFILE", chromosome = as.character(row$Chromosome))
    prefix <- file.path(run_dir, paste0("pwcoco_", row$Locus_File_Stem))
    coloc_path <- paste0(prefix, ".coloc")
    if (file.exists(coloc_path)) {
      stop(sprintf("Refusing to replace an existing PwCoCo output: %s", coloc_path), call. = FALSE)
    }
    status <- system2(
      pwcoco,
      args = pwcoco_arguments(bfile, metabolite_input, outcome_input, prefix, row$Chromosome)
    )
    if (!identical(as.integer(status), 0L)) {
      stop(sprintf("PwCoCo failed for %s with exit status %s.", row$Locus_ID, status), call. = FALSE)
    }
    require_output(coloc_path, sprintf("PwCoCo .coloc output for %s", row$Locus_ID))
    data.frame(
      Candidate_ID = as.character(row$Candidate_ID),
      Association_ID = as.character(row$Association_ID),
      Metabolite = as.character(row$Metabolite),
      Outcome = as.character(row$Outcome),
      Instrument_Design = as.character(row$Instrument_Design),
      Locus_ID = as.character(row$Locus_ID),
      SNP = as.character(row$SNP),
      Chromosome = as.character(row$Chromosome),
      Position = as.numeric(row$Position),
      Relative_Coloc_Path = file.path("pwcoco", "runs", basename(coloc_path)),
      stringsAsFactors = FALSE
    )
  })
  manifest <- do.call(rbind, runs)
  require_unique_keys(manifest, "Locus_ID", "PwCoCo run manifest")
  manifest_path <- file.path(pwcoco_root, "pwcoco_run_manifest.tsv")
  readr::write_tsv(manifest, manifest_path)
  require_output(manifest_path, "PwCoCo run manifest")
}

if (sys.nframe() == 0L) main()
