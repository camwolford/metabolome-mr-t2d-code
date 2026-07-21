archive_root <- function() normalizePath(".", winslash = "/", mustWork = TRUE)

require_directory <- function(variable, root = archive_root()) {
  value <- Sys.getenv(variable, unset = "")
  if (!nzchar(value)) stop(sprintf("Set %s before running this script.", variable), call. = FALSE)
  if (!dir.exists(value)) stop(sprintf("%s must name an existing directory: %s", variable, value), call. = FALSE)
  path <- normalizePath(value, winslash = "/", mustWork = TRUE)
  if (identical(path, root) || startsWith(path, paste0(root, "/"))) {
    stop(sprintf("%s must point outside the repository.", variable), call. = FALSE)
  }
  path
}

require_bfile_prefix <- function(variable, chromosome, root = archive_root()) {
  supported <- "METABOLOME_MR_UKB_EUR_BFILE"
  if (!variable %in% supported) {
    stop(sprintf("Unknown bfile variable: %s", variable), call. = FALSE)
  }
  value <- Sys.getenv(variable, unset = "")
  if (!nzchar(value)) stop(sprintf("Set %s before running this script.", variable), call. = FALSE)
  token_positions <- gregexpr("{chr}", value, fixed = TRUE)[[1]]
  token_count <- if (token_positions[1L] == -1L) 0L else length(token_positions)
  if (token_count != 1L) {
    stop(sprintf("%s must contain exactly one {chr} token.", variable), call. = FALSE)
  }
  if (length(chromosome) != 1L || is.na(chromosome)) {
    stop("chromosome must be a non-empty path-safe chromosome value.", call. = FALSE)
  }
  chromosome <- as.character(chromosome)
  if (!nzchar(chromosome) || !grepl("^[A-Za-z0-9]+$", chromosome)) {
    stop("chromosome must be a non-empty path-safe chromosome value.", call. = FALSE)
  }
  prefix <- sub("{chr}", chromosome, value, fixed = TRUE)
  parent <- dirname(prefix)
  if (!dir.exists(parent)) {
    stop(sprintf("%s must resolve to a prefix in an existing directory: %s", variable, prefix), call. = FALSE)
  }
  path <- file.path(normalizePath(parent, winslash = "/", mustWork = TRUE), basename(prefix))
  if (identical(path, root) || startsWith(path, paste0(root, "/"))) {
    stop(sprintf("%s must point outside the repository.", variable), call. = FALSE)
  }
  components <- paste0(path, c(".bed", ".bim", ".fam"))
  missing <- components[!file.exists(components)]
  if (length(missing)) {
    stop(
      sprintf("%s is missing PLINK binary-reference files: %s", variable, paste(basename(missing), collapse = ", ")),
      call. = FALSE
    )
  }
  resolved_components <- normalizePath(components, winslash = "/", mustWork = TRUE)
  archive_components <- resolved_components == root | startsWith(resolved_components, paste0(root, "/"))
  if (any(archive_components)) {
    stop(
      sprintf("%s component targets must point outside the repository: %s", variable, paste(basename(components[archive_components]), collapse = ", ")),
      call. = FALSE
    )
  }
  directories <- components[dir.exists(components)]
  if (length(directories)) {
    stop(
      sprintf("%s PLINK binary-reference components must be files: %s", variable, paste(basename(directories), collapse = ", ")),
      call. = FALSE
    )
  }
  path
}

archive_paths <- function(required) {
  labels <- c(
    METABOLOME_MR_INPUT_DIR = "input_dir",
    METABOLOME_MR_FULL_METABOLITE_GWAS_DIR = "full_metabolite_gwas_dir",
    METABOLOME_MR_WORK_DIR = "work_dir",
    METABOLOME_MR_OUTPUT_DIR = "output_dir",
    METABOLOME_MR_EUR_PANEL_DIR = "eur_panel_dir"
  )
  unknown <- setdiff(required, names(labels))
  if (length(unknown)) stop(sprintf("Unknown path variable: %s", paste(unknown, collapse = ", ")), call. = FALSE)
  values <- vapply(required, require_directory, character(1), root = archive_root())
  paths <- stats::setNames(values, unname(labels[required]))
  complete_and_work <- c("METABOLOME_MR_FULL_METABOLITE_GWAS_DIR", "METABOLOME_MR_WORK_DIR")
  if (all(complete_and_work %in% required)) {
    sparse_metabolite_gwas_dir <- file.path(paths[["work_dir"]], "Individual_Metabolite_GWAS")
    if (dir.exists(sparse_metabolite_gwas_dir)) {
      sparse_metabolite_gwas_dir <- normalizePath(sparse_metabolite_gwas_dir, winslash = "/", mustWork = TRUE)
      if (identical(paths[["full_metabolite_gwas_dir"]], sparse_metabolite_gwas_dir)) {
        stop(
          "METABOLOME_MR_FULL_METABOLITE_GWAS_DIR must not resolve to the sparse Stage 1 directory under METABOLOME_MR_WORK_DIR.",
          call. = FALSE
        )
      }
    }
  }
  paths
}

require_executable <- function(variable, root = archive_root()) {
  supported <- c("METABOLOME_MR_PLINK", "METABOLOME_MR_PWCOCO")
  if (!variable %in% supported) {
    stop(sprintf("Unknown executable variable: %s", variable), call. = FALSE)
  }
  value <- Sys.getenv(variable, unset = "")
  if (!nzchar(value)) stop(sprintf("Set %s before running this script.", variable), call. = FALSE)
  if (!file.exists(value)) {
    stop(sprintf("%s must name an existing executable: %s", variable, value), call. = FALSE)
  }
  path <- normalizePath(value, winslash = "/", mustWork = TRUE)
  if (dir.exists(path) || file.access(path, mode = 1L) != 0L) {
    stop(sprintf("%s must name an executable: %s", variable, path), call. = FALSE)
  }
  if (identical(path, root) || startsWith(path, paste0(root, "/"))) {
    stop(sprintf("%s must point outside the repository.", variable), call. = FALSE)
  }
  path
}
