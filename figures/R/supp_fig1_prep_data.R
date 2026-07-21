# Supplementary Figure 1 locus-table preparation from external regional inputs.

suppressPackageStartupMessages({
  library(readr)
})

# External paths
resolve_eur_bfile <- function(panel_dir) {
  prefix <- file.path(panel_dir, "EUR")
  components <- paste0(prefix, c(".bed", ".bim", ".fam"))
  invalid <- components[!file.exists(components) | dir.exists(components)]
  if (length(invalid)) {
    stop(sprintf("EUR PLINK panel prefix is incomplete: %s", prefix), call. = FALSE)
  }
  prefix
}

select_stage07_locus <- function(cfg, eligibility, ready) {
  positions <- suppressWarnings(as.numeric(eligibility$Position))
  matched <- eligibility[
    which(
      as.character(eligibility$Metabolite) == cfg$metab &
        as.character(eligibility$Outcome) == cfg$outcome_code &
        as.character(eligibility$SNP) == cfg$rsid &
        as.character(eligibility$Chromosome) == as.character(cfg$chr) &
        !is.na(positions) & positions == cfg$idx_pos
    ),
    , drop = FALSE
  ]
  if (nrow(matched) != 1L) {
    stop(sprintf("%s: expected exactly one Stage 07 locus, found %d.", cfg$id, nrow(matched)), call. = FALSE)
  }
  if (!identical(toupper(as.character(matched$Colocalisation_Assessed[[1]])), "TRUE")) {
    stop(sprintf("%s: selected Stage 07 locus is not colocalisation-assessed.", cfg$id), call. = FALSE)
  }
  ready_row <- ready[as.character(ready$Locus_ID) == as.character(matched$Locus_ID[[1]]), , drop = FALSE]
  if (nrow(ready_row) != 1L) {
    stop(sprintf("%s: expected exactly one PwCoCo-ready manifest row, found %d.", cfg$id, nrow(ready_row)), call. = FALSE)
  }
  shared <- c("Candidate_ID", "Association_ID", "Metabolite", "Outcome", "Instrument_Design", "Locus_ID", "Locus_File_Stem", "SNP", "Chromosome")
  matched_position <- suppressWarnings(as.numeric(matched$Position[[1]]))
  ready_position <- suppressWarnings(as.numeric(ready_row$Position[[1]]))
  if (any(vapply(shared, function(column) {
    !identical(as.character(matched[[column]][[1]]), as.character(ready_row[[column]][[1]]))
  }, logical(1))) || !is.finite(matched_position) || !is.finite(ready_position) || matched_position != ready_position) {
    stop(sprintf("%s: Stage 07 eligibility and PwCoCo-ready metadata disagree.", cfg$id), call. = FALSE)
  }
  list(eligibility = matched, ready = ready_row)
}

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript.", call. = FALSE)
ARCHIVE_ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), "..", ".."), mustWork = TRUE)
source(file.path(ARCHIVE_ROOT, "config", "environment.R"))
paths <- archive_paths(c(
  "METABOLOME_MR_INPUT_DIR",
  "METABOLOME_MR_WORK_DIR",
  "METABOLOME_MR_OUTPUT_DIR",
  "METABOLOME_MR_EUR_PANEL_DIR"
))
COL <- file.path(paths[["work_dir"]], "07_colocalisation")
REGIONS <- file.path(COL, "regions")
ELIGIBILITY_FILE <- file.path(REGIONS, "locus_eligibility.tsv")
READY_MANIFEST_FILE <- file.path(COL, "pwcoco", "pwcoco_ready_manifest.tsv")
WORK <- file.path(paths[["work_dir"]], "figures", "supp_fig1")
OUT <- WORK
EUR <- resolve_eur_bfile(paths[["eur_panel_dir"]])
PLINK <- require_executable("METABOLOME_MR_PLINK")
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Required work inputs
stopifnot_loud <- function(cond, msg) if (!isTRUE(cond)) stop(msg, call. = FALSE)
stopifnot_loud(dir.exists(COL),
  sprintf("Colocalisation work inputs not found at: %s", COL))
stopifnot_loud(file.exists(ELIGIBILITY_FILE) && !dir.exists(ELIGIBILITY_FILE),
  sprintf("Stage 07 locus eligibility is missing: %s", ELIGIBILITY_FILE))
stopifnot_loud(file.exists(READY_MANIFEST_FILE) && !dir.exists(READY_MANIFEST_FILE),
  sprintf("Stage 07 PwCoCo-ready manifest is missing: %s", READY_MANIFEST_FILE))

eligibility <- readr::read_tsv(ELIGIBILITY_FILE, col_types = cols(), progress = FALSE)
ready <- readr::read_tsv(READY_MANIFEST_FILE, col_types = cols(), progress = FALSE)
eligibility_columns <- c(
  "Candidate_ID", "Association_ID", "Metabolite", "Outcome", "Instrument_Design", "Locus_ID", "Locus_File_Stem",
  "SNP", "Chromosome", "Position", "Colocalisation_Assessed", "Metabolite_Region_File", "Outcome_Region_File"
)
ready_columns <- c(
  "Candidate_ID", "Association_ID", "Metabolite", "Outcome", "Instrument_Design", "Locus_ID", "Locus_File_Stem",
  "SNP", "Chromosome", "Position", "Metabolite_Input_File", "Outcome_Input_File"
)
stopifnot_loud(nrow(eligibility) > 0L && nrow(ready) > 0L,
  "Stage 07 locus manifests must not be empty.")
stopifnot_loud(all(eligibility_columns %in% names(eligibility)),
  "Stage 07 locus eligibility has an invalid schema.")
stopifnot_loud(all(ready_columns %in% names(ready)),
  "Stage 07 PwCoCo-ready manifest has an invalid schema.")
stopifnot_loud(!anyDuplicated(eligibility$Locus_ID) && !anyDuplicated(ready$Locus_ID),
  "Stage 07 locus manifests contain duplicate Locus_ID values.")

# Six plotted loci
loci <- list(
  list(id="1_sphinganine_t2d",       metab="sphinganine",                outcome_code="T2DM",  outcome_label="Type 2 Diabetes",       gene="ABO",     rsid="rs676457",   chr=9,  idx_pos=136146227),
  list(id="2_palmitoleoylGPC_t2d",   metab="1-palmitoleoyl-GPC (16_1)*",  outcome_code="T2DM",  outcome_label="Type 2 Diabetes",       gene="PBX4",    rsid="rs73004967", chr=19, idx_pos=19717056),
  list(id="3_X11849_t2d",            metab="X - 11849",                   outcome_code="T2DM",  outcome_label="Type 2 Diabetes",       gene="COMT",    rsid="rs4633",     chr=22, idx_pos=19950235),
  list(id="4_X11849_hba1c",          metab="X - 11849",                   outcome_code="HBA1C", outcome_label="HbA1c Levels",          gene="COMT",    rsid="rs4633",     chr=22, idx_pos=19950235),
  list(id="5_hydroxypalmitate_fg",   metab="2-hydroxypalmitate",          outcome_code="FG",    outcome_label="Fasting Glucose Levels", gene="TMEM45A", rsid="rs59771628", chr=3,  idx_pos=100219086),
  list(id="6_trimethylurate_t2d",    metab="1,3,7-trimethylurate",        outcome_code="T2DM",  outcome_label="Type 2 Diabetes",       gene="GSTA5",   rsid="rs4144185",  chr=6,  idx_pos=52702948)
)

# LD r-squared to the index SNP.
ld_r2 <- function(chr, rsid) {
  tmp <- tempfile(tmpdir = WORK)
  args <- c(
    "--bfile", shQuote(EUR),
    "--chr", as.character(chr),
    "--ld-snp", rsid,
    "--r2",
    "--ld-window-kb", "2000",
    "--ld-window", "999999",
    "--ld-window-r2", "0",
    "--out", shQuote(tmp)
  )
  st <- system2(PLINK, args = args, stdout = FALSE, stderr = FALSE)
  stopifnot_loud(st == 0 && file.exists(paste0(tmp, ".ld")),
                 sprintf("PLINK --r2 failed for %s (chr%d).", rsid, chr))
  ld <- readr::read_table(paste0(tmp, ".ld"), col_types = cols(), progress = FALSE)
  dplyr::transmute(ld, pos = BP_B, r2 = R2)
}

build_side <- function(selected, cfg, side) {
  region_name <- if (identical(side, "metabolite")) {
    as.character(selected$eligibility$Metabolite_Region_File[[1]])
  } else {
    as.character(selected$eligibility$Outcome_Region_File[[1]])
  }
  ready_name <- if (identical(side, "metabolite")) {
    as.character(selected$ready$Metabolite_Input_File[[1]])
  } else {
    as.character(selected$ready$Outcome_Input_File[[1]])
  }
  stopifnot_loud(!is.na(region_name) && nzchar(region_name) && basename(region_name) == region_name,
    sprintf("%s: Stage 07 region filename is invalid.", cfg$id))
  stopifnot_loud(!is.na(ready_name) && nzchar(ready_name) && !grepl("^/|(^|/)\\.\\.(/|$)", ready_name),
    sprintf("%s: Stage 07 PwCoCo input path is invalid.", cfg$id))
  region_path <- file.path(REGIONS, region_name)
  ready_path <- file.path(COL, ready_name)
  stopifnot_loud(file.exists(region_path) && !dir.exists(region_path),
    sprintf("%s: Stage 07 region file is missing: %s", cfg$id, region_path))
  stopifnot_loud(file.exists(ready_path) && !dir.exists(ready_path),
    sprintf("%s: Stage 07 PwCoCo input is missing: %s", cfg$id, ready_path))
  region_root <- normalizePath(REGIONS, winslash = "/", mustWork = TRUE)
  ready_root <- normalizePath(COL, winslash = "/", mustWork = TRUE)
  stopifnot_loud(startsWith(normalizePath(region_path, winslash = "/", mustWork = TRUE), paste0(region_root, "/")),
    sprintf("%s: Stage 07 region file escapes its directory.", cfg$id))
  stopifnot_loud(startsWith(normalizePath(ready_path, winslash = "/", mustWork = TRUE), paste0(ready_root, "/")),
    sprintf("%s: Stage 07 PwCoCo input escapes its work directory.", cfg$id))
  region <- readr::read_tsv(region_path, col_types = cols(), progress = FALSE)
  pwcoco <- readr::read_tsv(ready_path, col_types = cols(), progress = FALSE)
  region_columns <- c("SNP", "Source_SNP", "Chromosome", "Position", "EffectAllele", "NonEffectAllele", "EAF", "Beta", "SE", "Pval", "N")
  ready_columns <- c("SNP", "A1", "A2", "A1_freq", "beta", "se", "p", "n")
  if (identical(side, "outcome") && identical(cfg$outcome_code, "T2DM")) ready_columns <- c(ready_columns, "ncase")
  stopifnot_loud(nrow(region) > 0L && nrow(pwcoco) > 0L && all(region_columns %in% names(region)) && all(ready_columns %in% names(pwcoco)),
    sprintf("%s: Stage 07 %s inputs have an invalid schema.", cfg$id, side))
  stopifnot_loud(!anyDuplicated(region$SNP) && !anyDuplicated(pwcoco$SNP),
    sprintf("%s: Stage 07 %s inputs contain duplicate SNPs.", cfg$id, side))
  ready_coordinate <- sub("_[ACGT]+_[ACGT]+$", "", as.character(pwcoco$SNP))
  stopifnot_loud(!any(ready_coordinate == as.character(pwcoco$SNP)) && !anyDuplicated(ready_coordinate),
    sprintf("%s: Stage 07 %s PwCoCo SNP identifiers are invalid.", cfg$id, side))
  keep <- as.character(region$SNP) %in% ready_coordinate
  stopifnot_loud(any(keep), sprintf("%s: Stage 07 %s region and PwCoCo input have no shared SNPs.", cfg$id, side))
  dplyr::transmute(region[keep, , drop = FALSE],
                    chrom = as.character(Chromosome), pos = as.numeric(Position),
                    snp = paste0("chr", Chromosome, ":", Position), source_snp = as.character(Source_SNP),
                    a1 = EffectAllele, a2 = NonEffectAllele,
                    b = as.numeric(Beta), se = as.numeric(SE), p = as.numeric(Pval), freq = as.numeric(EAF)) |>
    dplyr::distinct(pos, .keep_all = TRUE)
}

manifest <- list()
for (cfg in loci) {
  selected <- select_stage07_locus(cfg, eligibility, ready)
  message(sprintf(">> %s  (%s x %s, %s / %s)", cfg$id, cfg$metab, cfg$outcome_code, cfg$gene, cfg$rsid))
  ex <- build_side(selected, cfg, "metabolite")
  oc <- build_side(selected, cfg, "outcome")

  common <- intersect(ex$snp, oc$snp)
  stopifnot_loud(length(common) > 0L, sprintf("%s: exposure and outcome have no shared plotting SNPs.", cfg$id))
  ex <- ex[ex$snp %in% common, ]; oc <- oc[oc$snp %in% common, ]

  r2 <- ld_r2(cfg$chr, cfg$rsid)
  ex <- dplyr::left_join(ex, r2, by = "pos")
  oc <- dplyr::left_join(oc, r2, by = "pos")

  exposure_index <- ex$source_snp == cfg$rsid
  outcome_index <- oc$chrom == as.character(cfg$chr) & oc$pos == cfg$idx_pos
  stopifnot_loud(sum(exposure_index) == 1L && sum(outcome_index) == 1L,
    sprintf("%s: Stage 07 index SNP is not uniquely represented in plotting inputs.", cfg$id))
  ex$snp[exposure_index] <- cfg$rsid
  oc$snp[outcome_index] <- cfg$rsid
  ex$source_snp <- NULL
  oc$source_snp <- NULL

  write_tsv(ex, file.path(OUT, sprintf("%s_exposure.tsv", cfg$id)))
  write_tsv(oc, file.path(OUT, sprintf("%s_outcome.tsv",  cfg$id)))
  manifest[[length(manifest) + 1]] <- data.frame(
    id = cfg$id, metabolite = cfg$metab, outcome = cfg$outcome_label, gene = cfg$gene,
    index_rsid = cfg$rsid, chr = cfg$chr, index_pos = cfg$idx_pos,
    n_snps = length(common),
    n_r2_ge_0.8 = sum(ex$r2 >= 0.8, na.rm = TRUE), stringsAsFactors = FALSE)
}
man <- dplyr::bind_rows(manifest)
write_tsv(man, file.path(OUT, "manifest.tsv"))
message("\nLocus tables written to: ", OUT)
print(man)
