# Supplementary Figure 1: build the six regional-plot tables from the PwCoCo colocalisation inputs,
# adding LD r2 with the lead variant from the 1000 Genomes European panel.

suppressPackageStartupMessages({
  library(readr)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
REPO_ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), "..", ".."), mustWork = TRUE)
source(file.path(REPO_ROOT, "R", "utils.R"))
paths <- data_paths(c("METABOLOME_MR_WORK_DIR", "METABOLOME_MR_EUR_PANEL_DIR"))
COLOC_DIR <- file.path(paths[["work_dir"]], "07_colocalisation")
OUT <- file.path(paths[["work_dir"]], "figures", "supp_fig1")
EUR <- file.path(paths[["eur_panel_dir"]], "EUR")
PLINK <- Sys.getenv("METABOLOME_MR_PLINK")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

loci_table <- read_tsv(file.path(COLOC_DIR, "loci.tsv"), col_types = cols())

# Six plotted loci
loci <- list(
  list(id="1_sphinganine_t2d",       metab="sphinganine",                outcome_code="T2DM",  outcome_label="Type 2 Diabetes",       gene="ABO",     rsid="rs676457",   chr=9,  idx_pos=136146227),
  list(id="2_palmitoleoylGPC_t2d",   metab="1-palmitoleoyl-GPC (16_1)*",  outcome_code="T2DM",  outcome_label="Type 2 Diabetes",       gene="PBX4",    rsid="rs73004967", chr=19, idx_pos=19717056),
  list(id="3_X11849_t2d",            metab="X - 11849",                   outcome_code="T2DM",  outcome_label="Type 2 Diabetes",       gene="COMT",    rsid="rs4633",     chr=22, idx_pos=19950235),
  list(id="4_X11849_hba1c",          metab="X - 11849",                   outcome_code="HBA1C", outcome_label="HbA1c Levels",          gene="COMT",    rsid="rs4633",     chr=22, idx_pos=19950235),
  list(id="5_hydroxypalmitate_fg",   metab="2-hydroxypalmitate",          outcome_code="FG",    outcome_label="Fasting Glucose Levels", gene="TMEM45A", rsid="rs59771628", chr=3,  idx_pos=100219086),
  list(id="6_trimethylurate_t2d",    metab="1,3,7-trimethylurate",        outcome_code="T2DM",  outcome_label="Type 2 Diabetes",       gene="GSTA5",   rsid="rs4144185",  chr=6,  idx_pos=52702948)
)

# LD r2 with the lead variant.
ld_r2 <- function(chr, rsid) {
  tmp <- tempfile(tmpdir = OUT)
  system2(PLINK, c("--bfile", shQuote(EUR), "--chr", chr, "--ld-snp", rsid, "--r2", "--ld-window-kb", "2000",
                   "--ld-window", "999999", "--ld-window-r2", "0", "--out", shQuote(tmp)), stdout = FALSE, stderr = FALSE)
  ld <- read_table(paste0(tmp, ".ld"), col_types = cols(), progress = FALSE)
  unlink(paste0(tmp, c(".ld", ".log", ".nosex")))
  dplyr::transmute(ld, pos = BP_B, r2 = R2)
}

# One side of a PwCoCo input (SNP IDs are chr_pos_A1_A2) as a plotting table.
read_side <- function(file, outcome) {
  d <- read_tsv(file.path(COLOC_DIR, "pwcoco_inputs", outcome, file), col_types = cols(), progress = FALSE)
  parts <- do.call(rbind, strsplit(d$SNP, "_", fixed = TRUE))
  data.frame(chrom = parts[, 1], pos = as.numeric(parts[, 2]), snp = paste0("chr", parts[, 1], ":", parts[, 2]),
             a1 = tolower(d$A1), a2 = tolower(d$A2), b = d$beta, se = d$se, p = d$p, freq = d$A1_freq)
}

manifest <- list()
for (cfg in loci) {
  locus <- loci_table[loci_table$Metabolite == cfg$metab & loci_table$Outcome == cfg$outcome_code & loci_table$SNP == cfg$rsid, ]
  ex <- read_side(locus$Metabolite_File, cfg$outcome_code)
  oc <- read_side(locus$Outcome_File, cfg$outcome_code)
  r2 <- ld_r2(cfg$chr, cfg$rsid)
  ex <- dplyr::left_join(ex, r2, by = "pos")
  oc <- dplyr::left_join(oc, r2, by = "pos")
  ex$snp[ex$pos == cfg$idx_pos] <- cfg$rsid
  oc$snp[oc$pos == cfg$idx_pos] <- cfg$rsid
  write_tsv(ex, file.path(OUT, sprintf("%s_exposure.tsv", cfg$id)))
  write_tsv(oc, file.path(OUT, sprintf("%s_outcome.tsv", cfg$id)))
  manifest[[length(manifest) + 1]] <- data.frame(
    id = cfg$id, metabolite = cfg$metab, outcome = cfg$outcome_label, gene = cfg$gene,
    index_rsid = cfg$rsid, chr = cfg$chr, index_pos = cfg$idx_pos, n_snps = nrow(ex),
    n_r2_ge_0.8 = sum(ex$r2 >= 0.8, na.rm = TRUE))
}
write_tsv(dplyr::bind_rows(manifest), file.path(OUT, "manifest.tsv"))
