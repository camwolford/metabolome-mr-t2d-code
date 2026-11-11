source("R/utils.R")

# Run PwCoCo for every instrument locus. PwCoCo appends each result to one .coloc file per outcome.
# METABOLOME_MR_PWCOCO is the PwCoCo executable; METABOLOME_MR_LD_BFILE is the LD reference PLINK prefix,
# with {chr} standing in for the chromosome number.
paths <- data_paths("METABOLOME_MR_WORK_DIR")
coloc_dir <- file.path(paths[["work_dir"]], "07_colocalisation")
loci <- readr::read_tsv(file.path(coloc_dir, "loci.tsv"), show_col_types = FALSE)

unlink(file.path(coloc_dir, sprintf("pwcoco_out_%s.coloc", c("t2dm", "fg", "hba1c"))))
for (i in seq_len(nrow(loci))) {
  locus <- loci[i, ]
  input_dir <- file.path(coloc_dir, "pwcoco_inputs", locus$Outcome)
  system2(Sys.getenv("METABOLOME_MR_PWCOCO"), c(
    "--bfile", shQuote(sub("{chr}", locus$Chromosome, Sys.getenv("METABOLOME_MR_LD_BFILE"), fixed = TRUE)),
    "--sum_stats1", shQuote(file.path(input_dir, locus$Metabolite_File)),
    "--sum_stats2", shQuote(file.path(input_dir, locus$Outcome_File)),
    "--out", shQuote(file.path(coloc_dir, paste0("pwcoco_out_", tolower(locus$Outcome)))),
    "--chr", locus$Chromosome,
    "--maf", "0.01"
  ))
}
