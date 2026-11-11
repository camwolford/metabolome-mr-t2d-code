source("R/utils.R")

paths <- data_paths(c("METABOLOME_MR_WORK_DIR", "METABOLOME_MR_EUR_PANEL_DIR"))
plink_bin <- Sys.getenv("METABOLOME_MR_PLINK")
work_dir <- paths[["work_dir"]]
eur_bfile <- file.path(paths[["eur_panel_dir"]], "EUR")

# For each variant, count the metabolites associated (at the instrument set's p threshold) with the variant
# itself or with a variant in LD with it (r^2 >= 0.8 in 1000 Genomes EUR, within 10 Mb, the clumping window).
metabolite_gwas <- readr::read_tsv(file.path(work_dir, "metabolite_gwas_associations_cleaned.tsv"), show_col_types = FALSE)
bim_snps <- utils::read.table(paste0(eur_bfile, ".bim"), colClasses = c("NULL", "character", rep("NULL", 4)))[[1]]

snp_file <- tempfile()
writeLines(intersect(unique(metabolite_gwas$SNP), bim_snps), snp_file)
ld_prefix <- tempfile()
status <- system2(plink_bin, c(
  "--bfile", shQuote(eur_bfile), "--extract", snp_file, "--r2",
  "--ld-window-kb", 10000, "--ld-window", 99999, "--ld-window-r2", 0.8, "--out", ld_prefix
), stdout = FALSE)
if (status != 0) stop("plink LD calculation failed")
ld <- utils::read.table(paste0(ld_prefix, ".ld"), header = TRUE, stringsAsFactors = FALSE)
ld_pairs <- rbind(data.frame(SNP = ld$SNP_A, partner = ld$SNP_B), data.frame(SNP = ld$SNP_B, partner = ld$SNP_A))

for (design in c("conservative", "liberal")) {
  p_threshold <- c(conservative = 5e-8, liberal = 5e-6)[[design]]
  associations <- unique(metabolite_gwas[metabolite_gwas$Pval < p_threshold, c("SNP", "Metabolite")])
  links <- unique(rbind(
    data.frame(SNP = unique(associations$SNP), partner = unique(associations$SNP)),
    ld_pairs[ld_pairs$SNP %in% associations$SNP, , drop = FALSE]
  ))
  linked <- merge(links, associations, by.x = "partner", by.y = "SNP")
  counts <- stats::aggregate(Metabolite ~ SNP, linked, function(x) length(unique(x)))
  names(counts)[2] <- "n_metab"
  readr::write_tsv(counts, file.path(work_dir, paste0("LD_Linked_Counts_", design, ".tsv")))
}
