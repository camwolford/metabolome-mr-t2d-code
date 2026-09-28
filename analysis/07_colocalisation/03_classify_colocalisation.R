source("R/utils.R")

# An instrument colocalises if any PwCoCo configuration reaches PP.H4 >= 0.80; conditioned configurations
# also need at least 100 variants. An association colocalises if all its assessed instruments do, and
# partially colocalises if only some do.
paths <- data_paths(c("METABOLOME_MR_WORK_DIR", "METABOLOME_MR_OUTPUT_DIR"))
coloc_dir <- file.path(paths[["work_dir"]], "07_colocalisation")
output_dir <- file.path(paths[["output_dir"]], "07_colocalisation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

loci <- as.data.frame(readr::read_tsv(file.path(coloc_dir, "loci.tsv"), show_col_types = FALSE))
associations <- as.data.frame(readr::read_tsv(file.path(coloc_dir, "associations.tsv"), show_col_types = FALSE))

coloc <- do.call(rbind, lapply(c("T2DM", "FG", "HBA1C"), function(outcome) {
  path <- file.path(coloc_dir, sprintf("pwcoco_out_%s.coloc", tolower(outcome)))
  if (!file.exists(path)) return(NULL)
  data <- as.data.frame(readr::read_tsv(path, show_col_types = FALSE))
  data$Outcome <- outcome
  data
}))
conditioned <- !(coloc$SNP1 == "unconditioned" & coloc$SNP2 == "unconditioned")
coloc$Pass <- coloc$H4 >= 0.80 & (!conditioned | coloc$nsnps >= 100)

rows_for <- function(locus) coloc[coloc$Outcome == locus$Outcome & basename(coloc$Dataset1) == locus$Metabolite_File, , drop = FALSE]
loci$Max_H4 <- vapply(seq_len(nrow(loci)), function(i) max(c(rows_for(loci[i, ])$H4, NA), na.rm = TRUE), numeric(1))
loci$Colocalised <- vapply(seq_len(nrow(loci)), function(i) any(rows_for(loci[i, ])$Pass), logical(1))

key <- function(data) paste(data$Metabolite, data$Outcome)
associations$Number_Colocalised <- vapply(key(associations), function(k) sum(loci$Colocalised[key(loci) == k]), integer(1))
associations$Colocalisation <- ifelse(
  associations$Number_Assessed == 0, "not_assessed",
  ifelse(associations$Number_Colocalised == associations$Number_Assessed, "colocalised",
         ifelse(associations$Number_Colocalised > 0, "partial", "not_colocalised"))
)

readr::write_tsv(loci, file.path(output_dir, "coloc_loci.tsv"))
readr::write_tsv(associations, file.path(output_dir, "coloc_associations.tsv"))
