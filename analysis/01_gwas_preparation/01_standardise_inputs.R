source("R/utils.R")

paths <- data_paths(c(
  "METABOLOME_MR_INPUT_DIR",
  "METABOLOME_MR_WORK_DIR",
  "METABOLOME_MR_OUTPUT_DIR"
))

input_dir <- paths[["input_dir"]]
work_dir <- paths[["work_dir"]]

metabolite_input <- file.path(input_dir, "metabolite_sup_associations.txt")
fasting_glucose_input <- file.path(input_dir, "fasting_glucose_gwas_cleaned.csv")
hba1c_input <- file.path(input_dir, "hbA1c_gwas_cleaned.csv")

metabolite_gwas <- utils::read.delim(metabolite_input, stringsAsFactors = FALSE)

colnames(metabolite_gwas) <- metabolite_gwas[1, ]
metabolite_gwas <- metabolite_gwas[-1, ]
metabolite_gwas <- metabolite_gwas[, -c(
  1, 2, 3, 11, 12, 13, 14, 15, 17, 18, 19, 20, 22, 23, 24,
  25, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41
)]

colnames(metabolite_gwas)[4] <- "SNP"
colnames(metabolite_gwas)[5] <- "EffectAllele"
colnames(metabolite_gwas)[6] <- "NonEffectAllele"
colnames(metabolite_gwas)[7] <- "Metabolite"
colnames(metabolite_gwas)[8] <- "EAF"
colnames(metabolite_gwas)[9] <- "Beta"
colnames(metabolite_gwas)[10] <- "SE"
colnames(metabolite_gwas)[11] <- "Pval"

metabolite_gwas$Chromosome <- as.numeric(metabolite_gwas$Chromosome)
metabolite_gwas$Position <- as.numeric(metabolite_gwas$Position)
metabolite_gwas$EAF <- as.numeric(metabolite_gwas$EAF)
metabolite_gwas$Beta <- as.numeric(metabolite_gwas$Beta)
metabolite_gwas$SE <- as.numeric(metabolite_gwas$SE)
metabolite_gwas$Pval <- as.numeric(metabolite_gwas$Pval)
metabolite_gwas$Pval_new <- 10^(-1 * metabolite_gwas$Pval)
metabolite_gwas$Pval <- metabolite_gwas$Pval_new
metabolite_gwas <- metabolite_gwas[, -12]
metabolite_gwas <- metabolite_gwas[stats::complete.cases(metabolite_gwas), ]
metabolite_gwas$Metabolite <- gsub(":", "_", metabolite_gwas$Metabolite)

metabolite_output <- file.path(work_dir, "metabolite_gwas_associations_cleaned.tsv")
utils::write.table(metabolite_gwas, metabolite_output, sep = "\t", row.names = FALSE)

fasting_glucose_gwas <- utils::read.csv(fasting_glucose_input, stringsAsFactors = FALSE)
colnames(fasting_glucose_gwas) <- c(
  "Chromosome", "Position", "EffectAllele", "NonEffectAllele", "EAF", "Beta", "SE", "Pval", "SampleSize"
)

fasting_glucose_output <- file.path(work_dir, "fasting_glucose_gwas_cleaned.tsv")
utils::write.table(fasting_glucose_gwas, fasting_glucose_output, sep = "\t", row.names = FALSE)

hba1c_gwas <- utils::read.csv(hba1c_input, stringsAsFactors = FALSE)
colnames(hba1c_gwas) <- c(
  "Chromosome", "Position", "EffectAllele", "NonEffectAllele", "EAF", "Beta", "SE", "Pval", "SampleSize"
)

hba1c_output <- file.path(work_dir, "hbA1c_gwas_cleaned.tsv")
utils::write.table(hba1c_gwas, hba1c_output, sep = "\t", row.names = FALSE)
