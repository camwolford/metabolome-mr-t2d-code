source("R/utils.R")

# For each metabolite-outcome association retained after the reverse-direction screen, take +-500 kb around
# every non-proxy instrument from the metabolite and outcome GWAS, harmonise the two, and write PwCoCo inputs.
paths <- data_paths(c("METABOLOME_MR_INPUT_DIR", "METABOLOME_MR_FULL_METABOLITE_GWAS_DIR", "METABOLOME_MR_WORK_DIR", "METABOLOME_MR_OUTPUT_DIR"))
sensitivity_dir <- file.path(paths[["output_dir"]], "05_sensitivity_analysis")
coloc_dir <- file.path(paths[["work_dir"]], "07_colocalisation")
outcomes <- data.frame(
  code = c("T2DM", "FG", "HBA1C"),
  prefix = c("t2dm", "fg", "hba1c"),
  gwas = c(
    file.path(paths[["input_dir"]], "t2dm_gwas_cleaned.tsv"),
    file.path(paths[["work_dir"]], "fasting_glucose_gwas_cleaned.tsv"),
    file.path(paths[["work_dir"]], "hbA1c_gwas_cleaned.tsv")
  )
)

# Keep one row per position: multi-allelic positions cannot be matched between studies by position.
unique_positions <- function(data) {
  key <- paste(data$Chromosome, data$Position)
  data[!(duplicated(key) | duplicated(key, fromLast = TRUE)), , drop = FALSE]
}

# Rows within +-500 kb of any instrument, as a filter for the chunked readers.
near_instruments <- function(instruments) function(d) {
  Reduce(`|`, lapply(seq_len(nrow(instruments)), function(j) {
    d$Chromosome == instruments$Chromosome[j] & abs(d$Position - instruments$Position[j]) <= 5e5
  }), FALSE)
}

# Stream the outcome GWAS in chunks, keeping only the instrument windows.
read_outcome_gwas <- function(path, keep) {
  columns <- readr::cols_only(Chromosome = "c", Position = "d", EffectAllele = "c", NonEffectAllele = "c", EAF = "d", Beta = "d",
                              SE = "d", Pval = "d", Ncases = "d", Ncontrols = "d", SampleSize = "d")
  filter_chunk <- function(chunk, index) {
    chunk <- chunk[!is.na(chunk$EAF) & !is.na(suppressWarnings(as.numeric(chunk$Chromosome))), ]  # autosomes only
    chunk$Chromosome <- as.numeric(chunk$Chromosome)
    chunk[keep(chunk), ]
  }
  gwas <- readr::read_tsv_chunked(path, readr::DataFrameCallback$new(filter_chunk), chunk_size = 1e6, col_types = columns)
  unique_positions(as.data.frame(gwas))
}

in_window <- function(data, chromosome, position) {
  data[data$Chromosome == chromosome & data$Position >= position - 5e5 & data$Position <= position + 5e5, , drop = FALSE]
}

# PwCoCo input with alphabetically ordered alleles in the SNP ID, matching the LD reference panel.
pwcoco_table <- function(snp, a1, a2, freq, beta, se, p, n) {
  flip <- a1 > a2
  data.frame(
    SNP = paste(snp, ifelse(flip, a2, a1), ifelse(flip, a1, a2), sep = "_"),
    A1 = ifelse(flip, a2, a1),
    A2 = ifelse(flip, a1, a2),
    A1_freq = ifelse(flip, 1 - freq, freq),
    beta = ifelse(flip, -beta, beta),
    se = se,
    p = p,
    n = n
  )
}

# Metabolites retained after the reverse-direction screen, and the instrument design used for each association.
screen <- readr::read_tsv(file.path(paths[["output_dir"]], "06_reverse_mr", "reverse_direction_screen.tsv"), show_col_types = FALSE)
associations <- do.call(rbind, lapply(outcomes$code, function(outcome) {
  do.call(rbind, lapply(c("conservative", "liberal"), function(design) {
    data <- readr::read_tsv(filtered_results_path(sensitivity_dir, design, outcome), show_col_types = FALSE)
    if (nrow(data)) data.frame(Metabolite = data$Metabolite, Outcome = outcome, Design = design)
  }))
}))
associations <- associations[associations$Metabolite %in% screen$Metabolite[screen$Retained], ]

comp_ids <- metabolite_comp_ids(paths[["input_dir"]])

loci <- list()
associations$Number_of_IVs <- NA_integer_
associations$Number_Assessed <- NA_integer_
for (outcome in outcomes$code) {
  outcome_info <- outcomes[outcomes$code == outcome, ]
  rows <- which(associations$Outcome == outcome)
  if (!length(rows)) next
  input_dir <- file.path(coloc_dir, "pwcoco_inputs", outcome)
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

  # Non-proxy instruments for each association (proxy instruments are not assessed).
  instrument_sets <- lapply(rows, function(i) {
    readr::read_tsv(instrument_files(harmonised_directory(paths[["output_dir"]], associations$Design[i], outcome))[[associations$Metabolite[i]]], show_col_types = FALSE)
  })
  associations$Number_of_IVs[rows] <- vapply(instrument_sets, nrow, integer(1))
  instrument_sets <- lapply(instrument_sets, function(x) x[is.na(x$proxy), ])
  associations$Number_Assessed[rows] <- vapply(instrument_sets, nrow, integer(1))
  outcome_gwas <- read_outcome_gwas(outcome_info$gwas, near_instruments(do.call(rbind, lapply(instrument_sets, `[`, c("Chromosome", "Position")))))

  for (k in seq_along(rows)) {
    association <- associations[rows[k], ]
    instruments <- instrument_sets[[k]]
    if (!nrow(instruments)) next
    metabolite_gwas <- unique_positions(read_metabolite_gwas(comp_ids[[association$Metabolite]], near_instruments(instruments)))

    for (j in seq_len(nrow(instruments))) {
      iv <- instruments[j, ]
      m <- in_window(metabolite_gwas, iv$Chromosome, iv$Position)
      o <- in_window(outcome_gwas, iv$Chromosome, iv$Position)
      exposure <- data.frame(
        SNP = paste(m$Chromosome, m$Position, sep = "_"), beta.exposure = m$Beta, se.exposure = m$SE,
        effect_allele.exposure = m$EffectAllele, other_allele.exposure = m$NonEffectAllele,
        eaf.exposure = m$EAF, pval.exposure = m$Pval, exposure = "metabolite", id.exposure = "metabolite"
      )
      outcome_data <- data.frame(
        SNP = paste(o$Chromosome, o$Position, sep = "_"), beta.outcome = o$Beta, se.outcome = o$SE,
        effect_allele.outcome = toupper(normalise_allele(o$EffectAllele)), other_allele.outcome = toupper(normalise_allele(o$NonEffectAllele)),
        eaf.outcome = o$EAF, pval.outcome = o$Pval, outcome = outcome, id.outcome = outcome
      )
      snv <- function(a) grepl("^[ACGT]$", a)
      exposure <- exposure[snv(exposure$effect_allele.exposure) & snv(exposure$other_allele.exposure), ]
      outcome_data <- outcome_data[snv(outcome_data$effect_allele.outcome) & snv(outcome_data$other_allele.outcome), ]
      h <- TwoSampleMR::harmonise_data(exposure, outcome_data, action = 2)
      h <- h[!h$remove, ]  # drop variants with incompatible alleles

      metabolite_n <- m$N[match(h$SNP, paste(m$Chromosome, m$Position, sep = "_"))]
      outcome_rows <- match(h$SNP, paste(o$Chromosome, o$Position, sep = "_"))
      metabolite_table <- pwcoco_table(h$SNP, h$effect_allele.exposure, h$other_allele.exposure, h$eaf.exposure,
                                       h$beta.exposure, h$se.exposure, h$pval.exposure, metabolite_n)
      if (outcome == "T2DM") {
        outcome_table <- pwcoco_table(h$SNP, h$effect_allele.outcome, h$other_allele.outcome, h$eaf.outcome,
                                      h$beta.outcome, h$se.outcome, h$pval.outcome, o$Ncases[outcome_rows] + o$Ncontrols[outcome_rows])
        outcome_table$ncase <- o$Ncases[outcome_rows]
      } else {
        outcome_table <- pwcoco_table(h$SNP, h$effect_allele.outcome, h$other_allele.outcome, h$eaf.outcome,
                                      h$beta.outcome, h$se.outcome, h$pval.outcome, o$SampleSize[outcome_rows])
      }

      metabolite_file <- paste0(association$Metabolite, "_", iv$SNP, "_metabolite_snps_pwcoco.tsv")
      outcome_file <- paste0(association$Metabolite, "_", iv$SNP, "_", outcome_info$prefix, "_snps_pwcoco.tsv")
      readr::write_tsv(metabolite_table, file.path(input_dir, metabolite_file))
      readr::write_tsv(outcome_table, file.path(input_dir, outcome_file))
      loci[[length(loci) + 1]] <- data.frame(
        Metabolite = association$Metabolite, Outcome = outcome, Design = association$Design,
        SNP = iv$SNP, Chromosome = iv$Chromosome, Position = iv$Position,
        Metabolite_File = metabolite_file, Outcome_File = outcome_file
      )
    }
  }
}

readr::write_tsv(associations, file.path(coloc_dir, "associations.tsv"))
readr::write_tsv(do.call(rbind, loci), file.path(coloc_dir, "loci.tsv"))
