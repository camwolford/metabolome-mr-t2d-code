# Supplementary Figure 1 locus plots from external region and colocalisation data.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr)
  library(locuszoomr); library(cowplot); library(ggplot2)
  library(EnsDb.Hsapiens.v75)
})

set.seed(20260717)

edb <- EnsDb.Hsapiens.v75

# External paths
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript.", call. = FALSE)
REPO_ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), "..", ".."), mustWork = TRUE)
source(file.path(REPO_ROOT, "R", "utils.R"))
paths <- data_paths(c("METABOLOME_MR_INPUT_DIR", "METABOLOME_MR_WORK_DIR", "METABOLOME_MR_OUTPUT_DIR"))
REG <- file.path(paths[["work_dir"]], "figures", "supp_fig1")
WORKBOOK <- file.path(paths[["input_dir"]], "tables", "Metabolome-wide MR Tables.xlsx")
OUT <- file.path(paths[["output_dir"]], "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

manifest <- read_tsv(file.path(REG, "manifest.tsv"), col_types = cols())

# Display labels
outcome_label <- c("Type 2 Diabetes"        = "Type 2 diabetes",
                   "Fasting Glucose Levels" = "Fasting glucose",
                   "HbA1c Levels"           = "HbA1c")
# PwCoCo results are Supplementary Tables 13-15 of the workbook.
coloc_sheet <- c("Type 2 Diabetes"        = "Supplementary Table 13",
                 "Fasting Glucose Levels" = "Supplementary Table 14",
                 "HbA1c Levels"           = "Supplementary Table 15")
# Restore display formatting from file-safe metabolite names.
pretty_metab <- function(x) {
  x <- str_replace(x, "^X - 11849$", "X-11849")
  str_replace(x, "\\((\\d+)_(\\d+)\\)", "(\\1:\\2)")
}
# Prefer an unconditioned colocalisation row; otherwise use the highest-H4 row.
norm_name <- function(x) gsub("[:/_ ]", "", x)  # match file-safe and workbook names
get_pp <- function(metab, outcome, rsid) {
  cf <- readxl::read_excel(WORKBOOK, sheet = coloc_sheet[[outcome]])
  rows <- cf %>% dplyr::filter(norm_name(Metabolite) == norm_name(metab),
                               `SNP IV` == rsid)
  if (nrow(rows) == 0)
    stop(sprintf("Supp Fig 1: no coloc row for %s / %s / %s", metab, outcome, rsid))
  uncond <- rows %>% dplyr::filter(`Condition SNP 1` == "unconditioned",
                                   `Condition SNP 2` == "unconditioned")
  row <- if (nrow(uncond) >= 1) uncond[1, ]
         else dplyr::slice(dplyr::arrange(rows, dplyr::desc(H4)), 1)
  cond <- !(row$`Condition SNP 1`[1] == "unconditioned" &&
            row$`Condition SNP 2`[1] == "unconditioned")
  c(h4 = row$H4[1], h3 = row$H3[1], nsnps = row$nsnps[1], conditioned = as.numeric(cond))
}

# Locus panel from one external region table
make_locus <- function(id, side, rsid) {
  f  <- file.path(REG, sprintf("%s_%s.tsv", id, side))
  df <- read_tsv(f, col_types = cols())
  locus(data = df, index_snp = rsid, LD = "r2", ens_db = edb, fix_window = 1e6,
        chrom = "chrom", pos = "pos", p = "p", labs = "snp")
}

# Name Helvetica explicitly: quartz renders an unnamed family as Arial.
for (g in list("text", "label", ggrepel::GeomTextRepel))
  update_geom_defaults(g, list(family = "Helvetica"))
helv <- theme(text = element_text(family = "Helvetica"))

strip <- function(g, txt) g + ggtitle(txt) + helv +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0,
                                  margin = margin(b = 2)))

# Extend the y-axis so index-SNP labels stay inside the panel.
label_ylim <- function(loc) {
  y <- loc$data[[loc$yvar]]
  c(min(y, na.rm = TRUE), max(y, na.rm = TRUE) + 0.7)
}

# One locus block
build_block <- function(m, letter, show_legend) {
  metab_pretty <- pretty_metab(m$metabolite)
  out_pretty   <- outcome_label[[m$outcome]]
  pp  <- get_pp(m$metabolite, m$outcome, m$index_rsid)
  ppf <- function(x) if (is.na(x)) "NA" else if (x < 0.001) "< 0.001" else sprintf("%.2f", x)

  loc_ex <- make_locus(m$id, "exposure", m$index_rsid)
  loc_oc <- make_locus(m$id, "outcome",  m$index_rsid)

  pA <- strip(gg_scatter(loc_ex, labels = "index", cex.axis = 0.9, cex.lab = 1.0,
                         ylim = label_ylim(loc_ex),
                         legend_pos = if (show_legend) "topleft" else NULL,
                         nudge_x = 0.025, nudge_y = 0.35, hjust = 0, vjust = 0,
                         force = 2, box.padding = 0.25, max.overlaps = Inf,
                         seed = 20260717),
              sprintf("Metabolite: %s", metab_pretty))
  pB <- strip(gg_scatter(loc_oc, labels = "index", cex.axis = 0.9, cex.lab = 1.0,
                         ylim = label_ylim(loc_oc),
                         legend_pos = NULL,
                         nudge_x = 0.025, nudge_y = 0.35, hjust = 0, vjust = 0,
                         force = 2, box.padding = 0.25, max.overlaps = Inf,
                         seed = 20260717),
              sprintf("Outcome: %s", out_pretty))
  condtxt <- if (isTRUE(pp["conditioned"] == 1)) "\n(conditional analysis)" else ""
  pB <- pB + annotate("label", x = Inf, y = Inf, hjust = 1.04, vjust = 1.1,
                      label = sprintf("%d SNPs\nPP.H4 = %s\nPP.H3 = %s%s",
                                      as.integer(pp["nsnps"]), ppf(pp["h4"]), ppf(pp["h3"]), condtxt),
                      size = 2.6, label.size = 0.25, fill = "grey96", lineheight = 0.95)
  pG <- gg_genetracks(loc_ex, italics = TRUE, filter_gene_biotype = "protein_coding",
                      cex.text = 0.6, cex.axis = 0.85, cex.lab = 0.95, maxrows = 4) + helv

  body <- plot_grid(pA, pB, pG, ncol = 1, rel_heights = c(1, 1, 0.55),
                    align = "v", axis = "lr")

  title <- bquote(atop("Colocalisation of variant" ~ .(m$index_rsid) ~ "in the" ~
                         italic(.(m$gene)) ~ "region",
                       "between" ~ .(metab_pretty) ~ "and" ~ .(out_pretty)))
  hdr <- ggdraw() +
    draw_label(letter, x = 0.01, y = 0.5, fontface = "bold", size = 14, hjust = 0,
               fontfamily = "Helvetica") +
    draw_label(title, x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
               size = 10, lineheight = 0.78, fontfamily = "Helvetica")

  plot_grid(hdr, body, ncol = 1, rel_heights = c(0.1, 1))
}

# Render three two-locus pages
pdf_out <- file.path(OUT, "supp_fig1_locuszoom.pdf")
# macOS quartz PDF keeps bold/italic faces; base pdf() + pdftocairo flattens them, so it is
# only the fallback off macOS.
use_quartz <- isTRUE(capabilities("aqua"))
if (use_quartz) {
  grDevices::quartz(file = pdf_out, type = "pdf", width = 13.5, height = 9.6, family = "Helvetica")
} else {
  grDevices::pdf(pdf_out, width = 13.5, height = 9.6, onefile = TRUE)
}
letters6 <- LETTERS[seq_len(nrow(manifest))]
pages <- list(c(1, 2), c(3, 4), c(5, 6))
for (pg in pages) {
  message(sprintf(">> page: %s", paste(letters6[pg], collapse = "|")))
  blocks <- lapply(pg, function(i)
    build_block(manifest[i, ], letters6[i], show_legend = (i == pg[1])))
  print(plot_grid(plotlist = blocks, ncol = 2))
}
invisible(dev.off())

# Embed fonts when pdftocairo is available (base pdf() path only).
if (!use_quartz) local({
  ptc <- Sys.which("pdftocairo")
  if (nzchar(ptc)) {
    tmp <- paste0(pdf_out, ".embed.tmp.pdf")
    rc <- tryCatch(system2(ptc, c("-pdf", shQuote(pdf_out), shQuote(tmp)),
                           stdout = FALSE, stderr = FALSE), error = function(e) 1L)
    if (isTRUE(rc == 0) && file.exists(tmp) && file.info(tmp)$size > 2000)
      file.rename(tmp, pdf_out) else unlink(tmp)
  }
})
message("\nWrote ", pdf_out)
