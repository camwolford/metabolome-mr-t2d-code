# Figure 4 forest plot from the published Table 2 (Additional file workbook).
# Type 2 diabetes is plotted on an odds-ratio scale; glycaemic traits use beta estimates.
# Rows keep Table 2 order and are grouped by colocalisation locus.

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(stringr); library(purrr)
  library(ggplot2)
})
set.seed(20260707)

# External paths
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript.", call. = FALSE)
REPO_ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), "..", ".."), mustWork = TRUE)
source(file.path(REPO_ROOT, "R", "utils.R"))
paths <- data_paths(c("METABOLOME_MR_INPUT_DIR", "METABOLOME_MR_OUTPUT_DIR"))
WORKBOOK <- file.path(paths[["input_dir"]], "tables", "Metabolome-wide MR Tables.xlsx")
OUT <- file.path(paths[["output_dir"]], "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(REPO_ROOT, "figures", "R", "theme.R"))

OC_LAB <- c(T2D = "Type 2 diabetes", FG = "Fasting glucose", HbA1c = "HbA1c")
OC_UNIT <- c("Type 2 diabetes" = "Type 2 diabetes\n(OR per 1-SD metabolite)",
             "Fasting glucose" = "Fasting glucose\n(mmol/l per 1-SD metabolite)",
             "HbA1c"           = "HbA1c\n(% per 1-SD metabolite)")
OC_NULL <- c("Type 2 diabetes" = 1, "Fasting glucose" = 0, "HbA1c" = 0)

# Parse Table 2; footnote rows have no locus group.
if (!file.exists(WORKBOOK)) stop(sprintf("Table workbook is missing: %s", WORKBOOK), call. = FALSE)
t2 <- read_excel(WORKBOOK, sheet = "Table 2")
required <- c("Metabolite", "Colocalisation locus group", "MR Estimate (95% CI)", "Number of SNP IVs")
absent <- setdiff(required, names(t2))
if (length(absent)) stop(sprintf("Table 2 is missing: %s", paste(absent, collapse = ", ")), call. = FALSE)
squish_ws <- function(x) str_squish(gsub("\n", " ", x))

t2 <- t2 |>
  filter(!is.na(`Colocalisation locus group`)) |>
  mutate(disp = squish_ws(gsub("X - 11849", "X-11849", Metabolite, fixed = TRUE)),
         locus = `Colocalisation locus group`,
         est_str = squish_ws(`MR Estimate (95% CI)`),
         iv_str  = squish_ws(`Number of SNP IVs`))
if (nrow(t2) != 19L || anyDuplicated(t2$disp)) stop("Table 2 must list 19 unique metabolites.", call. = FALSE)

est_pat <- "(T2D|FG|HbA1c):\\s*(-?[0-9.]+)\\s*\\(\\s*(-?[0-9.]+)\\s*-\\s*(-?[0-9.]+)\\s*\\)"
iv_pat  <- "(T2D|FG|HbA1c):\\s*([0-9]+)"

parse_row <- function(i) {
  em <- str_match_all(t2$est_str[i], est_pat)[[1]]
  im <- str_match_all(t2$iv_str[i], iv_pat)[[1]]
  if (nrow(em) == 0) stop(sprintf("no estimate parsed for '%s'", t2$disp[i]), call. = FALSE)
  n_ivs <- setNames(as.integer(im[, 3]), im[, 2])
  tibble(disp = t2$disp[i], locus = t2$locus[i], table_order = i,
         oc_code = em[, 2],
         est = as.numeric(em[, 3]), lo = as.numeric(em[, 4]), hi = as.numeric(em[, 5]),
         n_ivs = unname(n_ivs[em[, 2]]))
}
long <- map_dfr(seq_len(nrow(t2)), parse_row) |>
  mutate(outcome = factor(OC_LAB[oc_code], levels = OC_LAB),
         null_x  = OC_NULL[as.character(outcome)],
         direction = if_else(est > null_x, "risk", "protective"),
         single_iv = n_ivs == 1L)

stopifnot("expected 19 prioritised metabolites" = n_distinct(long$disp) == 19)
stopifnot("expected 28 final associations" = nrow(long) == 28)
stopifnot("expected every final estimate to use one variant" = all(long$single_iv %in% TRUE))

# Row order: Table 2 order, first row at the top; locus groups must be contiguous.
metab_levels <- rev(t2$disp)
locus_runs <- rle(t2$locus)
if (anyDuplicated(locus_runs$values)) stop("Locus groups are not contiguous in Table 2 order.", call. = FALSE)
long <- long |>
  mutate(disp = factor(disp, levels = metab_levels),
         locus = factor(locus, levels = locus_runs$values))

message(sprintf("Forest: 19 hits | %d locus groups | single-IV rows: %d",
                length(locus_runs$values), sum(long$single_iv)))

# Forest rendering
build_ggplot <- function() {
  ref_df <- tibble(outcome = factor(OC_LAB, levels = OC_LAB), x0 = OC_NULL)
  strip_lab <- function(x) OC_UNIT[as.character(x)]

  ggplot(long, aes(x = est, y = disp)) +
    geom_vline(data = ref_df, aes(xintercept = x0),
               linetype = "dashed", colour = fig_colours$reference, linewidth = 0.35) +
    geom_errorbarh(aes(xmin = lo, xmax = hi, colour = direction), height = 0.32, linewidth = 0.5) +
    geom_point(aes(colour = direction, fill = direction), shape = 21, size = 2.4, stroke = 0.3) +
    scale_colour_manual(values = c(risk = fig_colours$risk, protective = fig_colours$protective),
                        guide = "none") +
    scale_fill_manual(values = c(risk = fig_colours$risk, protective = fig_colours$protective),
                      guide = "none") +
    facet_grid(locus ~ outcome, scales = "free", space = "free_y", switch = "y",
               labeller = labeller(outcome = strip_lab,
                                   locus = function(x) sub("–", "–\n", x, fixed = TRUE))) +
    scale_x_continuous(breaks = scales::breaks_extended(n = 4)) +
    labs(x = NULL, y = NULL) +
    theme_paper() +
    theme(panel.spacing.x = unit(6, "mm"),
          panel.spacing.y = unit(1.2, "mm"),
          strip.text.x = element_text(size = rel(0.95), lineheight = 0.9),
          strip.background = element_rect(fill = "#F0F0F0", colour = NA),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 0.5, face = "italic", size = rel(0.8),
                                           lineheight = 0.9),
          strip.background.y = element_rect(fill = "#F0F0F0", colour = NA),
          axis.text.y = element_text(size = rel(0.82)),
          axis.text.x = element_text(size = rel(1.0)),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
          legend.position = "none")
}

gg <- build_ggplot()
save_figure(gg, file.path(OUT, "fig4_forest.pdf"), width_mm = 200, height_mm = 190)
message("Wrote: ", file.path(OUT, "fig4_forest.pdf"))

# Reconciliation output
cat("\n=== Figure 4 forest — reconciliation vs Table 2 ===\n")
cat(sprintf("Outcome rows: T2D %d | FG %d | HbA1c %d (total %d point estimates)\n",
            sum(long$oc_code == "T2D"), sum(long$oc_code == "FG"),
            sum(long$oc_code == "HbA1c"), nrow(long)))
cat("Locus groups, top->bottom:\n  ", paste(locus_runs$values, collapse = "\n  "), "\n")
