# Figure 4 forest plot from the authored Table 2 seed.
# Type 2 diabetes is plotted on an odds-ratio scale; glycaemic traits use beta estimates.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(purrr)
  library(ggplot2)
})
set.seed(20260707)

# External paths
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript.", call. = FALSE)
ARCHIVE_ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), "..", ".."), mustWork = TRUE)
source(file.path(ARCHIVE_ROOT, "config", "environment.R"))
paths <- archive_paths(c("METABOLOME_MR_INPUT_DIR", "METABOLOME_MR_OUTPUT_DIR"))
TABLES <- file.path(paths[["input_dir"]], "tables")
SIGNIFICANT <- file.path(paths[["output_dir"]], "04_significance_filtering", "combined")
OUT <- file.path(paths[["output_dir"]], "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(ARCHIVE_ROOT, "figures", "R", "theme_cvdiab.R"))

# En dash for displayed confidence intervals.
DASH <- "–"

OC_LAB <- c(T2D = "Type 2 diabetes", FG = "Fasting glucose", HbA1c = "HbA1c")
OC_UNIT <- c("Type 2 diabetes" = "Type 2 diabetes\n(OR per 1-SD metabolite)",
             "Fasting glucose" = "Fasting glucose\n(SD per 1-SD metabolite)",
             "HbA1c"           = "HbA1c\n(SD per 1-SD metabolite)")
OC_NULL <- c("Type 2 diabetes" = 1, "Fasting glucose" = 0, "HbA1c" = 0)
OC_CODE <- c(T2D = "T2DM", FG = "FG", HbA1c = "HBA1C")

# Parse Table 2
t2 <- read_csv(file.path(TABLES, "Partial Table 2.csv"), show_col_types = FALSE)
conf_col <- grep("Confidence", names(t2), value = TRUE)[1]
squish_ws <- function(x) str_squish(gsub("\n", " ", x))

t2 <- t2 |>
  mutate(disp = squish_ws(gsub("X - 11849", "X-11849", Metabolite, fixed = TRUE)),
         key  = gsub("[:/]", "_", Metabolite),
         tier = trimws(.data[[conf_col]]),
         est_str = squish_ws(`MR Estimate (95% CI)`),
         p_str   = squish_ws(`MR P-value`))

est_pat <- "(T2D|FG|HbA1c):\\s*(-?[0-9.]+)\\s*\\(\\s*(-?[0-9.]+)\\s*-\\s*(-?[0-9.]+)\\s*\\)"
p_pat   <- "(T2D|FG|HbA1c):\\s*([0-9.eE+-]+)"

parse_row <- function(i) {
  em <- str_match_all(t2$est_str[i], est_pat)[[1]]
  pm <- str_match_all(t2$p_str[i],  p_pat)[[1]]
  if (nrow(em) == 0) stop(sprintf("no estimate parsed for '%s'", t2$disp[i]))
  pvals <- setNames(suppressWarnings(as.numeric(pm[, 3])), pm[, 2])
  tibble(
    disp    = t2$disp[i], key = t2$key[i], tier = t2$tier[i],
    oc_code = em[, 2],
    est = as.numeric(em[, 3]), lo = as.numeric(em[, 4]), hi = as.numeric(em[, 5]),
    p   = unname(pvals[em[, 2]])
  )
}
long <- map_dfr(seq_len(nrow(t2)), parse_row) |>
  mutate(outcome = factor(OC_LAB[oc_code], levels = OC_LAB),
         null_x  = OC_NULL[as.character(outcome)],
         direction = if_else(est > null_x, "risk", "protective"))

stopifnot("expected 19 prioritised metabolites" = n_distinct(long$disp) == 19)

# Join instrument counts
sig <- read_tsv(file.path(SIGNIFICANT, "Full_Significant_Results_Manuscript.tsv"), show_col_types = FALSE)
iv_for <- function(key, oc_code) {
  col <- paste0("Number_of_IVs_", OC_CODE[[oc_code]])
  r <- which(sig$Metabolite == key)
  if (length(r) == 1 && col %in% names(sig)) sig[[col]][r] else NA_integer_
}
long <- long |>
  mutate(n_ivs = map2_int(key, oc_code, ~ as.integer(iv_for(.x, .y))),
         estimator = case_when(is.na(n_ivs) ~ "IVW",
                               n_ivs == 1   ~ "Wald",
                               n_ivs <= 3   ~ "fixed",
                               TRUE         ~ "random"),
         single_iv = !is.na(n_ivs) & n_ivs == 1)

# Order by confidence tier and type 2 diabetes distance from the null.
tier_lv <- c("High", "Medium", "Low")
sort_key <- long |>
  group_by(disp) |>
  summarise(tier = first(tier),
            t2d = { v <- est[outcome == "Type 2 diabetes"]; if (length(v)) v[1] else NA_real_ },
            anch = { v <- est[outcome == "Type 2 diabetes"]
                     if (length(v)) abs(log(v[1]))
                     else max(abs(est[outcome != "Type 2 diabetes"]), na.rm = TRUE) },
            .groups = "drop") |>
  mutate(tier = factor(tier, levels = tier_lv)) |>
  arrange(tier, desc(!is.na(t2d)), desc(anch))
metab_levels <- rev(sort_key$disp)
long <- long |>
  mutate(disp = factor(disp, levels = metab_levels),
         tier = factor(tier, levels = tier_lv))

tier_counts <- sort_key |> count(tier)
message(sprintf("Forest: 19 hits | tiers %s | single-IV rows: %d",
                paste(sprintf("%s=%d", tier_counts$tier, tier_counts$n), collapse = " "),
                sum(long$single_iv)))

# Forest rendering
build_ggplot <- function() {
  ref_df <- tibble(outcome = factor(OC_LAB, levels = OC_LAB), x0 = OC_NULL)
  strip_lab <- function(x) OC_UNIT[as.character(x)]

  ggplot(long, aes(x = est, y = disp)) +
    geom_vline(data = ref_df, aes(xintercept = x0),
               linetype = "dashed", colour = cvd_pal$reference, linewidth = 0.35) +
    geom_errorbarh(aes(xmin = lo, xmax = hi, colour = direction), height = 0.32, linewidth = 0.5) +
    geom_point(aes(colour = direction, fill = direction), shape = 21, size = 2.3, stroke = 0.3) +
    scale_colour_manual(values = c(risk = cvd_pal$risk, protective = cvd_pal$protective),
                        guide = "none") +
    scale_fill_manual(values = c(risk = cvd_pal$risk, protective = cvd_pal$protective),
                      guide = "none") +
    facet_grid(tier ~ outcome, scales = "free", space = "free_y",
               labeller = labeller(outcome = strip_lab,
                                   tier = c(High = "High\nconfidence",
                                            Medium = "Medium\nconfidence",
                                            Low = "Low\nconfidence"))) +
    labs(x = NULL, y = NULL) +
    theme_cvdiab() +
    theme(panel.spacing.x = unit(6, "mm"),
          panel.spacing.y = unit(1.2, "mm"),
          strip.text.x = element_text(size = rel(0.95), lineheight = 0.9),
          strip.text.y = element_text(size = rel(0.95), face = "bold", angle = 0),
          strip.background = element_rect(fill = "#F0F0F0", colour = NA),
          axis.text.y = element_text(size = rel(0.82)),
          axis.text.x = element_text(size = rel(1.0)),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
          legend.position = "none")
}

gg <- build_ggplot()
save_cvd_fig(gg, file.path(OUT, "fig4_forest.pdf"), width_mm = 200, height_mm = 190)
message("Wrote: ", file.path(OUT, "fig4_forest.pdf"))

# Instrument-count report for Table 2.
iv_table <- long |>
  arrange(disp, outcome) |>
  group_by(disp, tier) |>
  summarise(ivs = paste(sprintf("%s: %d", oc_code, n_ivs), collapse = ", "), .groups = "drop") |>
  mutate(disp = factor(disp, levels = rev(metab_levels))) |>
  arrange(disp)
cat("\n=== SNP-IV counts for Table 2 (top->bottom; add as one 'Number of SNP IVs' column) ===\n")
for (i in seq_len(nrow(iv_table)))
  cat(sprintf("  %-52s %s\n", as.character(iv_table$disp[i]), iv_table$ivs[i]))

# Reconciliation output
cat("\n=== Figure 4 forest — reconciliation vs Table 2 ===\n")
cat(sprintf("Prioritised hits: %d | tiers: %s\n", n_distinct(long$disp),
            paste(sprintf("%s=%d", tier_counts$tier, tier_counts$n), collapse = " ")))
cat(sprintf("Outcome rows: T2D %d | FG %d | HbA1c %d (total %d point estimates)\n",
            sum(long$oc_code == "T2D"), sum(long$oc_code == "FG"),
            sum(long$oc_code == "HbA1c"), nrow(long)))
cat(sprintf("Single-instrument (Wald) rows flagged: %d\n", sum(long$single_iv)))
spot <- long |> filter(disp == "sphinganine", oc_code == "T2D")
if (nrow(spot) == 1)
  cat(sprintf("Spot-check sphinganine T2D: OR %.3f (%.3f%s%.3f)  [Table 2: 1.595 (1.481-1.719)]\n",
              spot$est, spot$lo, DASH, spot$hi))
cat("Named hits, top->bottom:\n  ", paste(rev(metab_levels), collapse = "\n  "), "\n")
