# Figure 2 volcano plots from external MR summaries.
# Fixed IVW is used for <=3 instruments; otherwise random IVW is used.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr)
  library(ggplot2); library(patchwork)
})

# Reproducible label placement
set.seed(20260707)

# External paths
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript.", call. = FALSE)
ARCHIVE_ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), "..", ".."), mustWork = TRUE)
source(file.path(ARCHIVE_ROOT, "config", "environment.R"))
paths <- archive_paths(c("METABOLOME_MR_INPUT_DIR", "METABOLOME_MR_OUTPUT_DIR"))
FORWARD <- file.path(paths[["output_dir"]], "03_forward_mr", "liberal")
SIGNIFICANT <- file.path(paths[["output_dir"]], "04_significance_filtering", "combined")
TABLES <- file.path(paths[["input_dir"]], "tables")
OUT <- file.path(paths[["output_dir"]], "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(ARCHIVE_ROOT, "figures", "R", "theme_cvdiab.R"))

# Plot constants
BONF      <- 0.05 / 1769                 # 2.826456e-05
YTHRESH   <- -log10(BONF)                # 4.5486
EXPECT    <- c(T2DM = 78L, FG = 46L, HBA1C = 47L)  # Expected coloured counts.
# Caps preserve readable panels; off-scale points are labelled at the edge.
YCAP      <- c(T2DM = 35, FG = 38, HBA1C = 13)
# Limits are symmetric about each outcome null.
XLIM      <- list(T2DM = c(0.38, 1.62), FG = c(-0.35, 0.35), HBA1C = c(-0.13, 0.13))
# Label the most significant on-scale coloured metabolites.
TOP_N     <- c(T2DM = 12L, FG = 8L, HBA1C = 8L)

# External data
full_lib <- read_tsv(file.path(FORWARD, "Full_MR_Results_Liberal.tsv"),
                     show_col_types = FALSE)
sig94    <- read_tsv(file.path(SIGNIFICANT, "Full_Significant_Results_Manuscript.tsv"),
                     show_col_types = FALSE)
tab2     <- read_csv(file.path(TABLES, "Partial Table 2.csv"), show_col_types = FALSE)

# Validate the externally authored confidence tiers.
conf_col <- grep("Confidence", names(tab2), value = TRUE)[1]
labels19 <- tibble(
  key     = gsub("[:/]", "_", tab2$Metabolite),
  conf    = trimws(tab2[[conf_col]])
)
stopifnot("all 19 prioritised names present in the 94-significant set" =
            all(labels19$key %in% sig94$Metabolite))
local({
  tiers <- table(labels19$conf)
  stopifnot("expected 5 high + 3 medium + 11 low confidence tiers" =
              tiers[["High"]] == 5L && tiers[["Medium"]] == 3L && tiers[["Low"]] == 11L)
})

# Restore lipid chain notation from file-safe metabolite names.
clean_name <- function(x) {
  x <- gsub("X - ", "X-", x, fixed = TRUE)
  # Alternate chain length and double-bond separators.
  reformat_grp <- function(sub) {
    if (!str_detect(sub, "[0-9]_[0-9]")) return(sub)
    m <- str_match(sub, "^([^0-9]*)(.+)$"); prefix <- m[1, 2]; toks <- str_split(m[1, 3], "_")[[1]]
    if (length(toks) < 2L) return(sub)
    out <- toks[1]
    for (i in 2:length(toks)) out <- paste0(out, if (i %% 2L == 0L) ":" else "/", toks[i])
    paste0(prefix, out)
  }
  str_replace_all(x, "\\(([^)]*)\\)", function(g) {
    inner <- str_sub(g, 2L, -2L)
    parts <- vapply(str_split(inner, ", ")[[1]], reformat_grp, character(1))
    paste0("(", paste(parts, collapse = ", "), ")")
  })
}

# Outcome-specific plotting frame
rule_p <- function(d, oc) {
  ivs  <- d[[paste0("Number_of_IVs_", oc)]]
  fixp <- d[[paste0("Fixed_IVW_Pval_", oc)]]
  ranp <- d[[paste0("Random_IVW_Pval_", oc)]]
  ifelse(ivs <= 3, fixp, ranp)
}

build_outcome <- function(oc) {
  est_c <- paste0("Fixed_IVW_Estimate_", oc)
  sg <- tibble(name = sig94$Metabolite, est = sig94[[est_c]], p = rule_p(sig94, oc)) |>
    filter(p < BONF) |>
    mutate(direction = if_else(est > 0, "risk", "protective"))
  bg <- tibble(name = full_lib$Metabolite, est = full_lib[[est_c]], p = rule_p(full_lib, oc)) |>
    filter(!is.na(p), !is.na(est), !name %in% sg$name) |>
    mutate(direction = "ns")
  # Non-concordant points remain grey.
  above <- bg |> filter(p < BONF)
  if (nrow(above) > 0)
    message(sprintf("[%s] %d point(s) clear the p-line but fail the robust filter (grey by design): %s",
                    oc, nrow(above), paste(above$name, collapse = "; ")))
  if (nrow(sg) != EXPECT[[oc]])
    stop(sprintf("[%s] significant count %d != expected %d", oc, nrow(sg), EXPECT[[oc]]))
  bind_rows(bg, sg)
}

# T2D uses odds ratios; glycaemic traits use beta estimates.
prep_panel <- function(df, oc) {
  cap  <- YCAP[[oc]]
  xlim <- XLIM[[oc]]
  d <- df |>
    mutate(
      x = if (oc == "T2DM") exp(est) else est,
      y = -log10(p),
      x_off = if (is.null(xlim)) FALSE else (x < xlim[1] | x > xlim[2]),
      y_off = y > cap,
      offscale = x_off | y_off,
      x_plot = if (is.null(xlim)) x else pmin(pmax(x, xlim[1]), xlim[2]),
      y_plot = pmin(y, cap),
      direction = factor(direction, levels = cvd_direction_levels),
      disp = str_wrap(clean_name(name), width = 20)
    )
  chosen <- d |>
    filter(direction != "ns", !offscale) |>
    slice_max(y, n = TOP_N[[oc]], with_ties = FALSE) |>
    pull(name)
  d |> mutate(lab = if_else(direction != "ns" & name %in% chosen, disp, NA_character_))
}

d_t2d <- prep_panel(build_outcome("T2DM"),  "T2DM")
d_fg  <- prep_panel(build_outcome("FG"),    "FG")
d_hb  <- prep_panel(build_outcome("HBA1C"), "HBA1C")

union_names <- unique(c(
  d_t2d$name[d_t2d$direction != "ns"],
  d_fg$name[d_fg$direction != "ns"],
  d_hb$name[d_hb$direction != "ns"]
))
if (length(union_names) != 94L)
  stop(sprintf("union of significant metabolites = %d, expected 94", length(union_names)))

# Panel builder
volcano_panel <- function(df, oc, xlab, null_x, title) {
  cap   <- YCAP[[oc]]
  xlim  <- XLIM[[oc]]
  xview <- if (is.null(xlim)) NULL else xlim + c(-1, 1) * diff(xlim) * 0.06   # small edge margin
  ns   <- df |> filter(direction == "ns", !offscale)
  hit  <- df |> filter(direction != "ns", !offscale)
  off  <- df |> filter(offscale) |>
    mutate(nm  = gsub("X - ", "X-", ifelse(is.na(disp), name, gsub("\n", " ", disp))),
           val = if (oc == "T2DM") sprintf("OR %.2f", x) else sprintf("effect %.2f", x),
           offlab = sprintf("%s (%s, p = %s)", nm, val, formatC(p, format = "e", digits = 1)))

  # Place all labels in one repel pass.
  lab_df <- bind_rows(
    hit |> filter(!is.na(lab)) |> transmute(x_plot, y_plot, labtext = lab,    face = "plain"),
    off |> transmute(x_plot, y_plot, labtext = offlab, face = "italic")
  )

  p <- ggplot(mapping = aes(x = x, y = y_plot)) +
    geom_vline(xintercept = null_x, linetype = "dashed",
               colour = cvd_pal$reference, linewidth = 0.35) +
    geom_hline(yintercept = YTHRESH, linetype = "dashed",
               colour = cvd_pal$reference, linewidth = 0.35) +
    geom_point(data = ns, aes(fill = direction, shape = direction),
               colour = "grey55", stroke = 0.1, size = 1.0, alpha = 0.45) +
    geom_point(data = hit, aes(fill = direction, shape = direction),
               colour = "grey20", stroke = 0.15, size = 1.7, alpha = 0.9) +
    cvd_direction_scales(which = c("fill", "shape")) +
    scale_x_continuous(expand = expansion(mult = 0.06)) +
    labs(x = xlab, y = expression(-log[10](italic(p))), title = title) +
    theme_cvdiab() +
    theme(legend.position = "none")

  # Mark off-scale points at the plot boundary.
  if (nrow(off) > 0) {
    p <- p +
      geom_point(data = off, aes(x = x_plot, y = y_plot, shape = direction, colour = direction),
                 fill = NA, size = 2.4, stroke = 0.6) +
      scale_colour_manual(values = c(risk = cvd_pal$risk, protective = cvd_pal$protective,
                                     ns = cvd_pal$ns), guide = "none")
  }
  p <- p + ggrepel::geom_text_repel(
    data = lab_df, aes(x = x_plot, y = y_plot, label = labtext),
    fontface = lab_df$face, size = 1.75, colour = cvd_pal$text, lineheight = 0.8,
    box.padding = 0.32, point.padding = 0.1, min.segment.length = 0.05,
    segment.size = 0.18, segment.colour = "grey60", segment.alpha = 0.85,
    max.overlaps = Inf, force = 2.6, force_pull = 0.85,
    max.iter = 80000, max.time = 4, seed = 20260707, na.rm = TRUE)

  p + coord_cartesian(xlim = xview, ylim = c(0, cap * 1.02), clip = "off")
}

pA <- volcano_panel(d_t2d, "T2DM",
                    xlab = "Odds ratio for Type 2 diabetes per SD",
                    null_x = 1,
                    title = "Type 2 diabetes")
pB <- volcano_panel(d_fg, "FG",
                    xlab = "Effect on Fasting glucose (SD units)",
                    null_x = 0,
                    title = "Fasting glucose")
pC <- volcano_panel(d_hb, "HBA1C",
                    xlab = "Effect on HbA1c (SD units)",
                    null_x = 0,
                    title = "HbA1c")

# Compose panels with a shared direction legend.
legend_src <- ggplot(d_t2d, aes(x, y_plot, fill = direction, shape = direction)) +
  geom_point(size = 1.9, colour = "grey20", stroke = 0.15) +
  cvd_direction_scales(which = c("fill", "shape")) +
  theme_cvdiab() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(size = 2.4)))
legend <- cowplot::get_legend(legend_src)

composite <- (pA / (pB | pC)) +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"),
        plot.margin = margin(3, 7, 4, 7))

final <- cowplot::plot_grid(composite, legend, ncol = 1, rel_heights = c(1, 0.045))

save_cvd_fig(final, file.path(OUT, "fig2_volcano.pdf"), width_mm = 185, height_mm = 205)

# Reconciliation output
cat("\n=== Figure 2 volcano — reconciliation ===\n")
cat(sprintf("Bonferroni threshold: 0.05/1769 = %.4e  (-log10 = %.3f)\n", BONF, YTHRESH))
cat(sprintf("Significant (coloured): T2D %d | FG %d | HbA1c %d | union %d\n",
            sum(d_t2d$direction != "ns"), sum(d_fg$direction != "ns"),
            sum(d_hb$direction != "ns"), length(union_names)))
cat(sprintf("Direction split  T2D risk/prot: %d/%d  FG: %d/%d  HbA1c: %d/%d\n",
            sum(d_t2d$direction=="risk"), sum(d_t2d$direction=="protective"),
            sum(d_fg$direction=="risk"),  sum(d_fg$direction=="protective"),
            sum(d_hb$direction=="risk"),  sum(d_hb$direction=="protective")))
cat(sprintf("Labelled prioritised hits per panel: T2D %d | FG %d | HbA1c %d\n",
            sum(!is.na(d_t2d$lab)), sum(!is.na(d_fg$lab)), sum(!is.na(d_hb$lab))))
named <- function(d) paste(gsub("\n", " ", na.omit(unique(d$lab))), collapse = "; ")
cat("  Named T2D:  ", named(d_t2d), "\n")
cat("  Named FG:   ", named(d_fg), "\n")
cat("  Named HbA1c:", named(d_hb), "\n")
cat(sprintf("Off-scale points (drawn at cap): T2D %d | FG %d | HbA1c %d\n",
            sum(d_t2d$offscale), sum(d_fg$offscale), sum(d_hb$offscale)))
sphg <- d_t2d |> filter(name == "sphinganine")
if (nrow(sphg) == 1)
  cat(sprintf("Spot-check sphinganine T2D OR = %.3f (Table 2: 1.595)\n", sphg$x))
cat("Wrote:", file.path(OUT, "fig2_volcano.pdf"), "\n")
