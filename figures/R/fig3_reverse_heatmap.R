# Figure 3 reverse-MR heatmap from external Stage 06 summaries.
# Filled and ring markers denote nominal IVW and MR-Egger signals.

suppressPackageStartupMessages({
  library(readr); library(dplyr)
  library(ggplot2)
})

# Reproducibility
set.seed(20260707)

# External paths
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript.", call. = FALSE)
ARCHIVE_ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), "..", ".."), mustWork = TRUE)
source(file.path(ARCHIVE_ROOT, "config", "environment.R"))
paths <- archive_paths("METABOLOME_MR_OUTPUT_DIR")
REVERSE_RESULTS <- file.path(paths[["output_dir"]], "06_reverse_mr", "results")
ANNOTATION_FILE <- file.path(paths[["output_dir"]], "08_pathways_and_prioritisation",
                             "significant_metabolites_data.tsv")
OUT <- file.path(paths[["output_dir"]], "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(ARCHIVE_ROOT, "figures", "R", "theme_cvdiab.R"))

# Outcome metadata
OC <- tibble(
  file  = c("t2dm_reverse_mr_raw.tsv", "fg_reverse_mr_raw.tsv", "hba1c_reverse_mr_raw.tsv"),
  code  = c("T2DM", "FG", "HBA1C"),
  label = c("Type 2 diabetes", "Fasting glucose", "HbA1c")
)
OC$label <- factor(OC$label, levels = OC$label)

select_reverse_ivw <- function(data, label) {
  required <- c(
    "Metabolite", "Number_of_IVs", "Fixed_IVW_Estimate", "Fixed_IVW_Pval",
    "Random_IVW_Estimate", "Random_IVW_Pval", "Egger_Pval"
  )
  absent <- setdiff(required, names(data))
  if (length(absent)) {
    stop(sprintf("%s is missing required columns: %s", label, paste(absent, collapse = ", ")), call. = FALSE)
  }
  n_ivs <- suppressWarnings(as.numeric(data$Number_of_IVs))
  if (any(!is.finite(n_ivs) | n_ivs <= 0 | n_ivs != floor(n_ivs))) {
    stop(sprintf("%s has invalid positive integer Number_of_IVs values.", label), call. = FALSE)
  }
  estimate_columns <- c("Fixed_IVW_Estimate", "Random_IVW_Estimate")
  if (any(vapply(estimate_columns, function(column) any(!is.finite(as.numeric(data[[column]]))), logical(1)))) {
    stop(sprintf("%s has non-finite IVW estimates.", label), call. = FALSE)
  }
  p_columns <- c("Fixed_IVW_Pval", "Random_IVW_Pval")
  if (any(vapply(p_columns, function(column) {
    values <- as.numeric(data[[column]])
    any(!is.finite(values) | values < 0 | values > 1)
  }, logical(1)))) {
    stop(sprintf("%s has invalid p-values.", label), call. = FALSE)
  }
  egger_p <- as.numeric(data$Egger_Pval)
  if (any(!is.na(egger_p) & (!is.finite(egger_p) | egger_p < 0 | egger_p > 1))) {
    stop(sprintf("%s has invalid MR-Egger p-values.", label), call. = FALSE)
  }
  use_fixed <- n_ivs <= 3L
  data.frame(
    effect = ifelse(use_fixed, as.numeric(data$Fixed_IVW_Estimate), as.numeric(data$Random_IVW_Estimate)),
    ivw_p = ifelse(use_fixed, as.numeric(data$Fixed_IVW_Pval), as.numeric(data$Random_IVW_Pval)),
    stringsAsFactors = FALSE
  )
}

# Reverse-MR results
read_rev <- function(fp, lab) {
  path <- file.path(REVERSE_RESULTS, fp)
  if (!file.exists(path) || dir.exists(path)) stop(sprintf("Reverse-MR result is missing: %s", path), call. = FALSE)
  raw <- read_tsv(path, show_col_types = FALSE)
  if (!nrow(raw)) stop(sprintf("Reverse-MR result has zero rows: %s", path), call. = FALSE)
  selected <- select_reverse_ivw(raw, basename(path))
  tibble(Metabolite = raw$Metabolite, outcome = lab,
         effect = selected$effect, ivw_p = selected$ivw_p, egger_p = raw$Egger_Pval)
}
rev_long <- do.call(bind_rows, Map(read_rev, OC$file, as.character(OC$label)))
rev_long$outcome <- factor(rev_long$outcome, levels = levels(OC$label))

# Stage 08 superpathway annotation
norm_name <- function(x) gsub("[:/]", "_", x)
if (!file.exists(ANNOTATION_FILE) || dir.exists(ANNOTATION_FILE)) {
  stop(sprintf("Stage 08 annotation is missing: %s", ANNOTATION_FILE), call. = FALSE)
}
cls <- read_tsv(ANNOTATION_FILE, show_col_types = FALSE)
if (!nrow(cls) || !all(c("name", "superpathway") %in% names(cls))) {
  stop("Stage 08 annotation has an invalid schema.", call. = FALSE)
}
cls <- cls |>
  transmute(key = name, superpathway)
rev_long <- rev_long |>
  mutate(key = norm_name(Metabolite)) |>
  left_join(cls, by = "key")
stopifnot("every reverse-MR metabolite has a superpathway class" =
            all(!is.na(rev_long$superpathway)))

# Order rows by superpathway and mean absolute effect.
class_order <- c("Lipid", "Amino Acid", "Peptide", "Nucleotide", "Xenobiotics", "Unknown")
present_classes <- intersect(class_order, unique(rev_long$superpathway))
metab_order <- rev_long |>
  group_by(Metabolite, superpathway) |>
  summarise(rank_val = mean(abs(effect), na.rm = TRUE), .groups = "drop") |>
  mutate(superpathway = factor(superpathway, levels = present_classes)) |>
  arrange(superpathway, desc(rank_val))
rev_long <- rev_long |>
  mutate(Metabolite = factor(Metabolite, levels = rev(metab_order$Metabolite)),
         superpathway = factor(superpathway, levels = present_classes))

n_metab <- nlevels(rev_long$Metabolite)
message(sprintf("Reverse-MR heatmap: %d metabolites x %d outcomes; classes: %s",
                n_metab, nlevels(rev_long$outcome),
                paste(sprintf("%s=%d", present_classes,
                              tapply(metab_order$Metabolite, metab_order$superpathway, length)[present_classes]),
                      collapse = ", ")))

# Signed square-root colour scale preserves the full effect range.
FCAP <- ceiling(max(abs(rev_long$effect), na.rm = TRUE) / 0.05) * 0.05
sqrt_rescaler <- function(x, to = c(0, 1), from = c(-FCAP, FCAP)) {
  scales::rescale(sign(x) * sqrt(abs(x)), to = to, from = sign(from) * sqrt(abs(from)))
}
message(sprintf("Diverging fill: full range +/- %.2f (max |effect| = %.3f), signed-sqrt mapping, 0 = white",
                FCAP, max(abs(rev_long$effect), na.rm = TRUE)))
PUOR <- c(cvd_pal$diverging[["low"]], cvd_pal$diverging[["mid"]], cvd_pal$diverging[["high"]])

# Nominal significance flags
n_ivw   <- sum(rev_long$ivw_p   < 0.05, na.rm = TRUE)
n_egger <- sum(rev_long$egger_p < 0.05, na.rm = TRUE)
message(sprintf("Nominal reverse-causation flags: IVW p<0.05 in %d/%d cells; Egger p<0.05 in %d/%d cells",
                n_ivw, nrow(rev_long), n_egger, nrow(rev_long)))

# Categorical superpathway palette
class_pal <- c(
  "Lipid"       = "#E69F00",
  "Amino Acid"  = "#56B4E9",
  "Peptide"     = "#CC79A7",
  "Nucleotide"  = "#0072B2",
  "Xenobiotics" = "#F0E442",
  "Unknown"     = "#999999"
)[present_classes]

# Apply per-class strip tints after ggplot builds the gtable.
tint_toward_white <- function(hex, p = 0.72) {
  cc <- grDevices::col2rgb(hex) / 255
  w  <- cc * (1 - p) + p
  grDevices::rgb(w[1], w[2], w[3])
}
set_strip_fill <- function(strip_grob, fill) {
  strip_grob$grobs <- lapply(strip_grob$grobs, function(gr) {
    if (inherits(gr, "gTree") && !is.null(gr$children)) {
      gr$children <- lapply(gr$children, function(ch) {
        if (inherits(ch, "rect")) { ch$gp$fill <- fill; ch$gp$col <- NA }
        ch
      })
    } else if (inherits(gr, "rect")) { gr$gp$fill <- fill; gr$gp$col <- NA }
    gr
  })
  strip_grob
}
apply_strip_tint <- function(gg, classes, pal) {
  g  <- ggplot2::ggplotGrob(gg)
  sl <- which(grepl("^strip-l", g$layout$name))
  sl <- sl[order(g$layout$t[sl])]
  stopifnot("one left strip per class block" = length(sl) == length(classes))
  for (k in seq_along(sl))
    g$grobs[[sl[k]]] <- set_strip_fill(g$grobs[[sl[k]]], tint_toward_white(pal[[classes[k]]]))
  g
}

# Heatmap rendering
ivw_hits   <- rev_long |> filter(ivw_p   < 0.05)
egger_hits <- rev_long |> filter(egger_p < 0.05)

build_ggplot <- function() {
  # Map constant shapes so the legend matches the plotted markers.
  ggplot(rev_long, aes(x = outcome, y = Metabolite, fill = effect)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    geom_point(data = egger_hits, aes(shape = "Egger"), size = 3.0, stroke = 0.55,
               fill = NA, colour = "#111111") +
    geom_point(data = ivw_hits, aes(shape = "IVW"), size = 1.5, stroke = 0.35,
               fill = "#111111", colour = "white") +
    scale_shape_manual(
      name = NULL, breaks = c("IVW", "Egger"), values = c(IVW = 21, Egger = 21),
      labels = c(IVW = "Reverse-MR IVW p < 0.05", Egger = "MR-Egger p < 0.05"),
      guide = guide_legend(order = 1, ncol = 1, override.aes = list(
        fill = c("#111111", NA), colour = c("white", "#111111"),
        size = c(2.0, 2.7), stroke = c(0.4, 0.55)))) +
    scale_fill_gradientn(colours = PUOR, limits = c(-FCAP, FCAP), rescaler = sqrt_rescaler,
                         oob = scales::squish, breaks = c(-0.5, -0.2, 0, 0.2, 0.5),
                         na.value = "#DDDDDD",
                         name = "Reverse MR effect (IVW)\nsigned square-root colour scale",
                         guide = guide_colourbar(order = 2, theme = theme(
                           legend.title.position = "top",
                           legend.title = element_text(hjust = 0.5, lineheight = 0.9,
                                                        margin = margin(b = 3)),
                           legend.key.width  = unit(34, "mm"),
                           legend.key.height = unit(3,  "mm")))) +
    scale_x_discrete(position = "top", expand = c(0, 0),
                     labels = c("Type 2 diabetes" = "Type 2\ndiabetes",
                                "Fasting glucose" = "Fasting\nglucose",
                                "HbA1c" = "HbA1c")) +
    scale_y_discrete(expand = c(0, 0)) +
    facet_grid(rows = vars(superpathway), scales = "free_y", space = "free_y", switch = "y") +
    labs(x = NULL, y = NULL) +
    theme_cvdiab() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x.top = element_text(face = "plain", size = rel(1.0), lineheight = 0.9, vjust = 0.5),
          axis.text.y = element_text(size = rel(0.62)),
          panel.spacing.y = unit(1.1, "mm"),
          strip.placement = "outside",
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text.y.left = element_text(angle = 0, hjust = 0.5, face = "bold",
                                           size = rel(0.88), colour = cvd_pal$text,
                                           margin = margin(0, 4, 0, 4)),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.spacing.x = unit(8, "mm"),
          legend.key.spacing.y = unit(0.4, "mm"),
          legend.location = "plot",
          legend.box.margin = margin(t = 0, b = 0),
          plot.margin = margin(2, 6, 2, 2))
}

gg <- build_ggplot()
g  <- tryCatch(apply_strip_tint(gg, present_classes, class_pal),
               error = function(e) { message("strip tint skipped: ", conditionMessage(e)); gg })
save_cvd_fig(g, file.path(OUT, "fig3_reverse_heatmap.pdf"),
             width_mm = 150, height_mm = 222)
message("Wrote: ", file.path(OUT, "fig3_reverse_heatmap.pdf"))

# Reconciliation output
cat("\n=== Figure 3 reverse-MR heatmap — reconciliation ===\n")
cat(sprintf("Metabolites (rows): %d | outcomes (cols): %d | cells: %d\n",
            n_metab, nlevels(rev_long$outcome), nrow(rev_long)))
cat(sprintf("Class counts: %s\n",
            paste(sprintf("%s=%d", present_classes,
                          tapply(metab_order$Metabolite, factor(metab_order$superpathway, levels=present_classes), length)),
                  collapse = ", ")))
cat(sprintf("Diverging fill: full range +/- %.2f (signed-sqrt), max|effect| %.3f\n",
            FCAP, max(abs(rev_long$effect), na.rm = TRUE)))
cat(sprintf("Nominal reverse flags: IVW p<0.05 %d cells | Egger p<0.05 %d cells | NA (untested) %d cells\n",
            n_ivw, n_egger, sum(is.na(rev_long$effect))))
