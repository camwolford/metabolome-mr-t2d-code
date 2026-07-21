# Shared theme, direction encoding, and vector-PDF writer.
# Direction uses reserved green/red hues and redundant triangle shapes.

suppressPackageStartupMessages({
  library(ggplot2)
})

# Reserved direction palette
cvd_pal <- list(
  protective = "#1E8E4E",
  risk       = "#C0392B",
  ns         = "#BFBFBF",
  categorical = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7", "#000000"),
  sequential = "viridis",
  diverging  = c(low = "#5E3C99", mid = "#FFFFFF", high = "#E66101"),
  reference  = "#444444",
  axis       = "#333333",
  text       = "#1A1A1A"
)

# Filled triangles encode direction; circles encode non-significance.
cvd_shape <- c(protective = 25L, risk = 24L, ns = 21L)

cvd_direction_levels <- c("risk", "protective", "ns")
cvd_direction_labels <- c(risk = "Higher risk", protective = "Lower risk", ns = "Not significant")

#' theme_cvdiab
#' Shared ggplot theme.
#' @param base_size base font size in points
#' @param base_family font family
#' @return a ggplot2 theme object
theme_cvdiab <- function(base_size = 9, base_family = "Helvetica") {
  half <- base_size / 2
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line       = element_line(colour = cvd_pal$axis, linewidth = 0.4),
      axis.ticks      = element_line(colour = cvd_pal$axis, linewidth = 0.4),
      axis.text       = element_text(colour = cvd_pal$text, size = rel(0.9)),
      axis.title      = element_text(colour = cvd_pal$text, size = rel(1.0)),
      axis.title.x    = element_text(margin = margin(t = half)),
      axis.title.y    = element_text(margin = margin(r = half), angle = 90),
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.key       = element_blank(),
      legend.background = element_blank(),
      legend.title     = element_text(colour = cvd_pal$text, size = rel(0.9)),
      legend.text      = element_text(colour = cvd_pal$text, size = rel(0.85)),
      legend.position  = "bottom",
      plot.tag         = element_text(face = "bold", size = rel(1.3), colour = cvd_pal$text),
      plot.tag.position = c(0.01, 0.99),
      plot.title       = element_text(face = "plain", size = rel(1.0), hjust = 0.5,
                                      colour = cvd_pal$text,
                                      margin = margin(b = half)),
      plot.margin      = margin(half, half, half, half),
      complete = TRUE
    )
}

#' cvd_direction_scales
#' Matched colour, fill, and shape scales for effect direction.
#' @param which aesthetics to build
#' @param drop retain absent direction levels when FALSE
#' @return a list of ggplot scales to add to a plot
cvd_direction_scales <- function(which = c("colour", "fill", "shape"), drop = FALSE) {
  which <- match.arg(which, several.ok = TRUE)
  vals_col <- c(risk = cvd_pal$risk, protective = cvd_pal$protective, ns = cvd_pal$ns)
  out <- list()
  if ("colour" %in% which)
    out <- c(out, list(scale_colour_manual(values = vals_col, breaks = cvd_direction_levels,
                        labels = cvd_direction_labels, name = NULL, drop = drop)))
  if ("fill" %in% which)
    out <- c(out, list(scale_fill_manual(values = vals_col, breaks = cvd_direction_levels,
                      labels = cvd_direction_labels, name = NULL, drop = drop)))
  if ("shape" %in% which)
    out <- c(out, list(scale_shape_manual(values = cvd_shape, breaks = cvd_direction_levels,
                      labels = cvd_direction_labels, name = NULL, drop = drop)))
  out
}

#' embed_fonts_pdftocairo
#' Embed PDF fonts with pdftocairo when available.
#' @param file path to a PDF to rewrite in place
#' @return TRUE if the file was re-rendered with embedded fonts, FALSE otherwise
embed_fonts_pdftocairo <- function(file) {
  ptc <- Sys.which("pdftocairo")
  if (!nzchar(ptc)) return(FALSE)
  tmp <- paste0(file, ".embed.tmp.pdf")
  rc <- tryCatch(system2(ptc, c("-pdf", shQuote(file), shQuote(tmp)),
                         stdout = FALSE, stderr = FALSE),
                 error = function(e) 1L)
  if (isTRUE(rc == 0) && file.exists(tmp) && file.info(tmp)$size > 2000) {
    file.rename(tmp, file); TRUE
  } else { unlink(tmp); FALSE }
}

#' save_cvd_fig
#' Write a plot to vector PDF.
#' @param plot a ggplot / patchwork object
#' @param file output path (.pdf)
#' @param width_mm figure width in millimetres
#' @param height_mm figure height in millimetres
#' @return the file path (invisibly)
save_cvd_fig <- function(plot, file, width_mm = 170, height_mm = 170) {
  # Prefer Cairo and fall back to base vector PDF.
  unlink(file)
  ok <- FALSE
  if (isTRUE(capabilities("cairo"))) {
    ok <- tryCatch({
      suppressWarnings(ggsave(filename = file, plot = plot,
                              width = width_mm, height = height_mm, units = "mm",
                              device = grDevices::cairo_pdf))
      file.exists(file) && file.info(file)$size > 2000
    }, error = function(e) FALSE, warning = function(w) FALSE)
  }
  if (!ok) {
    unlink(file)
    # Preserve typographic glyphs in the fallback PDF.
    pdf_winansi <- function(filename, ...) grDevices::pdf(filename, ..., encoding = "WinAnsi.enc")
    ggsave(filename = file, plot = plot,
           width = width_mm, height = height_mm, units = "mm", device = pdf_winansi)
    embed_fonts_pdftocairo(file)
  }
  invisible(file)
}
