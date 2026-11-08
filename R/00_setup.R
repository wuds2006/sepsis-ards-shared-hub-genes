# ---------------------------------------------------------------------------
# 00_setup.R -- shared graphics constants and export helpers
#
# All figure scripts source this file. It fixes the physical layout
# (180 mm column), the type specification (Helvetica; 7 pt body, 9 pt bold
# panel tags), the qualitative palette, and the four output formats
# (PNG 300 dpi, PDF vector, LZW TIFF 300 dpi, SVG with live text nodes).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(gridExtra)
  library(scales)
})

# --- geometry --------------------------------------------------------------

MM_PER_IN <- 25.4
FIG_W_MM   <- 180                      # journal column width
FIG_H_MIN  <- 200                      # preferred height window
FIG_H_MAX  <- 220

mm <- function(x) x / MM_PER_IN        # mm -> inches, for ggsave-style devices

# --- type ------------------------------------------------------------------

FONT_BODY  <- "Helvetica"
FONT_TAG   <- "Helvetica"
PT_BODY    <- 7
PT_TAG     <- 9
PT_SMALL   <- 6                        # only for axis tick annotations that
                                       # cannot be shortened; used sparingly

# --- palette ---------------------------------------------------------------
# Low-chroma land/sea palette: no pure primaries, all series separated in both
# hue and lightness so that they survive greyscale conversion.

pal <- c(
  brick   = "#B05B5B",   # disease / ARDS
  steel   = "#4C6E91",   # comparison group / sepsis
  sage    = "#7E9E95",   # control / healthy
  ochre   = "#C99A5B",   # third category, highlight
  plum    = "#8A7391",   # fourth category
  olive   = "#7D8C5C",   # fifth category
  slate   = "#6F7B85",   # neutral reference line, annotation
  sand    = "#D9C9B4",   # fill for non-significant features
  ink     = "#2B2B2B"    # text and rules
)

# six hub genes: fixed order and fixed colour so every panel agrees
HUB6 <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")

hub_col <- c(
  CDK5RAP2 = "#4C6E91",
  CD177    = "#B05B5B",
  GRB10    = "#7E9E95",
  UPP1     = "#C99A5B",
  EXOSC4   = "#8A7391",
  VNN1     = "#7D8C5C"
)

# diverging ramp for scaled expression / correlation (blue - white - red)
ramp_div <- c("#3E5C7E", "#7E9CB8", "#DCE3EA", "#F5F3EF", "#E2C6C2",
              "#C08080", "#A04848")

# sequential ramp for activity scores and -log10 P
ramp_seq <- c("#F5F3EF", "#DCE3EA", "#A9BFD0", "#7E9CB8", "#56789E", "#3E5C7E")

# group colours used across bulk, single-cell and pharmacology panels
grp_col <- c(Control = "#7E9E95", Healthy = "#7E9E95", HC = "#7E9E95",
             ARDS = "#B05B5B", Sepsis = "#4C6E91", Disease = "#B05B5B")

# --- theme -----------------------------------------------------------------

theme_ajr <- function(base_size = PT_BODY, base_family = FONT_BODY) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text              = element_text(colour = pal[["ink"]]),
      axis.line         = element_line(linewidth = 0.35, colour = pal[["ink"]]),
      axis.ticks        = element_line(linewidth = 0.35, colour = pal[["ink"]]),
      axis.ticks.length = unit(1.1, "pt"),
      axis.text         = element_text(size = base_size, colour = pal[["ink"]]),
      axis.title        = element_text(size = base_size, colour = pal[["ink"]]),
      plot.title        = element_text(size = base_size + 1, face = "bold",
                                       hjust = 0, colour = pal[["ink"]]),
      plot.subtitle     = element_text(size = base_size, hjust = 0,
                                       colour = pal[["slate"]]),
      plot.caption      = element_text(size = base_size - 1, hjust = 0,
                                       colour = pal[["slate"]]),
      strip.background  = element_blank(),
      strip.text        = element_text(size = base_size, face = "bold",
                                       colour = pal[["ink"]]),
      legend.title      = element_text(size = base_size, colour = pal[["ink"]]),
      legend.text       = element_text(size = base_size, colour = pal[["ink"]]),
      legend.key.size   = unit(7, "pt"),
      legend.key.height = unit(7, "pt"),
      legend.key.width  = unit(9, "pt"),
      legend.background = element_blank(),
      legend.margin     = margin(0, 0, 0, 0),
      panel.grid        = element_blank(),
      plot.margin       = margin(2, 3, 2, 3)
    )
}

# theme_void inherits the 0.8 x rel() default for legend text from theme_grey,
# which would drop legend labels to 5.6 pt; this variant pins them to 7 pt.
theme_void_ajr <- function(base_size = PT_BODY, base_family = FONT_BODY) {
  theme_void(base_size = base_size, base_family = base_family) +
    theme(text = element_text(colour = pal[["ink"]], family = base_family),
          legend.text = element_text(size = base_size, colour = pal[["ink"]]),
          legend.title = element_text(size = base_size, colour = pal[["ink"]]),
          plot.title = element_text(size = base_size, face = "bold", hjust = 0),
          plot.subtitle = element_text(size = base_size, hjust = 0,
                                       colour = pal[["slate"]]))
}

theme_set(theme_ajr())

# heatmaps are drawn with ComplexHeatmap; keep one parameter block for them
hm_par <- list(
  fontsize         = PT_BODY,
  fontsize_row     = PT_BODY - 1,
  fontsize_col     = PT_BODY - 1,
  fontsize_legend  = PT_SMALL,
  border           = TRUE,
  rect_gp          = gpar(col = "white", lwd = 0.4),
  row_names_gp     = gpar(fontfamily = FONT_BODY, fontsize = PT_BODY - 1),
  column_names_gp  = gpar(fontfamily = FONT_BODY, fontsize = PT_BODY - 1),
  heatmap_legend_param = list(
    labels_gp  = gpar(fontfamily = FONT_BODY, fontsize = PT_SMALL),
    title_gp   = gpar(fontfamily = FONT_BODY, fontsize = PT_SMALL),
    grid_border = "transparent"
  )
)

# --- panel tags ------------------------------------------------------------
# Tags sit outside the panel, top-left, and are the only bold 9 pt element on
# the page. tag_plot() is for single panels, tag_fig() for an assembled figure.

tag_plot <- function(p, lab, x = 0.004, y = 1.03) {
  p + labs(tag = lab) +
    theme(plot.tag = element_text(family = FONT_TAG, size = PT_TAG,
                                  face = "bold", colour = pal[["ink"]]),
          plot.tag.position = c(x, y))
}

tag_fig <- function(p, labels = NULL, x = 0.008, y = 0.995) {
  if (is.null(labels)) return(p)
  n <- length(labels)
  xs <- if (length(x) == 1) rep(x, n) else x
  ys <- if (length(y) == 1) rep(y, n) else y
  for (i in seq_len(n)) {
    p <- p + cowplot::draw_label(labels[i], x = xs[i], y = ys[i], hjust = 0, vjust = 1,
                        fontfamily = FONT_TAG, size = PT_TAG, fontface = "bold",
                        colour = pal[["ink"]])
  }
  p
}

# --- export ----------------------------------------------------------------

out_dir <- function(kind = c("main", "supplementary", "render")) {
  kind <- match.arg(kind)
  root <- Sys.getenv("FIG_ROOT", unset = file.path("..", "Figures"))
  d <- if (kind == "render") file.path("..", "_work", "render") else file.path(root, kind)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# Draw `p` once per device so that every format carries identical geometry.
# `p` may be a ggplot or a grid grob / gTree (ComplexHeatmap draw output).
save_fig <- function(p, name, width_mm = FIG_W_MM, height_mm,
                     dir = out_dir("main"), formats = c("pdf", "svg", "png", "tiff")) {
  stopifnot(length(height_mm) == 1)
  w <- mm(width_mm); h <- mm(height_mm)

  draw <- function() {
    if (inherits(p, "ggplot") || inherits(p, "patchwork")) {
      print(p)
    } else {
      grid.newpage(); grid.draw(p)
    }
  }

  if ("pdf" %in% formats) {
    # grDevices::pdf keeps text as text and honours semi-transparent fills;
    # cairo_pdf is unavailable on this machine, and the SVG/PNG devices below
    # reproduce the same geometry independently.
    grDevices::pdf(file.path(dir, paste0(name, ".pdf")),
                   width = w, height = h, family = FONT_BODY,
                   useDingbats = FALSE, onefile = FALSE, bg = "white",
                   colormodel = "srgb")
    draw(); grDevices::dev.off()
  }
  if ("svg" %in% formats) {
    svglite::svglite(file.path(dir, paste0(name, ".svg")), width = w, height = h,
                     bg = "white")
    draw(); grDevices::dev.off()
  }
  png_file <- file.path(dir, paste0(name, ".png"))
  if (any(c("png", "tiff") %in% formats)) {
    ragg::agg_png(png_file, width = w, height = h, units = "in", res = 300,
                  background = "white")
    draw(); grDevices::dev.off()
    # ragg does not always write a pHYs chunk; stamp the resolution so that
    # downstream converters and the journal's checker read 300 dpi
    img <- magick::image_read(png_file)
    magick::image_write(img, png_file, format = "png", density = "300x300")
  }
  if ("tiff" %in% formats) {
    img <- magick::image_read(png_file)
    magick::image_write(img, file.path(dir, paste0(name, ".tiff")),
                        format = "tiff", compression = "lzw", density = "300x300")
  }
  invisible(file.path(dir, paste0(name, c(".pdf", ".png"))))
}

# --- small utilities -------------------------------------------------------
# Helvetica on this platform carries no subscript or superscript digits and no
# Greek letters, so axis labels are written in plain text ("log2 fold change",
# "-log10 adjusted P") rather than with plotmath, which would pull in the
# Symbol font and shrink glyphs to 4.9 pt.

fmt_p <- function(p, digits = 2) {
  vapply(p, function(x) {
    if (is.na(x)) return("NA")
    if (x < 1e-4) formatC(x, format = "e", digits = digits)
    else formatC(x, format = "f", digits = if (x < 0.01) 4 else 3)
  }, character(1))
}

star_p <- function(p) {
  cut(p, c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf),
      labels = c("****", "***", "**", "*", "ns"), right = FALSE)
}

read_tab <- function(...) {
  utils::read.csv(file.path("data", ...), stringsAsFactors = FALSE,
                  check.names = FALSE)
}

msg <- function(...) cat(sprintf(...), "\n", sep = "")
