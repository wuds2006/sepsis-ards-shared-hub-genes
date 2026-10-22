#!/usr/bin/env Rscript
# AJRCCM figure theme + ggsci NPG palette + Illustrator-editable PDF export
# Shared by all ggplot2 figure-redraw scripts (step30_*.R)
# Conventions:
#   - Vector PDF via cairo_pdf: text stays as text (no outlines), fonts embedded
#   - useDingbats = FALSE for Illustrator compatibility
#   - AJRCCM widths: single column 89 mm, 1.5 column 140 mm, double column 178 mm
#   - Base font Arial 7 pt (final-size rule: >= 6 pt at print size)
#   - Global categorical palette: ggsci NPG (npg)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
  library(patchwork)
})

AJRCCM_W_SINGLE <- 89
AJRCCM_W_ONEHALF <- 140
AJRCCM_W_DOUBLE <- 178
AJRCCM_FONT <- "Helvetica"

# NPG palette (10 colors, recycled if more needed)
pal_npg <- function(n) {
  cols <- ggsci::pal_npg("nrc")(10)
  if (n <= 10) cols[seq_len(n)] else colorRampPalette(cols)(n)
}

# Disease-vs-control two-group colors: disease = NPG red, control = NPG blue
COL_DISEASE <- "#E64B35FF"
COL_CONTROL <- "#4DBBD5FF"
COL_UP      <- "#E64B35FF"   # upregulated
COL_DOWN    <- "#3C5488FF"   # downregulated
COL_NS      <- "grey75"

theme_ajrccm <- function(base_size = 7, base_family = AJRCCM_FONT) {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid       = element_blank(),
      panel.border     = element_rect(color = "black", linewidth = 0.4),
      axis.line        = element_blank(),
      axis.ticks       = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(1.2, "mm"),
      axis.text        = element_text(color = "black", size = base_size),
      axis.title       = element_text(color = "black", size = base_size),
      plot.title       = element_text(size = base_size, face = "bold", hjust = 0),
      plot.subtitle    = element_text(size = base_size * 0.9),
      legend.text      = element_text(size = base_size * 0.9),
      legend.title     = element_text(size = base_size * 0.9),
      legend.key.size  = unit(3, "mm"),
      legend.margin    = margin(0, 0, 0, 0),
      plot.margin      = margin(2, 2, 2, 2),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(size = base_size, face = "bold", margin = margin(1, 1, 1, 1))
    )
}

# Save ggplot object as Illustrator-editable vector PDF (width/height in mm)
# Uses the base pdf() device with Helvetica (PDF core font): text stays as text,
# never outlined, and Illustrator maps Helvetica to the local font on open.
# (cairo_pdf is avoided: CID-font lookup for Arial fails on this system.)
save_pdf <- function(plot, filename, width_mm, height_mm, ...) {
  grDevices::pdf(filename, width = width_mm / 25.4, height = height_mm / 25.4,
                 family = "Helvetica", useDingbats = FALSE, ...)
  print(plot)
  dev.off()
  invisible(filename)
}

# Significance stars from two-sided Wilcoxon P value
sig_stars <- function(p) {
  ifelse(p < 1e-4, "****", ifelse(p < 1e-3, "***",
         ifelse(p < 1e-2, "**", ifelse(p < 0.05, "*", "ns"))))
}

# Italic gene labels helper: returns parse-able expression vector
gene_italic <- function(g) paste0("italic('", g, "')")

# Project paths
PROJ <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT_MAIN <- file.path(PROJ, "AJRCCM_submission/03_figures/Figures_ggplot2_NPG/main_panels")
OUT_SUPP <- file.path(PROJ, "AJRCCM_submission/03_figures/Figures_ggplot2_NPG/supp_panels")
OUT_COMP <- file.path(PROJ, "AJRCCM_submission/03_figures/Figures_ggplot2_NPG/composites")
dir.create(OUT_MAIN, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_SUPP, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_COMP, recursive = TRUE, showWarnings = FALSE)

HUB6 <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")

cat("ajrccm_theme.R loaded. Outputs:\n ", OUT_MAIN, "\n ", OUT_SUPP, "\n ", OUT_COMP, "\n")
