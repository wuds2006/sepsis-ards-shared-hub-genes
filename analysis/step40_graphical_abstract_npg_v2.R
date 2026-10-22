# ============================================================
# Graphical Abstract v2.1 (AJRCCM) - ggplot2 + ggsci NPG
# Style: multi-panel data-rich layout after reference GA
#   left panel  : Identification of potential core hub genes
#   right-top   : Exploration of molecular mechanisms
#   right-bottom: Drug prediction
#   center      : vertical ellipse - six shared hub genes
# Font: Helvetica | Palette: ggsci NPG | >=300 dpi, no title inside
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggforce)
  library(ggsci)
})

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "AJRCCM_submission/04_graphical_abstract")

npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
         "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")
RED <- npg[1]; CYAN <- npg[2]; GREEN <- npg[3]; BLUE <- npg[4]
CORAL <- npg[5]; GREY <- npg[6]; LTEAL <- npg[7]
TEXT <- "grey10"
FONT <- "Helvetica"

tint <- function(col, p = 0.90) {
  rgb <- col2rgb(col)[, 1]
  rgb(255 - (255 - rgb[1]) * (1 - p),
      255 - (255 - rgb[2]) * (1 - p),
      255 - (255 - rgb[3]) * (1 - p),
      maxColorValue = 255)
}

arc_pts <- function(cx, cy, r, from_deg, to_deg, n = 8) {
  th <- seq(from_deg, to_deg, length.out = n) * pi / 180
  data.frame(x = cx + r * cos(th), y = cy + r * sin(th))
}
rrpoly <- function(id, x, y, w, h, r) {
  data.frame(id = id, rbind(
    data.frame(x = x + r, y = y),
    arc_pts(x + w - r, y + r, r, 270, 360),
    data.frame(x = x + w, y = y + h - r),
    arc_pts(x + w - r, y + h - r, r, 0, 90),
    data.frame(x = x + r, y = y + h),
    arc_pts(x + r, y + h - r, r, 90, 180),
    data.frame(x = x, y = y + r),
    arc_pts(x + r, y + r, r, 180, 270)))
}

W <- 100; H <- 70
arrow_spec <- arrow(length = unit(2.4, "mm"), type = "closed")

p <- ggplot() +
  coord_fixed(xlim = c(0, W), ylim = c(0, H), expand = FALSE) +
  theme_void(base_family = FONT) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(0, 0, 0, 0))

# ============================================================
# panels
# ============================================================
pl <- rrpoly("pl", 1, 1.5, 44.5, 67, 1.2)
rt <- rrpoly("rt", 52, 37, 47, 31.5, 1.2)
rb <- rrpoly("rb", 52, 1.5, 47, 31.5, 1.2)
pans <- rbind(pl, rt, rb)
pans$fill <- NA
pans$fill[pans$id == "pl"] <- tint(GREEN, 0.90)
pans$fill[pans$id == "rt"] <- tint(CORAL, 0.86)
pans$fill[pans$id == "rb"] <- tint(BLUE, 0.90)
pans$edge <- NA
pans$edge[pans$id == "pl"] <- GREEN
pans$edge[pans$id == "rt"] <- CORAL
pans$edge[pans$id == "rb"] <- BLUE

p <- p +
  geom_polygon(data = pans, aes(x = x, y = y, group = id),
               fill = pans$fill, color = pans$edge, linewidth = 0.8)

p <- p +
  annotate("text", x = 23.25, y = 66.2, label = "Identification of potential core hub genes",
           family = FONT, fontface = "bold", size = 4.5, color = GREEN) +
  annotate("text", x = 75.5, y = 66.4, label = "Exploration of molecular mechanisms",
           family = FONT, fontface = "bold", size = 4.5, color = RED) +
  annotate("text", x = 75.5, y = 30.7, label = "Drug prediction",
           family = FONT, fontface = "bold", size = 4.5, color = BLUE)

# ============================================================
# LEFT PANEL - identification
# ============================================================
p <- p +
  geom_circle(data = data.frame(x0 = 5.4, y0 = 61.6, r = 1.5),
              aes(x0 = x0, y0 = y0, r = r), fill = RED, color = NA) +
  annotate("text", x = 5.4, y = 61.6, label = "L", size = 2.6,
           color = "white", fontface = "bold", family = FONT) +
  annotate("text", x = 8.2, y = 62.9, label = "GSE76293", hjust = 0,
           family = FONT, fontface = "bold.italic", size = 3.3, color = TEXT) +
  annotate("text", x = 8.2, y = 60.9, label = "ARDS \u00b7 PMNs (n = 12 vs 12)", hjust = 0,
           family = FONT, size = 2.8, color = TEXT) +
  annotate("text", x = 4.2, y = 58.2, label = "557 DEGs", hjust = 0,
           family = FONT, size = 2.8, color = TEXT) +
  annotate("text", x = 4.2, y = 56.6, label = "(291 up / 266 down)", hjust = 0,
           family = FONT, size = 2.5, color = "grey30") +
  annotate("text", x = 4.2, y = 55.0, label = "13 WGCNA modules", hjust = 0,
           family = FONT, size = 2.5, color = "grey30") +
  geom_circle(data = data.frame(x0 = 22.4, y0 = 61.6, r = 1.5),
              aes(x0 = x0, y0 = y0, r = r), fill = CYAN, color = NA) +
  annotate("text", x = 22.4, y = 61.6, label = "B", size = 2.6,
           color = "white", fontface = "bold", family = FONT) +
  annotate("text", x = 25.2, y = 62.9, label = "GSE69528", hjust = 0,
           family = FONT, fontface = "bold.italic", size = 3.3, color = TEXT) +
  annotate("text", x = 25.2, y = 60.9, label = "Sepsis \u00b7 whole blood (n = 83 vs 28)", hjust = 0,
           family = FONT, size = 2.8, color = TEXT) +
  annotate("text", x = 21.2, y = 58.2, label = "1,458 DEGs", hjust = 0,
           family = FONT, size = 2.8, color = TEXT) +
  annotate("text", x = 21.2, y = 56.6, label = "(658 up / 800 down)", hjust = 0,
           family = FONT, size = 2.5, color = "grey30") +
  annotate("text", x = 21.2, y = 55.0, label = "11 WGCNA modules", hjust = 0,
           family = FONT, size = 2.5, color = "grey30") +
  # Venn
  geom_circle(data = data.frame(x0 = 39.6, y0 = 61.0, r = 2.4),
              aes(x0 = x0, y0 = y0, r = r),
              fill = alpha(RED, 0.45), color = RED, linewidth = 0.5) +
  geom_circle(data = data.frame(x0 = 42.6, y0 = 61.0, r = 2.4),
              aes(x0 = x0, y0 = y0, r = r),
              fill = alpha(CYAN, 0.45), color = CYAN, linewidth = 0.5) +
  annotate("text", x = 41.1, y = 56.9, label = "DEGs & WGCNA",
           family = FONT, size = 2.8, color = TEXT, fontface = "bold")

# --- 179 shared candidate genes -> ML + SHAP ---
p <- p +
  annotate("segment", x = 23.25, y = 53.4, xend = 23.25, yend = 51.4,
           color = GREY, linewidth = 0.9, arrow = arrow_spec) +
  annotate("text", x = 23.25, y = 50.0, label = "179 shared candidate genes",
           family = FONT, fontface = "bold", size = 3.7, color = GREEN) +
  annotate("text", x = 23.25, y = 48.2, label = "\u2229 DEGs & disease-associated modules, both cohorts",
           family = FONT, size = 2.6, color = "grey30") +
  annotate("text", x = 23.25, y = 46.0,
           label = "11 machine-learning classifiers + SHAP (AUC 0.75\u20131.00)",
           family = FONT, fontface = "bold", size = 3.2, color = TEXT)

# --- model selection bar chart ---
ms_x <- c(8.0, 12.3, 16.6, 20.9)
cls <- c("LASSO", "RF", "XGBoost", "SVM")
acc <- c(0.82, 0.88, 0.95, 0.90); f1 <- c(0.80, 0.86, 0.94, 0.88); auc <- c(0.86, 0.92, 1.00, 0.93)
base_ms <- 32.5; top_ms <- 41.5; bw <- 1.05
yb <- function(v) base_ms + (v - 0.6) / 0.4 * (top_ms - base_ms)
err <- c(0.04, 0.03, 0.02, 0.03)
p <- p +
  annotate("text", x = 14.0, y = 44.0, label = "Model selection",
           family = FONT, fontface = "bold", size = 3.2, color = TEXT) +
  annotate("segment", x = 5.8, y = base_ms, xend = 23.9, yend = base_ms,
           color = "grey35", linewidth = 0.5)
for (i in seq_along(ms_x)) {
  p <- p +
    annotate("rect", xmin = ms_x[i] - 1.5*bw, xmax = ms_x[i] - 0.5*bw,
             ymin = base_ms, ymax = yb(acc[i]), fill = CYAN, color = NA) +
    annotate("rect", xmin = ms_x[i] - 0.5*bw, xmax = ms_x[i] + 0.5*bw,
             ymin = base_ms, ymax = yb(f1[i]), fill = GREY, color = NA) +
    annotate("rect", xmin = ms_x[i] + 0.5*bw, xmax = ms_x[i] + 1.5*bw,
             ymin = base_ms, ymax = yb(auc[i]), fill = BLUE, color = NA) +
    annotate("errorbar", x = ms_x[i] + 1.0*bw, ymin = yb(auc[i] - err[i]),
             ymax = yb(pmin(1, auc[i] + err[i])), width = 0.4,
             color = "grey25", linewidth = 0.4) +
    annotate("text", x = ms_x[i], y = base_ms - 1.2, label = cls[i],
             family = FONT, size = 2.4, color = "grey25")
}
p <- p +
  annotate("rect", xmin = 24.6, xmax = 25.7, ymin = 39.6, ymax = 40.7, fill = CYAN, color = NA) +
  annotate("text", x = 26.3, y = 40.2, label = "Accuracy", hjust = 0, family = FONT, size = 2.3, color = "grey25") +
  annotate("rect", xmin = 24.6, xmax = 25.7, ymin = 37.8, ymax = 38.9, fill = GREY, color = NA) +
  annotate("text", x = 26.3, y = 38.4, label = "F1 score", hjust = 0, family = FONT, size = 2.3, color = "grey25") +
  annotate("rect", xmin = 24.6, xmax = 25.7, ymin = 36.0, ymax = 37.1, fill = BLUE, color = NA) +
  annotate("text", x = 26.3, y = 36.6, label = "AUC", hjust = 0, family = FONT, size = 2.3, color = "grey25")

# --- feature selection stacked bars (SHAP-style) ---
fe_genes <- c("VNN1", "CD177", "GRB10", "UPP1", "EXOSC4", "CDK5RAP2")
fe_len <- list(c(4.2, 3.0, 2.4), c(3.8, 3.2, 2.0), c(3.4, 2.8, 2.2, 1.4),
               c(3.0, 2.6, 1.8), c(2.8, 2.2, 1.6), c(2.4, 2.0, 1.5))
fe_cols <- list(c(BLUE, CYAN, GREEN), c(RED, BLUE, GREEN), c(GREY, CORAL, CYAN, BLUE),
                c(GREEN, BLUE, RED), c(CYAN, GREY, GREEN), c(CORAL, BLUE, GREY))
p <- p +
  annotate("text", x = 36.2, y = 44.0, label = "Feature selection (SHAP)",
           family = FONT, fontface = "bold", size = 3.2, color = TEXT)
fy <- 41.6
for (g in seq_along(fe_genes)) {
  seg_x <- 31.4
  for (s in seq_along(fe_len[[g]])) {
    p <- p + annotate("rect",
      xmin = seg_x, xmax = seg_x + fe_len[[g]][s],
      ymin = fy - 0.72, ymax = fy + 0.72,
      fill = fe_cols[[g]][s], color = "white", linewidth = 0.3)
    seg_x <- seg_x + fe_len[[g]][s] + 0.22
  }
  p <- p + annotate("text", x = 30.8, y = fy, label = fe_genes[g], hjust = 1,
                    family = FONT, fontface = "bold.italic", size = 2.4, color = TEXT)
  fy <- fy - 1.85
}
p <- p +
  annotate("segment", x = 31.2, y = 31.0, xend = 44.4, yend = 31.0,
           color = "grey35", linewidth = 0.5) +
  annotate("text", x = 37.8, y = 30.0, label = "cumulative importance", family = FONT,
           size = 2.3, color = "grey30")

# --- external validation horizontal bars ---
p <- p +
  annotate("text", x = 4.6, y = 27.0, label = "External validation (six-gene model)", hjust = 0,
           family = FONT, fontface = "bold", size = 3.2, color = TEXT) +
  annotate("rect", xmin = 13.5, xmax = 13.5 + 0.984 * 24, ymin = 23.8, ymax = 25.6,
           fill = GREEN, color = NA) +
  annotate("text", x = 12.9, y = 24.7, label = "ARDS", hjust = 1, family = FONT, size = 2.8, color = "grey25") +
  annotate("text", x = 13.5 + 0.984 * 24 + 0.7, y = 24.7, label = "AUC 0.984", hjust = 0,
           family = FONT, fontface = "bold", size = 2.8, color = GREEN) +
  annotate("rect", xmin = 13.5, xmax = 13.5 + 0.708 * 24, ymin = 21.4, ymax = 23.2,
           fill = GREY, color = NA) +
  annotate("text", x = 12.9, y = 22.3, label = "sepsis", hjust = 1, family = FONT, size = 2.8, color = "grey25") +
  annotate("text", x = 13.5 + 0.708 * 24 + 0.7, y = 22.3, label = "AUC 0.708", hjust = 0,
           family = FONT, size = 2.8, color = "grey25") +
  annotate("text", x = 4.6, y = 19.4, hjust = 0,
           label = "single-gene AUC 0.753\u20130.998 (GSE243066, GSE100159)",
           family = FONT, size = 2.7, color = "grey30")

# --- hub-gene expression heat strip (six genes x two cohorts) ---
fe_genes_h <- c("VNN1", "CD177", "GRB10", "UPP1", "EXOSC4", "CDK5RAP2")
z_ards <- c(0.9, 0.6, -0.5, 1.0, -0.3, 0.7)
z_seps <- c(0.7, 1.0, -0.7, 0.9, -0.6, 0.5)
hx <- 12.5; hy0 <- 15.6
p <- p +
  annotate("text", x = 4.6, y = 18.4, label = "Hub-gene expression (scaled)", hjust = 0,
           family = FONT, fontface = "bold", size = 3.0, color = TEXT) +
  annotate("text", x = hx + 2.6, y = hy0 + 1.7, label = "ARDS", family = FONT, size = 2.4, color = "grey25") +
  annotate("text", x = hx + 8.0, y = hy0 + 1.7, label = "sepsis", family = FONT, size = 2.4, color = "grey25")
for (g in seq_along(fe_genes_h)) {
  yrow <- hy0 - (g - 1) * 1.7
  za <- z_ards[g]; zs <- z_seps[g]
  ca <- if (za > 0) alpha(RED, 0.35 + 0.5 * abs(za)) else alpha(BLUE, 0.35 + 0.5 * abs(za))
  cs <- if (zs > 0) alpha(RED, 0.35 + 0.5 * abs(zs)) else alpha(BLUE, 0.35 + 0.5 * abs(zs))
  p <- p +
    annotate("text", x = hx - 1.4, y = yrow, label = fe_genes_h[g], hjust = 1,
             family = FONT, fontface = "bold.italic", size = 2.4, color = TEXT) +
    annotate("rect", xmin = hx, xmax = hx + 5.2, ymin = yrow - 0.6, ymax = yrow + 0.6,
             fill = ca, color = "white", linewidth = 0.25) +
    annotate("rect", xmin = hx + 5.4, xmax = hx + 10.6, ymin = yrow - 0.6, ymax = yrow + 0.6,
             fill = cs, color = "white", linewidth = 0.25)
}
p <- p +
  annotate("text", x = 26.5, y = 10.0, hjust = 0,
           label = "CD177 and UPP1 up-regulated; GRB10 and\nEXOSC4 down-regulated in both diseases",
           family = FONT, size = 2.7, color = "grey30", lineheight = 0.95)

# ============================================================
# CENTER ELLIPSE - six shared hub genes
# ============================================================
p <- p +
  geom_ellipse(data = data.frame(x0 = 48.8, y0 = 35, a = 7.3, b = 11.6, angle = 0),
               aes(x0 = x0, y0 = y0, a = a, b = b, angle = angle),
               fill = tint(CORAL, 0.55), color = CORAL, linewidth = 1.1) +
  annotate("text", x = 48.8, y = 40.6, label = "Six shared",
           family = FONT, fontface = "bold", size = 3.8, color = TEXT) +
  annotate("text", x = 48.8, y = 38.6, label = "hub genes",
           family = FONT, fontface = "bold", size = 3.8, color = TEXT) +
  annotate("text", x = 48.8, y = 35.9, label = "CDK5RAP2   CD177",
           family = FONT, fontface = "bold.italic", size = 3.0, color = TEXT) +
  annotate("text", x = 48.8, y = 34.0, label = "GRB10   UPP1",
           family = FONT, fontface = "bold.italic", size = 3.0, color = TEXT) +
  annotate("text", x = 48.8, y = 32.1, label = "EXOSC4   VNN1",
           family = FONT, fontface = "bold.italic", size = 3.0, color = TEXT) +
  annotate("text", x = 48.8, y = 29.0, label = "markers of disease\nstate",
           family = FONT, size = 2.6, color = "grey30", lineheight = 0.95)

p <- p +
  annotate("segment", x = 56.4, y = 39.5, xend = 59.6, yend = 43.5,
           color = GREY, linewidth = 1.2, arrow = arrow_spec) +
  annotate("segment", x = 56.4, y = 30.5, xend = 59.6, yend = 26.5,
           color = GREY, linewidth = 1.2, arrow = arrow_spec)

# ============================================================
# RIGHT TOP - molecular mechanisms
# ============================================================
p <- p +
  annotate("text", x = 54.2, y = 63.8, label = "Immune infiltration & single-cell", hjust = 0,
           family = FONT, fontface = "bold", size = 3.2, color = RED)
fp_x0 <- 60.0
fp_rows <- c("M2 macrophages", "CD8\u207a T cells", "Monocytes")
fp_y <- c(60.4, 58.5, 56.6)
fp_sep <- c(0.32, -0.28, 0.22); fp_ard <- c(0.12, -0.05, 0.18)
ciw <- 0.09
for (i in seq_along(fp_rows)) {
  p <- p +
    annotate("text", x = fp_x0 - 1.2, y = fp_y[i], label = fp_rows[i], hjust = 1,
             family = FONT, size = 2.5, color = "grey20") +
    annotate("segment", x = fp_x0 + fp_sep[i] - ciw, xend = fp_x0 + fp_sep[i] + ciw,
             y = fp_y[i] + 0.30, yend = fp_y[i] + 0.30, color = RED, linewidth = 0.7) +
    annotate("point", x = fp_x0 + fp_sep[i], y = fp_y[i] + 0.30, size = 1.8, color = RED) +
    annotate("segment", x = fp_x0 + fp_ard[i] - ciw, xend = fp_x0 + fp_ard[i] + ciw,
             y = fp_y[i] - 0.30, yend = fp_y[i] - 0.30, color = GREEN, linewidth = 0.7) +
    annotate("point", x = fp_x0 + fp_ard[i], y = fp_y[i] - 0.30, size = 1.8, color = GREEN)
}
p <- p +
  annotate("segment", x = fp_x0, y = 55.8, xend = fp_x0, yend = 61.6,
           color = "grey40", linewidth = 0.4, linetype = "dashed") +
  annotate("segment", x = fp_x0 - 1.0, y = 55.4, xend = fp_x0 + 1.0, yend = 55.4,
           color = "grey35", linewidth = 0.5) +
  annotate("text", x = fp_x0, y = 54.4, label = "ssGSEA score (vs healthy)",
           family = FONT, size = 2.3, color = "grey30") +
  annotate("point", x = 62.6, y = 61.6, size = 1.8, color = RED) +
  annotate("text", x = 62.9, y = 61.6, label = "sepsis", hjust = 0,
           family = FONT, size = 2.3, color = "grey25") +
  annotate("point", x = 62.6, y = 60.3, size = 1.8, color = GREEN) +
  annotate("text", x = 62.9, y = 60.3, label = "ARDS", hjust = 0,
           family = FONT, size = 2.3, color = "grey25")

# --- single-gene GSEA curves ---
p <- p +
  annotate("text", x = 76.5, y = 63.8, label = "Single-gene GSEA", hjust = 0,
           family = FONT, fontface = "bold", size = 3.2, color = RED)
gx <- seq(0, 1, 0.02)
g1 <- gx * exp(-2.4 * gx); g1 <- g1 / max(g1)
g2 <- gx * exp(-3.2 * gx); g2 <- g2 / max(g2)
p <- p +
  annotate("text", x = 76.5, y = 61.6, parse = TRUE, hjust = 0,
           label = "UPP1*','~P~'='~5.8%*%10^-14~'(sepsis)'",
           family = FONT, size = 2.6, color = "grey25") +
  annotate("segment", x = 68.5, y = 55.4, xend = 84.5, yend = 55.4, color = "grey35", linewidth = 0.5) +
  annotate("path", x = 68.5 + gx * 16, y = 55.6 + g1 * 5.0, color = RED, linewidth = 0.9) +
  annotate("path", x = 68.5 + gx * 16, y = 55.6 + g2 * 5.0, color = GREEN, linewidth = 0.9) +
  annotate("text", x = 76.5, y = 53.4, label = "UPP1 enriched in both diseases",
           family = FONT, size = 2.4, color = "grey30") +
  annotate("text", x = 86.0, y = 59.8, label = "sepsis", color = RED, family = FONT, size = 2.4) +
  annotate("text", x = 86.0, y = 58.4, label = "ARDS", color = GREEN, family = FONT, size = 2.4)

# --- transcription & translation: CEBPB -> UPP1 + miRNAs ---
p <- p +
  annotate("text", x = 75.5, y = 50.2, label = "Transcription & translation", hjust = 0.5,
           family = FONT, fontface = "bold", size = 3.2, color = RED) +
  geom_circle(data = data.frame(x0 = 63.6, y0 = 42.6, r = 1.4),
              aes(x0 = x0, y0 = y0, r = r), fill = GREY, color = NA) +
  geom_circle(data = data.frame(x0 = 65.6, y0 = 41.2, r = 1.4),
              aes(x0 = x0, y0 = y0, r = r), fill = LTEAL, color = NA) +
  geom_circle(data = data.frame(x0 = 67.7, y0 = 42.7, r = 1.4),
              aes(x0 = x0, y0 = y0, r = r), fill = CORAL, color = NA) +
  annotate("segment", x = 64.6, y = 43.6, xend = 68.8, yend = 45.0, color = "grey45", linewidth = 0.4) +
  annotate("segment", x = 67.0, y = 42.5, xend = 69.3, yend = 44.4, color = "grey45", linewidth = 0.4) +
  annotate("text", x = 64.8, y = 44.7, label = "miRNAs", family = FONT, size = 2.3, color = "grey30") +
  geom_circle(data = data.frame(x0 = 71.0, y0 = 45.4, r = 2.5),
              aes(x0 = x0, y0 = y0, r = r), fill = CORAL, color = NA) +
  annotate("text", x = 71.0, y = 45.4, label = "CEBPB", family = FONT,
           fontface = "bold", size = 2.8, color = "white") +
  geom_circle(data = data.frame(x0 = 81.6, y0 = 45.4, r = 2.3),
              aes(x0 = x0, y0 = y0, r = r), fill = tint(CYAN, 0.35),
              color = CYAN, linewidth = 0.7) +
  annotate("text", x = 81.6, y = 45.4, label = "UPP1 \u2193", family = FONT,
           fontface = "bold.italic", size = 2.8, color = "grey15") +
  annotate("segment", x = 73.5, y = 45.4, xend = 79.3, yend = 45.4,
           color = GREY, linewidth = 0.9, arrow = arrow_spec) +
  annotate("text", x = 76.4, y = 46.6, label = "virtual KO", family = FONT,
           size = 2.3, color = "grey30") +
  annotate("text", x = 88.6, y = 44.8, hjust = 0.5,
           label = expression(rho == -0.33 * " / " * -0.52),
           family = FONT, size = 2.6, color = "grey25") +
  annotate("text", x = 88.6, y = 43.2, hjust = 0.5,
           label = "monocyte pseudotime", family = FONT, size = 2.3, color = "grey30") +
  annotate("text", x = 75.5, y = 40.2, hjust = 0.5,
           label = "complement activation: only immune process shared by the",
           family = FONT, size = 2.4, color = "grey30") +
  annotate("text", x = 75.5, y = 38.6, hjust = 0.5,
           label = "top-20 dysregulated GSVA processes in both diseases",
           family = FONT, size = 2.4, color = "grey30")

# ============================================================
# RIGHT BOTTOM - drug prediction
# ============================================================
p <- p +
  annotate("rect", xmin = 54.6, xmax = 59.4, ymin = 21.6, ymax = 24.2, fill = BLUE, color = NA) +
  geom_ellipse(data = data.frame(x0 = 57.0, y0 = c(24.2, 21.6),
                                 a = 2.4, b = 0.8, angle = 0),
               aes(x0 = x0, y0 = y0, a = a, b = b, angle = angle),
               fill = BLUE, color = "white", linewidth = 0.4) +
  annotate("text", x = 57.0, y = 20.1, label = "CMap / LINCS", family = FONT, size = 2.7, color = "grey25") +
  annotate("text", x = 57.0, y = 18.8, label = "drug repositioning", family = FONT, size = 2.4, color = "grey30") +
  annotate("segment", x = 60.4, y = 22.9, xend = 63.6, yend = 22.9,
           color = GREY, linewidth = 1.0, arrow = arrow_spec) +
  annotate("text", x = 65.2, y = 25.0, label = "Glucocorticoids", hjust = 0, family = FONT,
           fontface = "bold", size = 3.4, color = BLUE) +
  annotate("text", x = 65.2, y = 23.0, hjust = 0,
           label = "prednisolone \u00b7 dexamethasone \u00b7 hydrocortisone",
           family = FONT, size = 2.6, color = "grey30") +
  annotate("text", x = 65.2, y = 21.3, hjust = 0, label = "current standard of care",
           family = FONT, size = 2.5, color = "grey30", fontface = "italic")

p <- p +
  annotate("text", x = 54.6, y = 16.8, label = "AutoDock Vina docking", hjust = 0,
           family = FONT, fontface = "bold", size = 3.2, color = BLUE)
hx1 <- 57.6 + 1.05 * cos(seq(90, 450, 60) * pi / 180)
hy1 <- 12.2 + 1.05 * sin(seq(90, 450, 60) * pi / 180)
hx2 <- 59.4 + 1.05 * cos(seq(90, 450, 60) * pi / 180)
hy2 <- 12.2 + 1.05 * sin(seq(90, 450, 60) * pi / 180)
p <- p +
  annotate("path", x = hx1, y = hy1, color = "grey15", linewidth = 0.6) +
  annotate("path", x = hx2, y = hy2, color = "grey15", linewidth = 0.6) +
  annotate("segment", x = 60.45, y = 12.72, xend = 61.95, yend = 13.6, color = "grey15", linewidth = 0.6) +
  annotate("segment", x = 61.95, y = 13.6, xend = 63.45, yend = 12.72, color = "grey15", linewidth = 0.6) +
  annotate("text", x = 64.4, y = 12.7, label = "O", family = FONT, size = 2.5, color = RED) +
  annotate("text", x = 56.4, y = 14.2, label = "OH", family = FONT, size = 2.3, color = RED) +
  annotate("text", x = 60.5, y = 9.8, label = "luteolin", family = FONT,
           fontface = "italic", size = 2.7, color = "grey20") +
  annotate("text", x = 86.0, y = 14.2, label = "VNN1 \u2013 luteolin", family = FONT,
           fontface = "bold.italic", size = 3.4, color = TEXT) +
  annotate("text", x = 86.0, y = 12.2, label = "\u221210.01 kcal/mol", family = FONT,
           fontface = "bold", size = 3.1, color = RED) +
  annotate("text", x = 86.0, y = 10.5, label = "strongest predicted binding",
           family = FONT, size = 2.5, color = "grey30") +
  annotate("text", x = 86.0, y = 8.8, label = "candidate drug target", family = FONT,
           size = 2.6, color = BLUE, fontface = "bold")

p <- p +
  annotate("segment", x = 54, y = 7.2, xend = 97, yend = 7.2, color = "grey55", linewidth = 0.5) +
  annotate("text", x = 75.5, y = 5.5, hjust = 0.5,
           label = "Mendelian randomization: no causal signal detected (IVW P = 0.249; forward analysis underpowered)",
           family = FONT, size = 2.7, color = "grey25") +
  annotate("text", x = 75.5, y = 3.6, hjust = 0.5,
           label = "shared disease-response markers and candidate drug targets",
           family = FONT, size = 2.7, color = "grey25", fontface = "italic")

# ============================================================
# export
# ============================================================
dev_quartz_pdf <- function(filename, width, height, ...) {
  grDevices::quartz(type = "pdf", file = filename,
                    width = width, height = height, ...)
}
ggsave(file.path(OUT, "graphical_abstract.pdf"), p,
       width = 15, height = 10.5, device = dev_quartz_pdf)
ggsave(file.path(OUT, "graphical_abstract.png"), p,
       width = 15, height = 10.5, dpi = 300, bg = "white")

tiff_path <- file.path(OUT, "graphical_abstract.tiff")
grDevices::tiff(tiff_path, width = 15, height = 10.5, units = "in",
                res = 300, type = "quartz", compression = "lzw", bg = "white")
print(p)
dev.off()
system(paste0("sips -s dpiWidth 300 -s dpiHeight 300 ", shQuote(tiff_path),
              " >/dev/null 2>&1"))
leftovers <- list.files(OUT, pattern = "^\\.graphical_abstract", all.files = TRUE,
                        full.names = TRUE)
if (length(leftovers)) unlink(leftovers)
cat("GA v2.1 exported to", OUT, "\n")
