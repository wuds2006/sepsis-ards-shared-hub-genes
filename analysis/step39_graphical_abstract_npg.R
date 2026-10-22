# ============================================================
# Graphical Abstract (AJRCCM submission) - ggplot2 + ggsci NPG
# 三区左->右阅读动线: 1 Discovery -> 2 Shared immune
# microenvironment -> 3 Validation & therapeutic relevance
# Journal requirements applied:
#   - >=300 dpi, PDF vector + PNG + TIFF (LZW)
#   - Arial (standard sans-serif), minimal text, no title inside
#   - NPG palette only (ggsci::pal_npg), colorblind-safe accents
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggforce)
  library(ggsci)
})

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "AJRCCM_submission/04_graphical_abstract")

# ---- NPG palette (ggsci) ----
npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
         "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")
COL_P1   <- npg[4]  # #3C5488 dark blue  - discovery
COL_P2   <- npg[3]  # #00A087 teal       - mechanism
COL_P3   <- npg[1]  # #E64B35 red        - translation
COL_GRey <- npg[6]  # #8491B4
COL_CyAn <- npg[2]  # #4DBBD5
COL_CorAl<- npg[5]  # #F39B7F
COL_LT   <- npg[7]  # #91D1C2
TEXT     <- "grey12"

tint <- function(col, p = 0.88) {
  rgb <- col2rgb(col)[, 1]
  rgb(255 - (255 - rgb[1]) * (1 - p),
      255 - (255 - rgb[2]) * (1 - p),
      255 - (255 - rgb[3]) * (1 - p),
      maxColorValue = 255)
}

FONT <- "Arial"

# ---- canvas ----
W <- 100; H <- 60
p <- ggplot() +
  coord_fixed(xlim = c(0, W), ylim = c(0, H), expand = FALSE) +
  theme_void(base_family = FONT) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(0, 0, 0, 0))

# rounded-rect helper (proper polygon, chosen corners rounded)
arc_pts <- function(cx, cy, r, from_deg, to_deg, n = 8) {
  th <- seq(from_deg, to_deg, length.out = n) * pi / 180
  data.frame(x = cx + r * cos(th), y = cy + r * sin(th))
}
rrpoly <- function(id, x, y, w, h, r) {   # all 4 corners rounded
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
rrpoly_top <- function(id, x, y, w, h, r) {  # top corners only
  data.frame(id = id, rbind(
    data.frame(x = x, y = y),
    data.frame(x = x + w, y = y),
    data.frame(x = x + w, y = y + h - r),
    arc_pts(x + w - r, y + h - r, r, 0, 90),
    data.frame(x = x + r, y = y + h),
    arc_pts(x + r, y + h - r, r, 90, 180),
    data.frame(x = x, y = y + h - r)))
}

arrow_spec <- arrow(length = unit(2.6, "mm"), type = "closed",
                    ends = "last")

# ============================================================
# central finding line (single sentence, no article title)
# ============================================================
p <- p +
  annotate("text", x = W/2, y = 57.6, hjust = 0.5,
           label = "Six shared hub genes define the molecular intersection of sepsis and ARDS",
           family = FONT, fontface = "bold", size = 5.6, color = COL_P1)

# ============================================================
# panel frames
# ============================================================
P1 <- rrpoly("p1", 1.5, 2, 30, 51, 0.9)
P2 <- rrpoly("p2", 35, 2, 30, 51, 0.9)
P3 <- rrpoly("p3", 68.5, 2, 30, 51, 0.9)
panel_df <- rbind(P1, P2, P3)
panel_df$fill <- NA
panel_df$fill[panel_df$id == "p1"] <- tint(COL_P1)
panel_df$fill[panel_df$id == "p2"] <- tint(COL_P2)
panel_df$fill[panel_df$id == "p3"] <- tint(COL_P3)

p <- p +
  geom_polygon(data = panel_df,
               aes(x = x, y = y, group = id),
               fill = panel_df$fill, color = NA)

# header bands (rounded top corners, square bottom)
hdr <- rbind(
  rrpoly_top("p1", 1.5, 47.5, 30, 5.5, 0.9),
  rrpoly_top("p2", 35.0, 47.5, 30, 5.5, 0.9),
  rrpoly_top("p3", 68.5, 47.5, 30, 5.5, 0.9))
hdr$col <- NA
hdr$col[hdr$id == "p1"] <- COL_P1
hdr$col[hdr$id == "p2"] <- COL_P2
hdr$col[hdr$id == "p3"] <- COL_P3
p <- p +
  geom_polygon(data = hdr, aes(x = x, y = y, group = id),
               fill = hdr$col, color = NA)

hdr_lab <- data.frame(
  x = c(16.5, 50, 83.5),
  y = rep(50.25, 3),
  lab = c("1  |  Identification of shared hub genes",
          "2  |  Shared immune microenvironment",
          "3  |  Validation and therapeutic relevance"))
p <- p +
  geom_text(data = hdr_lab, aes(x = x, y = y, label = lab),
            family = FONT, fontface = "bold", size = 4.3,
            color = "white")

# ---- inter-panel flow arrows ----
p <- p +
  annotate("segment", x = 31.8, y = 28, xend = 34.7, yend = 28,
           color = COL_GRey, linewidth = 1.1, arrow = arrow_spec) +
  annotate("segment", x = 65.3, y = 28, xend = 68.2, yend = 28,
           color = COL_GRey, linewidth = 1.1, arrow = arrow_spec) +
  annotate("text", x = 33.25, y = 30.2, label = "shared\nhub genes",
           family = FONT, size = 2.7, color = COL_GRey, fontface = "bold") +
  annotate("text", x = 66.75, y = 30.2, label = "immune\ncontext",
           family = FONT, size = 2.7, color = COL_GRey, fontface = "bold")

# ============================================================
# PANEL 1 - discovery
# ============================================================
cx1 <- 16.5  # panel-1 center x

# cohort boxes
cohort <- rbind(
  rrpoly("c1", 3.5, 39, 12.5, 8, 0.5),
  rrpoly("c2", 17.5, 39, 12.5, 8, 0.5))
p <- p +
  geom_polygon(data = cohort, aes(x = x, y = y, group = id),
               fill = "white", color = COL_GRey, linewidth = 0.5) +
  annotate("text", x = 9.75, y = 45.3, label = "GSE76293",
           family = FONT, fontface = "bold.italic", size = 3.6, color = TEXT) +
  annotate("text", x = 23.75, y = 45.3, label = "GSE69528",
           family = FONT, fontface = "bold.italic", size = 3.6, color = TEXT) +
  annotate("text", x = 9.75, y = 43.3, label = "ARDS \u00b7 PMNs",
           family = FONT, size = 3.1, color = TEXT) +
  annotate("text", x = 23.75, y = 43.3, label = "Sepsis \u00b7 whole blood",
           family = FONT, size = 3.1, color = TEXT) +
  annotate("text", x = 9.75, y = 41.4, label = "n = 12 vs 12",
           family = FONT, size = 3.1, color = "grey35") +
  annotate("text", x = 23.75, y = 41.4, label = "n = 83 vs 28",
           family = FONT, size = 3.1, color = "grey35") +
  # DEG / module counts
  annotate("text", x = 9.75, y = 37.2, label = "557 DEGs · 13 modules",
           family = FONT, size = 2.9, color = "grey30") +
  annotate("text", x = 23.75, y = 37.2, label = "1,458 DEGs · 11 modules",
           family = FONT, size = 2.9, color = "grey30")

# flow: cohorts -> intersection -> 179 genes -> ML -> hub genes
p <- p +
  annotate("segment", x = 9.75, y = 36.4, xend = 13.9, yend = 34.4,
           color = COL_P1, linewidth = 0.7, arrow = arrow_spec) +
  annotate("segment", x = 23.75, y = 36.4, xend = 19.1, yend = 34.4,
           color = COL_P1, linewidth = 0.7, arrow = arrow_spec) +
  annotate("text", x = cx1, y = 33.1,
           label = "DEGs \u2229 disease-associated WGCNA modules",
           family = FONT, size = 3.0, color = TEXT) +
  annotate("segment", x = cx1, y = 32.2, xend = cx1, yend = 30.9,
           color = COL_P1, linewidth = 0.7, arrow = arrow_spec) +
  annotate("text", x = cx1, y = 29.6, label = "179 shared candidate genes",
           family = FONT, fontface = "bold", size = 3.5, color = COL_P1) +
  annotate("segment", x = cx1, y = 28.7, xend = cx1, yend = 27.4,
           color = COL_P1, linewidth = 0.7, arrow = arrow_spec) +
  annotate("text", x = cx1, y = 26.1, label = "11 machine-learning classifiers + SHAP",
           family = FONT, size = 3.1, color = TEXT) +
  annotate("segment", x = cx1, y = 25.2, xend = cx1, yend = 23.9,
           color = COL_P1, linewidth = 0.7, arrow = arrow_spec)

# hub-gene box (white, thick NPG border)
hub <- rrpoly("hub", 4.5, 15.6, 24, 7.3, 0.5)
p <- p +
  geom_polygon(data = hub, aes(x = x, y = y),
               fill = "white", color = COL_P1, linewidth = 0.9) +
  annotate("text", x = cx1, y = 21.7, label = "Six hub genes",
           family = FONT, fontface = "bold", size = 3.8, color = COL_P1) +
  annotate("text", x = cx1, y = 19.6,
           label = "CDK5RAP2    CD177    GRB10",
           family = FONT, fontface = "bold.italic", size = 3.4, color = TEXT) +
  annotate("text", x = cx1, y = 17.3,
           label = "UPP1    EXOSC4    VNN1",
           family = FONT, fontface = "bold.italic", size = 3.4, color = TEXT)

# mini bar chart: cross-validated AUC of ML models
bar_x <- c(8.2, 10.6, 13.0, 15.4, 17.8)
bar_h <- c(0.75, 0.85, 0.92, 0.97, 1.00)
base_y <- 5.2; bar_w <- 1.6; hscale <- 7
bar_cols <- c(COL_GRey, COL_GRey, COL_CyAn, COL_CyAn, COL_P1)
p <- p +
  annotate("segment", x = 6.6, y = base_y, xend = 19.4, yend = base_y,
           color = "grey45", linewidth = 0.5) +
  annotate("rect",
           xmin = bar_x - bar_w/2, xmax = bar_x + bar_w/2,
           ymin = base_y, ymax = base_y + (bar_h - 0.6) * hscale,
           fill = bar_cols, color = NA) +
  annotate("text", x = cx1, y = 14.0,
           label = "Cross-validated AUC 0.75-1.00 (both cohorts)",
           family = FONT, size = 2.8, color = "grey30") +
  annotate("text", x = cx1, y = 3.6,
           label = "confirmed in GSE243066 (ARDS) and GSE100159 (sepsis)",
           family = FONT, size = 2.8, color = "grey30", fontface = "italic")

# ============================================================
# PANEL 2 - shared immune microenvironment
# ============================================================
cx2 <- 50

p <- p +
  # M2 macrophages
  annotate("text", x = cx2, y = 45.3, label = "M2 macrophage expansion in sepsis blood",
           family = FONT, fontface = "bold", size = 3.3, color = COL_P2) +
  annotate("text", x = cx2, y = 43.2, parse = TRUE,
           label = "P~'='~5.8 %*% 10^-14~'(ssGSEA, GSE69528)'",
           family = FONT, size = 2.9, color = "grey30") +
  # complement
  annotate("segment", x = 39, y = 41.4, xend = 61, yend = 41.4,
           color = "white", linewidth = 1.4) +
  annotate("text", x = cx2, y = 39.6, label = "Complement activation",
           family = FONT, fontface = "bold", size = 3.3, color = COL_P2) +
  annotate("text", x = cx2, y = 37.6,
           label = "the only immune pathway shared by the top-20 dysregulated\nprocesses in both diseases (GSVA)",
           family = FONT, size = 2.9, color = "grey30", lineheight = 0.95) +
  annotate("segment", x = 39, y = 35.6, xend = 61, yend = 35.6,
           color = "white", linewidth = 1.4) +
  # single-cell block
  annotate("text", x = cx2, y = 33.8, label = "Single-cell transcriptomics (PBMC + CD14\u207a monocytes)",
           family = FONT, fontface = "bold", size = 3.2, color = COL_P2) +
  annotate("text", x = cx2, y = 31.6,
           label = expression(UPP1~declines~along~monocyte~pseudotime~(rho == -0.33 * "/" * -0.52)),
           family = FONT, size = 2.9, color = "grey25") +
  annotate("text", x = cx2, y = 29.6,
           label = expression(CEBPB~dominant~upstream~regulator * ";"~virtual~knockout%down%UPP1),
           family = FONT, size = 2.9, color = "grey25") +
  annotate("text", x = cx2, y = 27.6,
           label = "scTenifoldKnk + CellOracle · CellChat · slingshot",
           family = FONT, size = 2.7, color = "grey40", fontface = "italic")

# pseudotime mini-trajectory: UPP1 high -> low along differentiation
set.seed(11)
n <- 90
tt <- runif(n)
tx <- 40 + tt * 18
ty <- 22.5 - tt * 1.6 + sin(tt * 4.5) * 0.9 + rnorm(n, 0, 0.55)
cell_col <- colorRampPalette(c(COL_P2, COL_LT, "#D9EDE7"))(100)[pmin(100, pmax(1, round(1 + tt * 99)))]
p <- p +
  annotate("path", x = 40 + seq(0, 18, 0.2),
           y = 22.5 - seq(0, 18, 0.2)/11.25 + sin(seq(0, 18, 0.2) * 0.25) * 0.9,
           color = COL_P2, linewidth = 0.7, alpha = 0.55) +
  annotate("point", x = tx, y = ty, size = 1.15, color = cell_col, alpha = 0.85) +
  annotate("text", x = 40.4, y = 24.6, label = "UPP1 high",
           family = FONT, size = 2.7, color = COL_P2, fontface = "bold.italic") +
  annotate("text", x = 57.8, y = 20.2, label = "UPP1 low",
           family = FONT, size = 2.7, color = "grey35", fontface = "bold.italic") +
  annotate("text", x = cx2, y = 17.6, label = "monocyte differentiation \u2192",
           family = FONT, size = 2.7, color = "grey40", fontface = "italic")

# normal - sepsis - ARDS axis with inverted-U metabolic activity
ax_y <- 12.2
p <- p +
  annotate("segment", x = 39, y = ax_y, xend = 61.4, yend = ax_y,
           color = COL_GRey, linewidth = 1.0,
           arrow = arrow(length = unit(2.6, "mm"), type = "closed")) +
  annotate("text", x = 40.5, y = ax_y - 1.6, label = "health",
           family = FONT, size = 2.9, color = "grey35") +
  annotate("text", x = 50, y = ax_y - 1.6, label = "sepsis",
           family = FONT, size = 2.9, color = "grey35") +
  annotate("text", x = 59.9, y = ax_y - 1.6, label = "ARDS",
           family = FONT, size = 2.9, color = "grey35")
# inverted-U curve above axis
ux <- seq(39, 61, 0.25)
uy <- ax_y + 3.6 * exp(-((ux - 50)^2) / (2 * 6.2^2))
p <- p +
  annotate("path", x = ux, y = uy, color = COL_P2, linewidth = 1.1) +
  annotate("text", x = 50, y = ax_y + 4.6,
           label = "UPP1 program and myeloid metabolic activity peak at sepsis, partially recede in ARDS",
           family = FONT, size = 2.8, color = "grey30")

# ============================================================
# PANEL 3 - validation & therapeutic relevance
# ============================================================
cx3 <- 83.5

# external validation: transferred six-gene model AUCs (horizontal bars)
val_y <- c(44.6, 41.8); val_h <- 1.9
val_len <- c(0.984, 0.708) * 20
p <- p +
  annotate("text", x = cx3, y = 46.6, label = "External validation (six-gene model)",
           family = FONT, fontface = "bold", size = 3.2, color = COL_P3) +
  annotate("rect", xmin = 71, xmax = 71 + val_len[1],
           ymin = val_y[1] - val_h/2, ymax = val_y[1] + val_h/2,
           fill = COL_P3, color = NA) +
  annotate("rect", xmin = 71, xmax = 71 + val_len[2],
           ymin = val_y[2] - val_h/2, ymax = val_y[2] + val_h/2,
           fill = COL_CorAl, color = NA) +
  annotate("text", x = 70.4, y = val_y[1], label = "ARDS",
           family = FONT, size = 2.9, hjust = 1, color = "grey30") +
  annotate("text", x = 70.4, y = val_y[2], label = "sepsis",
           family = FONT, size = 2.9, hjust = 1, color = "grey30") +
  annotate("text", x = 71 + val_len[1] + 0.7, y = val_y[1],
           label = "AUC 0.984", family = FONT, size = 2.9, hjust = 0,
           fontface = "bold", color = COL_P3) +
  annotate("text", x = 71 + val_len[2] + 0.7, y = val_y[2],
           label = "AUC 0.708", family = FONT, size = 2.9, hjust = 0,
           color = "grey30") +
  annotate("text", x = cx3, y = 40.0,
           label = "single-gene AUC 0.753-0.998 across GSE243066 and GSE100159",
           family = FONT, size = 2.8, color = "grey30") +
  annotate("segment", x = 71.5, y = 38.4, xend = 95.5, yend = 38.4,
           color = "white", linewidth = 1.4)

# MR null
p <- p +
  annotate("text", x = cx3, y = 36.6, label = "Mendelian randomization: no causal effect",
           family = FONT, fontface = "bold", size = 3.2, color = COL_P3) +
  annotate("text", x = cx3, y = 34.7,
           label = "IVW P = 0.249; LDSC rg \u2248 0",
           family = FONT, size = 2.9, color = "grey30") +
  annotate("text", x = cx3, y = 32.9,
           label = "shared disease-response markers, not causal drivers",
           family = FONT, size = 2.9, color = "grey30", fontface = "italic") +
  annotate("segment", x = 71.5, y = 31.6, xend = 95.5, yend = 31.6,
           color = "white", linewidth = 1.4)

# drug repositioning
p <- p +
  annotate("text", x = cx3, y = 30.6, label = "CMap drug repositioning",
           family = FONT, fontface = "bold", size = 3.2, color = COL_P3) +
  annotate("text", x = cx3, y = 28.5, label = "glucocorticoids - current standard of care",
           family = FONT, size = 2.9, color = "grey30") +
  annotate("segment", x = 71.5, y = 26.4, xend = 95.5, yend = 26.4,
           color = "white", linewidth = 1.4)

# docking: VNN1 - luteolin
p <- p +
  annotate("text", x = cx3, y = 24.6, label = "AutoDock Vina docking",
           family = FONT, fontface = "bold", size = 3.2, color = COL_P3) +
  annotate("text", x = cx3, y = 22.5, label = "VNN1 - luteolin",
           family = FONT, fontface = "bold.italic", size = 3.3, color = TEXT) +
  annotate("text", x = cx3, y = 20.5, label = "\u221210.01 kcal/mol (strongest predicted binding)",
           family = FONT, size = 2.9, color = "grey30")

# conclusion box
concl <- rrpoly("concl", 71, 6.5, 25, 9.7, 0.6)
p <- p +
  geom_polygon(data = concl, aes(x = x, y = y),
               fill = COL_P1, color = NA) +
  annotate("text", x = cx3, y = 13.6, label = "Markers of disease state",
           family = FONT, fontface = "bold", size = 3.4, color = "white") +
  annotate("text", x = cx3, y = 11.2, label = "Candidate drug targets",
           family = FONT, fontface = "bold", size = 3.4, color = "white") +
  annotate("text", x = cx3, y = 8.4,
           label = "for a comorbidity that remains common and lethal",
           family = FONT, size = 2.8, color = "grey85", fontface = "italic")

# footer method chain
p <- p +
  annotate("text", x = W/2, y = 0.8,
           label = "Multi-omics integration  \u00b7  machine learning  \u00b7  single-cell transcriptomics  \u00b7  Mendelian randomization  \u00b7  drug repositioning and docking",
           family = FONT, size = 2.8, color = "grey45", fontface = "italic")

# ============================================================
# export: PDF vector + PNG 300 dpi + TIFF 300 dpi (LZW)
# ============================================================
dev_quartz_pdf <- function(filename, width, height, ...) {
  grDevices::quartz(type = "pdf", file = filename,
                    width = width, height = height, ...)
}
ggsave(file.path(OUT, "graphical_abstract.pdf"), p,
       width = 16, height = 9.6, device = dev_quartz_pdf)
ggsave(file.path(OUT, "graphical_abstract.png"), p,
       width = 16, height = 9.6, dpi = 300, bg = "white")
ggsave(file.path(OUT, "graphical_abstract.tiff"), p,
       width = 16, height = 9.6, dpi = 300, bg = "white",
       device = grDevices::tiff, compression = "lzw", type = "quartz")

cat("GA exported to", OUT, "\n")

# fix TIFF dpi metadata (quartz device writes 72; enforce 300 dpi)
tiff_path <- file.path(OUT, "graphical_abstract.tiff")
if (file.exists(tiff_path)) {
  system(paste0("sips -s dpiWidth 300 -s dpiHeight 300 ", shQuote(tiff_path),
                " >/dev/null 2>&1"))
}
