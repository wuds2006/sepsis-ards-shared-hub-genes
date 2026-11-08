# ---------------------------------------------------------------------------
# 45_supp_S25.R -- Supplementary Figure S25
#   A  2D protein-ligand interaction diagram, VNN1-luteolin (7SLV, chain A)
#   B  3D binding pose in the VNN1 pocket (cartoon + pocket sticks)
#
# The 2D skeleton comes from RDKit (CoordGen layout, Kekule bond orders) so that
# the structure is luteolin (C15H10O6) rather than a projection of 3D coordinates;
# every contact line is anchored on the ligand atom that PLIP assigned to it.
# Input  : R/figures/_figS25_work/{atoms2d,bonds2d,contacts2d}.csv, render3d.png
# Output : Figures/supplementary/FigureS25.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")
suppressPackageStartupMessages({library(png); library(grid)})

WK <- file.path("R", "figures", "_figS25_work")

A <- read.csv(file.path(WK, "atoms2d.csv"), stringsAsFactors = FALSE)
B <- read.csv(file.path(WK, "bonds2d.csv"), stringsAsFactors = FALSE)
C <- read.csv(file.path(WK, "contacts2d.csv"), stringsAsFactors = FALSE)
RG <- read.csv(file.path(WK, "rings2d.csv"), stringsAsFactors = FALSE)

# --- interaction colours: identical to Figure 10E ---------------------------
TYPE_COL <- c(`hydrogen bond` = pal[["brick"]],
              hydrophobic   = pal[["olive"]],
              `pi-stacking` = pal[["steel"]])
TYPE_SHAPE <- c(`hydrogen bond` = 16, hydrophobic = 15, `pi-stacking` = 17)

# --- hand-placed label anchors (molecule coordinates) -----------------------
# xy = where the residue name sits; ax/ay = the ligand atom it points to.
LBL <- data.frame(
  residue  = c("Glu164", "Lys157", "Met218", "Trp64", "Phe161", "Phe191", "Trp217", "Phe161"),
  type     = c("hydrogen bond", "hydrogen bond", "hydrogen bond",
               "hydrophobic", "hydrophobic", "hydrophobic",
               "pi-stacking", "pi-stacking"),
  dist     = c(2.87, 2.90, 3.42, 3.12, 3.72, 3.38, 4.57, 3.87),
  ax       = c(-2.22, -0.49, 3.84, -3.09, 0.37, 2.10, 2.10, 2.10),
  ay       = c(2.43, 2.43, -0.07, 0.93, 0.93, -2.07, -1.07, -1.07),
  x        = c(-3.15, 0.05, 4.75, -4.85, 1.05, 2.55, 4.10, 4.10),
  y        = c(3.35, 3.35, 1.00, 1.75, 2.05, -3.10, -2.35, -1.20),
  hjust    = c(0.5, 0.5, 0, 0.5, 0, 0.5, 0, 0),
  stringsAsFactors = FALSE
)
LBL$lab <- sprintf("%s\n%.2f \u00c5", LBL$residue, LBL$dist)

# --- panel A ----------------------------------------------------------------
ring_centroid <- function(rid) {
  ids <- RG$idx0[RG$ring == rid] + 1
  c(mean(A$x[ids]), mean(A$y[ids]))
}

panel_a <- function() {
  # 骨架：所有键先画一条完整实线，环因此闭合
  skel <- data.frame(x = A$x[B$i + 1], y = A$y[B$i + 1],
                     xend = A$x[B$j + 1], yend = A$y[B$j + 1])

  # 双键：第二条线画在环内侧（非环双键取任一侧），两端各缩进 18%
  dbl <- B[B$order == 2, ]
  inner <- do.call(rbind, lapply(seq_len(nrow(dbl)), function(k) {
    i0 <- dbl$i[k] + 1; i1 <- dbl$j[k] + 1
    x0 <- A$x[i0]; y0 <- A$y[i0]; x1 <- A$x[i1]; y1 <- A$y[i1]
    dx <- x1 - x0; dy <- y1 - y0; L <- sqrt(dx^2 + dy^2)
    px <- -dy / L; py <- dx / L                       # 单位法向

    # 该键所属的环（两原子同环）
    rid <- RG$ring[RG$idx0 == dbl$i[k]]
    rid <- rid[rid %in% RG$ring[RG$idx0 == dbl$j[k]]]
    if (length(rid)) {
      cen <- ring_centroid(rid[1])
      mx <- (x0 + x1) / 2; my <- (y0 + y1) / 2
      if ((cen[1] - mx) * px + (cen[2] - my) * py < 0) { px <- -px; py <- -py }
    }
    off <- 0.085
    tx <- dx / L * 0.18; ty <- dy / L * 0.18          # 端部缩进
    data.frame(x = x0 + px * off + tx, y = y0 + py * off + ty,
               xend = x1 + px * off - tx, yend = y1 + py * off - ty)
  }))

  het <- A[A$symbol != "C", ]
  het$lab <- ifelse(het$nH > 0, "OH", "O")

  ggplot() +
    # contacts first so the skeleton sits on top
    geom_segment(data = LBL, aes(x = ax, y = ay, xend = x, yend = y,
                                 colour = type),
                 linetype = "22", linewidth = 0.32, show.legend = FALSE) +
    geom_segment(data = skel, aes(x = x, y = y, xend = xend, yend = yend),
                 linewidth = 0.45, colour = pal[["ink"]], lineend = "round") +
    geom_segment(data = inner, aes(x = x, y = y, xend = xend, yend = yend),
                 linewidth = 0.45, colour = pal[["ink"]], lineend = "round") +
    geom_point(data = LBL, aes(x = ax, y = ay, colour = type, shape = type),
               size = 0.9, show.legend = FALSE) +
    geom_label(data = het, aes(x = x, y = y, label = lab), size = PT_BODY / .pt,
               fontface = "plain", label.size = 0, label.padding = unit(0.6, "pt"),
               fill = "white", colour = pal[["ink"]], family = FONT_BODY) +
    geom_label(data = LBL, aes(x = x, y = y, label = lab, colour = type),
               size = PT_BODY / .pt, lineheight = 0.95, hjust = LBL$hjust,
               label.size = 0, label.padding = unit(0.9, "pt"),
               fill = "white", family = FONT_BODY, show.legend = FALSE) +
    scale_colour_manual(values = TYPE_COL) +
    scale_shape_manual(values = TYPE_SHAPE) +
    coord_equal(clip = "off") +
    theme_void_ajr() +
    theme(plot.margin = margin(4, 6, 2, 6))
}

# --- panel B ----------------------------------------------------------------
panel_b <- function() {
  img <- png::readPNG(file.path(WK, "render3d.png"))
  g <- rasterGrob(img, interpolate = TRUE)
  ggplot() +
    annotation_custom(g, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void_ajr() +
    theme(plot.margin = margin(2, 2, 2, 2))
}

# --- assemble ---------------------------------------------------------------
# 两个面板都用 coord_equal，因此面板框的纵横比必须与内容一致，否则会出现
# 上下或左右留白（letterbox）。纵横比由内容决定：
#   A 约 9.9 x 6.5 分子坐标单位；B 由 3D 渲染画布 2400 x 1800 决定
ASP_A <- 1.52
ASP_B <- 2400 / 1800
H_A   <- FIG_W_MM / ASP_A
H_B   <- FIG_W_MM / ASP_B

fig <- cowplot::plot_grid(panel_a(), panel_b(), ncol = 1,
                          rel_heights = c(H_A, H_B),
                          labels = c("A", "B"),
                          label_size = PT_TAG, label_fontfamily = FONT_TAG,
                          label_fontface = "bold", label_colour = pal[["ink"]])

save_fig(fig, "FigureS25", width_mm = FIG_W_MM, height_mm = H_A + H_B + 3,
         dir = out_dir("supplementary"))
msg(sprintf("FigureS25 written (%.0f mm x %.0f mm)", FIG_W_MM, H_A + H_B + 3))
