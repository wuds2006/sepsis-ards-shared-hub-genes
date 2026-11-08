# ---------------------------------------------------------------------------
# 37_supp_S17.R -- Figure S17. scTenifoldKnk virtual knockout of the hub genes
#   A  number of significantly re-regulated genes per knockout
#   B  distribution of the differential-regulation statistic per knockout
#   C  genes with the strongest regulatory shift, across the four knockouts
# Input  : data/frozen/08_singlecell/GSE167363/knockout
# Output : Figures/supplementary/FigureS17.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

KO <- function(...) file.path("data", "frozen", "08_singlecell", "GSE167363", "knockout", ...)
GENES <- c("CDK5RAP2", "EXOSC4", "UPP1", "VNN1")
ko_col <- setNames(c(pal[["brick"]], pal[["steel"]], pal[["sage"]], pal[["ochre"]]), GENES)

d <- do.call(rbind, lapply(GENES, function(g) {
  x <- read.csv(KO(sprintf("KO_%s_diffRegulation.csv", g)), stringsAsFactors = FALSE)
  x$ko <- g
  x
}))
d$ko <- factor(d$ko, levels = GENES)
d$sig <- !is.na(d$p.adj) & d$p.adj < 0.05

# --- panel A ---------------------------------------------------------------

panel_a <- function() {
  a <- data.frame(ko = factor(GENES, levels = GENES),
                  n = as.numeric(tapply(d$sig, d$ko, sum)))
  ggplot(a, aes(n, ko, fill = ko)) +
    geom_col(width = 0.62, colour = NA) +
    geom_text(aes(label = n), hjust = -0.3, family = FONT_BODY,
              size = 7 / .pt, colour = pal[["ink"]]) +
    scale_fill_manual(values = ko_col, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = "Genes with adjusted P < 0.05", y = NULL) +
    theme_ajr() +
    theme(axis.text.y = element_text(face = "italic"))
}

# --- panel B ---------------------------------------------------------------

panel_b <- function() {
  ggplot(d, aes(ko, Z, fill = ko)) +
    geom_violin(width = 0.72, linewidth = 0.25, colour = "#4A4A4A", scale = "width") +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", linewidth = 0.25) +
    scale_fill_manual(values = ko_col, guide = "none") +
    scale_x_discrete(labels = GENES) +
    labs(x = NULL, y = "Differential-regulation statistic (Z)") +
    theme_ajr() +
    theme(axis.text.x = element_text(face = "italic"))
}

# --- panel C ---------------------------------------------------------------

panel_c <- function() {
  top <- d[order(-abs(d$Z)), ]
  keep <- unique(top$gene)[1:15]
  m <- d[d$gene %in% keep, c("gene", "ko", "Z")]
  m$Z <- pmax(pmin(m$Z, 8), -8)
  m$gene <- factor(m$gene, levels = rev(keep))

  ggplot(m, aes(ko, gene, fill = Z)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_gradientn(colours = ramp_div, limits = c(-8, 8),
                         oob = scales::squish, name = "Z") +
    scale_x_discrete(labels = GENES) +
    labs(x = NULL, y = NULL) +
    theme_ajr() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1,
                                     face = "italic"),
          axis.line = element_blank(), axis.ticks = element_blank(),
          legend.key.height = unit(6, "mm"), legend.key.width = unit(2.2, "mm"),
          legend.title = element_text(size = PT_BODY))
}

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid,
                c(list(panel_a(), panel_b(), ncol = 2, labels = c("A", "B")), lab_args))
row2 <- do.call(cowplot::plot_grid, c(list(panel_c(), labels = "C"), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(0.9, 1.1))
save_fig(fig, "FigureS17", dir = out_dir("supplementary"), height_mm = 204)
msg("FigureS17 written")
