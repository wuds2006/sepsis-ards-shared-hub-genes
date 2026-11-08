# ---------------------------------------------------------------------------
# 26_supp_S7.R -- Figure S7. Per-cohort principal-component analysis
#   batch-effect assessment within each bulk cohort
# Input  : data/frozen/06_external_validation/reviewer_w5/cohort_PCA_coords.csv
# Output : Figures/supplementary/FigureS7.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

p <- read.csv(file.path("data", "frozen", "06_external_validation", "reviewer_w5",
                        "cohort_PCA_coords.csv"), stringsAsFactors = FALSE)

coh <- c(GSE76293 = "ARDS discovery (GSE76293)",
         GSE69528 = "Sepsis discovery (GSE69528)",
         GSE243066 = "ARDS validation (GSE243066)",
         GSE100159 = "Sepsis validation (GSE100159)")
p$coh <- factor(coh[p$cohort], levels = coh)
p$grp <- factor(ifelse(p$group == "Control", "Control", "Disease"),
                levels = c("Control", "Disease"))

pca_panel <- function(k) {
  d <- p[p$coh == k, ]
  n <- table(d$grp)
  ggplot(d, aes(PC1, PC2, colour = grp)) +
    geom_point(size = 1.0) +
    scale_colour_manual(values = c(Control = pal[["sage"]], Disease = pal[["brick"]]),
                        labels = c(sprintf("Control (n = %d)", n[["Control"]]),
                                   sprintf("Disease (n = %d)", n[["Disease"]])),
                        name = NULL) +
    labs(x = "PC1", y = "PC2", subtitle = k) +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid,
                c(list(pca_panel(coh[1]), pca_panel(coh[2]), ncol = 2,
                       labels = c("A", "B")), lab_args))
row2 <- do.call(cowplot::plot_grid,
                c(list(pca_panel(coh[3]), pca_panel(coh[4]), ncol = 2,
                       labels = c("C", "D")), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1)
save_fig(fig, "FigureS7", dir = out_dir("supplementary"), height_mm = 200)
msg("FigureS7 written")
