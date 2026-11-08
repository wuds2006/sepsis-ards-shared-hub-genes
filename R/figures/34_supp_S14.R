# ---------------------------------------------------------------------------
# 34_supp_S14.R -- Figure S14. Monocle3 trajectory analysis in both cohorts
#   A  pseudotime distribution in the PBMC cohort
#   B  pseudotime distribution in the CD14+ monocyte cohort
#   C  Moran's I of each hub gene along pseudotime
# Input  : data/frozen/08_singlecell/monocle3
# Output : Figures/supplementary/FigureS14.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

MC <- function(...) file.path("data", "frozen", "08_singlecell", "monocle3", ...)

dist_panel <- function(cd, sub) {
  d <- read.csv(MC(sprintf("monocle3_pseudotime_%s.csv", cd)), stringsAsFactors = FALSE)
  w <- read.csv(MC(sprintf("monocle3_pseudotime_wilcox_%s.csv", cd)), stringsAsFactors = FALSE)
  d$grp <- factor(ifelse(d$group == "HC", "Healthy", d$group),
                  levels = c("Healthy", "ARDS", "Sepsis"))
  d$grp <- droplevels(d$grp)
  plab <- sprintf("P = %s", formatC(w$p.value[1], format = "e", digits = 1))

  ggplot(d, aes(pseudotime, grp, fill = grp)) +
    geom_violin(width = 0.72, linewidth = 0.25, colour = "#4A4A4A", scale = "width") +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.25) +
    annotate("text", x = 0, y = Inf, vjust = 1.8, hjust = 0, label = plab,
             family = FONT_BODY, size = 7 / .pt, colour = pal[["ink"]]) +
    scale_fill_manual(values = c(Healthy = pal[["sage"]], ARDS = pal[["steel"]],
                                 Sepsis = pal[["brick"]]), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.06))) +
    labs(x = "monocle3 pseudotime", y = NULL, subtitle = sub) +
    theme_ajr() +
    theme(plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

# --- panel C ---------------------------------------------------------------

panel_c <- function() {
  f <- function(cd, lab) {
    x <- read.csv(MC(sprintf("monocle3_hub_moransI_%s.csv", cd)), stringsAsFactors = FALSE)
    x$cohort <- lab
    x
  }
  d <- rbind(f("GSE151263", "PBMC (GSE151263)"),
             f("GSE167363", "CD14+ monocytes (GSE167363)"))
  d <- d[d$gene %in% HUB6, ]
  d$gene <- factor(d$gene, levels = rev(HUB6))

  ggplot(d, aes(morans_I, gene, fill = cohort)) +
    geom_col(position = position_dodge(width = 0.68), width = 0.62, colour = NA) +
    scale_fill_manual(values = c(`PBMC (GSE151263)` = pal[["steel"]],
                                 `CD14+ monocytes (GSE167363)` = pal[["brick"]]),
                      name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Moran's I along pseudotime", y = NULL) +
    theme_ajr() +
    theme(axis.text.y = element_text(face = "italic"),
          legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"))
}

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid,
                c(list(dist_panel("GSE151263", "PBMC (GSE151263)"),
                       dist_panel("GSE167363", "CD14+ monocytes (GSE167363)"),
                       ncol = 2, labels = c("A", "B")), lab_args))
row2 <- do.call(cowplot::plot_grid, c(list(panel_c(), labels = "C"), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(0.95, 1.05))
save_fig(fig, "FigureS14", dir = out_dir("supplementary"), height_mm = 202)
msg("FigureS14 written")
