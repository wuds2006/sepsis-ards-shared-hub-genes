# ---------------------------------------------------------------------------
# 32_supp_S12.R -- Figure S12. CellChat pathway-level comparison
#   A  pathway strength, sepsis versus healthy, CD14+ monocytes
#   B  summed interaction probability between populations, healthy
#   C  summed interaction probability between populations, sepsis
#   D  pathway strength, sepsis versus ARDS, PBMC cohort
# Input  : data/frozen/08_singlecell/GSE167363/cellchat, 08_singlecell/cellchat
# Output : Figures/supplementary/FigureS12.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

CC1 <- function(...) file.path("data", "frozen", "08_singlecell", "GSE167363", "cellchat", ...)
CC2 <- function(...) file.path("data", "frozen", "08_singlecell", "cellchat", ...)

# --- scatter panels --------------------------------------------------------

scatter_panel <- function(file, xlab, ylab, sub) {
  d <- read.csv(file, stringsAsFactors = FALSE)
  d$lr <- log2((d[[3]] + 1e-6) / (d[[2]] + 1e-6))
  top <- d[order(-abs(d$lr)), ][1:6, ]

  ggplot(d, aes(.data[[names(d)[2]]], .data[[names(d)[3]]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3,
                colour = pal[["slate"]]) +
    geom_point(size = 1.0, colour = pal[["slate"]]) +
    geom_point(data = top, size = 1.4, colour = pal[["brick"]]) +
    ggrepel::geom_text_repel(data = top, aes(label = pathway), size = 7 / .pt,
                             family = FONT_BODY, colour = pal[["ink"]],
                             min.segment.length = 0, segment.size = 0.2,
                             box.padding = 0.25, point.padding = 0.12,
                             direction = "both", max.overlaps = Inf, seed = 2) +
    scale_x_log10(expand = expansion(mult = c(0.04, 0.10))) +
    scale_y_log10(expand = expansion(mult = c(0.04, 0.10))) +
    labs(x = xlab, y = ylab, subtitle = sub) +
    theme_ajr() +
    theme(plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

# --- interaction heatmaps --------------------------------------------------

inter_panel <- function(file, sub) {
  d <- read.csv(file, stringsAsFactors = FALSE)
  ag <- aggregate(prob ~ source + target, data = d, FUN = sum)
  ag$source <- factor(ag$source, levels = sort(unique(ag$source)))
  ag$target <- factor(ag$target, levels = levels(ag$source))

  ggplot(ag, aes(target, source, fill = prob)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_gradientn(colours = ramp_seq, labels = function(v) formatC(v, format = "e", digits = 1),
                         name = "Total\nprobability") +
    scale_x_discrete(labels = function(v) sub("Mono$", " Mo", v)) +
    scale_y_discrete(labels = function(v) sub("Mono$", " Mo", v)) +
    labs(x = "Target", y = "Source", subtitle = sub) +
    theme_ajr() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
          axis.line = element_blank(), axis.ticks = element_blank(),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"),
          legend.key.height = unit(6, "mm"), legend.key.width = unit(2.2, "mm"),
          legend.title = element_text(size = PT_BODY))
}

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid,
                c(list(inter_panel(CC1("HC_interactions.csv"), "Healthy"),
                       inter_panel(CC1("Sepsis_interactions.csv"), "Sepsis"),
                       ncol = 2, labels = c("A", "B")), lab_args))
row2 <- do.call(cowplot::plot_grid,
                c(list(scatter_panel(CC1("pathway_strength_comparison.csv"),
                                     "Pathway strength, healthy",
                                     "Pathway strength, sepsis",
                                     "CD14+ monocytes (GSE167363)"),
                       scatter_panel(CC2("pathway_strength_comparison.csv"),
                                     "Pathway strength, ARDS",
                                     "Pathway strength, sepsis",
                                     "PBMC (GSE151263)"),
                       ncol = 2, labels = c("C", "D")), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 0.95))
save_fig(fig, "FigureS12", dir = out_dir("supplementary"), height_mm = 208)
msg("FigureS12 written")
