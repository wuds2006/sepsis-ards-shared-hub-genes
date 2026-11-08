# ---------------------------------------------------------------------------
# 31_supp_S11.R -- Figure S11. Single-cell atlas of the CD14+ monocyte cohort
#   A  UMAP coloured by monocyte subpopulation
#   B  UMAP split by group
#   C  marker score of each cluster against the annotated populations
# Input  : data/frozen/08_singlecell/GSE167363
# Output : Figures/supplementary/FigureS11.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

SC <- function(...) file.path("data", "frozen", "08_singlecell", "GSE167363", ...)

um <- read.csv(SC("input", "umap.csv"), row.names = 1, stringsAsFactors = FALSE)
mt <- read.csv(SC("input", "meta.csv"), row.names = 1, stringsAsFactors = FALSE)
um <- um[rownames(mt), ]
CT <- c("ClassicalMono", "InflammatoryMono", "IntermediateMono", "ISGMono",
        "MonoDC", "NonclassicalMono", "pDC")
ct_col <- setNames(c(pal[["brick"]], pal[["steel"]], pal[["sage"]], pal[["ochre"]],
                     pal[["plum"]], pal[["olive"]], pal[["slate"]]), CT)

# --- panel A ---------------------------------------------------------------

panel_a <- function() {
  d <- data.frame(um, celltype = factor(mt$celltype, levels = CT))
  cent <- aggregate(cbind(umap_1, umap_2) ~ celltype, data = d, FUN = median)

  ggplot(d, aes(umap_1, umap_2, colour = celltype)) +
    geom_point(size = 0.22, stroke = 0, alpha = 0.9) +
    ggrepel::geom_label_repel(data = cent, aes(umap_1, umap_2, label = celltype),
                              inherit.aes = FALSE, family = FONT_BODY,
                              size = 7 / .pt, fill = "white", alpha = 0.85,
                              label.size = 0, label.padding = unit(0.6, "mm"),
                              min.segment.length = 0, box.padding = 0.2,
                              segment.size = 0.2, max.overlaps = Inf, seed = 3,
                              colour = pal[["ink"]]) +
    scale_colour_manual(values = ct_col, guide = "none") +
    labs(x = "UMAP 1", y = "UMAP 2") +
    coord_fixed(ratio = 1) +
    theme_ajr() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_line(linewidth = 0.3))
}

# --- panel B ---------------------------------------------------------------

panel_b <- function() {
  d <- data.frame(um, grp = factor(ifelse(mt$group == "HC", "Healthy", "Sepsis"),
                                   levels = c("Healthy", "Sepsis")))
  ggplot(d, aes(umap_1, umap_2, colour = grp)) +
    geom_point(size = 0.18, stroke = 0, alpha = 0.9) +
    facet_wrap(~ grp, nrow = 1) +
    scale_colour_manual(values = c(Healthy = pal[["sage"]], Sepsis = pal[["brick"]]),
                        guide = "none") +
    labs(x = "UMAP 1", y = "UMAP 2") +
    coord_fixed(ratio = 1) +
    theme_ajr() +
    theme(strip.text = element_text(size = PT_BODY, hjust = 0),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_line(linewidth = 0.3),
          panel.spacing.x = unit(5, "mm"))
}

# --- panel C ---------------------------------------------------------------

panel_c <- function() {
  m <- read.csv(SC("GSE167363_cluster_marker_scores.csv"), row.names = 1,
                check.names = FALSE)
  d <- as.data.frame(as.table(as.matrix(m)))
  names(d) <- c("cluster", "population", "score")
  d$population <- factor(gsub("1$", "", d$population),
                         levels = sub("1$", "", colnames(m)))
  d$cluster <- factor(d$cluster, levels = sort(unique(as.integer(as.character(d$cluster)))))

  ggplot(d, aes(population, cluster, fill = score)) +
    geom_tile(colour = "white", linewidth = 0.45) +
    scale_fill_gradientn(colours = ramp_seq, name = "Marker\nscore") +
    scale_x_discrete(labels = function(v) stringr::str_wrap(v, 12)) +
    labs(x = NULL, y = "Cluster") +
    theme_ajr() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1,
                                     lineheight = 0.85),
          axis.line = element_blank(), axis.ticks = element_blank(),
          legend.key.height = unit(6, "mm"), legend.key.width = unit(2.2, "mm"),
          legend.title = element_text(size = PT_BODY))
}

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid,
                c(list(panel_a(), panel_b(), ncol = 2, rel_widths = c(0.5, 0.5),
                       labels = c("A", "B")), lab_args))
row2 <- do.call(cowplot::plot_grid, c(list(panel_c(), labels = "C"), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(1.15, 0.85))
save_fig(fig, "FigureS11", dir = out_dir("supplementary"), height_mm = 208)
msg("FigureS11 written")
