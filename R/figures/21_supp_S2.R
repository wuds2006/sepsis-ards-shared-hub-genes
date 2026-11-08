# ---------------------------------------------------------------------------
# 21_supp_S2.R -- Figure S2. Donor-level sensitivity of the single-cell contrasts
#   A  monocle3 pseudotime per donor
#   B  UPP1 programme score per donor
#   C  metabolic programme activity per donor
#   D  cells recovered per donor
# Input  : data/frozen/08_singlecell/donor_sensitivity, 08_singlecell/metabolism
# Output : Figures/supplementary/FigureS2.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

DN <- function(...) file.path("data", "frozen", "08_singlecell", "donor_sensitivity", ...)

d <- read.csv(DN("donor_level_summary.csv"), stringsAsFactors = FALSE)
d$grp <- factor(d$group, levels = c("HC", "ARDS", "Sepsis"))
d$dataset <- factor(d$dataset,
                    levels = c("GSE151263", "GSE167363"),
                    labels = c("PBMC\nGSE151263", "CD14+ monocytes\nGSE167363"))
gcol <- c(HC = pal[["sage"]], ARDS = pal[["steel"]], Sepsis = pal[["brick"]])

donor_panel <- function(metric, ylab, sub) {
  x <- d[d$metric == metric, ]
  ggplot(x, aes(stats::reorder(donor, value), value, colour = grp)) +
    geom_point(size = 1.3) +
    facet_wrap(~ dataset, scales = "free_x", nrow = 1) +
    scale_colour_manual(values = gcol, name = NULL) +
    labs(x = NULL, y = ylab, subtitle = sub) +
    theme_ajr() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
          legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

panel_a <- function() donor_panel("pseudotime", "Donor median pseudotime",
                                  "Monocle3 pseudotime")
panel_b <- function() donor_panel("UPP1_program", "Donor median AUCell score",
                                  "UPP1 co-expression programme")

# --- panel C ---------------------------------------------------------------

panel_c <- function() {
  x <- d[d$metric == "metabolism", ]
  x$pathway <- factor(x$pathway,
                      levels = c("Glycolysis", "Oxidative phosphorylation",
                                 "Fatty acid metabolism", "Xenobiotic metabolism"))
  ggplot(x, aes(value, pathway, colour = grp)) +
    geom_point(size = 1.2, position = position_jitter(height = 0.16, seed = 4)) +
    facet_wrap(~ dataset, nrow = 1) +
    scale_colour_manual(values = gcol, name = NULL) +
    scale_y_discrete(labels = function(v) stringr::str_wrap(v, 18)) +
    labs(x = "Donor median programme score", y = NULL, subtitle = "Metabolic programmes") +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"),
          axis.text.y = element_text(lineheight = 0.86),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

# --- panel D ---------------------------------------------------------------

panel_d <- function() {
  x <- unique(d[!is.na(d$n_cells), c("donor", "group", "n_cells", "dataset")])
  x$grp <- factor(x$group, levels = c("HC", "ARDS", "Sepsis"))
  ggplot(x, aes(stats::reorder(donor, n_cells), n_cells, fill = grp)) +
    geom_col(width = 0.66, colour = NA) +
    facet_wrap(~ dataset, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = gcol, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = NULL, y = "Cells per donor", subtitle = "Sampling balance") +
    theme_ajr() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
          legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

# --- assemble --------------------------------------------------------------

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid,
                c(list(panel_a(), panel_b(), ncol = 1, labels = c("A", "B")), lab_args))
row2 <- do.call(cowplot::plot_grid,
                c(list(panel_c(), panel_d(), ncol = 2, labels = c("C", "D")), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(1.1, 0.9))
save_fig(fig, "FigureS2", dir = out_dir("supplementary"), height_mm = 208)
msg("FigureS2 written")
