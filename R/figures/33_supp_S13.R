# ---------------------------------------------------------------------------
# 33_supp_S13.R -- Figure S13. Pseudotime, programme score and regulatory
# connectivity in the CD14+ monocyte cohort
#   A  UMAP coloured by monocle3 pseudotime
#   B  AUCell six-gene signature score by subpopulation
#   C  regulatory network connectivity per subpopulation
# Input  : data/frozen/08_singlecell (slingshot, monocle3, aucell, celloracle_stratified)
# Output : Figures/supplementary/FigureS13.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

SC <- function(...) file.path("data", "frozen", "08_singlecell", ...)
CT <- c("ClassicalMono", "InflammatoryMono", "IntermediateMono", "ISGMono",
        "MonoDC", "NonclassicalMono", "pDC")

# --- panel A ---------------------------------------------------------------

panel_a <- function() {
  um <- read.csv(SC("GSE167363", "slingshot", "mono_umap.csv"), stringsAsFactors = FALSE)
  pt <- read.csv(SC("monocle3", "monocle3_pseudotime_GSE167363.csv"),
                 stringsAsFactors = FALSE)
  d <- merge(um[, c("cell", "umap_1", "umap_2")], pt[, c("cell", "pseudotime")],
             by = "cell")

  ggplot(d, aes(umap_1, umap_2, colour = pseudotime)) +
    geom_point(size = 0.22, stroke = 0) +
    scale_colour_gradientn(colours = rev(ramp_seq), name = "Pseudotime") +
    labs(x = "UMAP 1", y = "UMAP 2") +
    coord_fixed(ratio = 1) +
    theme_ajr() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_line(linewidth = 0.3),
          legend.key.height = unit(6, "mm"), legend.key.width = unit(2.2, "mm"),
          legend.title = element_text(size = PT_BODY))
}

# --- panel B ---------------------------------------------------------------

panel_b <- function() {
  d <- read.csv(SC("GSE167363", "aucell", "GSE167363_AUCell_scores.csv"),
                stringsAsFactors = FALSE)
  d$ct <- factor(d$celltype, levels = CT)
  d$grp <- factor(ifelse(d$group == "HC", "Healthy", "Sepsis"),
                  levels = c("Healthy", "Sepsis"))

  ggplot(d, aes(ct, HubSignature, fill = grp)) +
    geom_violin(width = 0.78, linewidth = 0.25, colour = "#4A4A4A",
                scale = "width", position = position_dodge(0.82)) +
    geom_boxplot(width = 0.09, outlier.shape = NA, fill = "white",
                 linewidth = 0.25, position = position_dodge(0.82)) +
    scale_fill_manual(values = c(Healthy = pal[["sage"]], Sepsis = pal[["brick"]]),
                      name = NULL) +
    scale_x_discrete(labels = function(v) stringr::str_wrap(sub("Mono$", " Mo", v), 10)) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
    labs(x = NULL, y = "AUCell score, six-gene signature") +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1,
                                     lineheight = 0.85))
}

# --- panel C ---------------------------------------------------------------

panel_c <- function() {
  f <- function(g) {
    x <- read.csv(SC("celloracle_stratified",
                     sprintf("CellOracle_TRN_connectivity_GSE167363_%s.csv", g)),
                  stringsAsFactors = FALSE)
    data.frame(ct = x$celltype, group = g, mean_abs = x$mean_abs_weight,
               n_tf = x$n_tf, n_targets = x$n_targets)
  }
  d <- rbind(f("HC"), f("Sepsis"))
  d$ct <- factor(d$ct, levels = CT)
  d$grp <- factor(ifelse(d$group == "HC", "Healthy", "Sepsis"),
                  levels = c("Healthy", "Sepsis"))

  ggplot(d, aes(ct, mean_abs, fill = grp)) +
    geom_col(position = position_dodge(width = 0.68), width = 0.62, colour = NA) +
    scale_fill_manual(values = c(Healthy = pal[["sage"]], Sepsis = pal[["brick"]]),
                      name = NULL) +
    scale_x_discrete(labels = function(v) stringr::str_wrap(sub("Mono$", " Mo", v), 10)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(x = NULL, y = "Mean absolute regulatory weight") +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1,
                                     lineheight = 0.85))
}

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid,
                c(list(panel_a(), panel_b(), ncol = 2, rel_widths = c(0.5, 0.5),
                       labels = c("A", "B")), lab_args))
row2 <- do.call(cowplot::plot_grid, c(list(panel_c(), labels = "C"), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(1.15, 0.85))
save_fig(fig, "FigureS13", dir = out_dir("supplementary"), height_mm = 206)
msg("FigureS13 written")
