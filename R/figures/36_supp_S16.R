# ---------------------------------------------------------------------------
# 36_supp_S16.R -- Figure S16. Donor-level pseudobulk and sensitivity analyses
#   A  pseudobulk differential expression, ARDS versus sepsis
#   B  pseudobulk differential expression, sepsis versus healthy
#   C  donor-level CEBPB and UPP1 abundance
# Input  : data/frozen/08_singlecell/donor_level
# Output : Figures/supplementary/FigureS16.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

DL <- function(...) file.path("data", "frozen", "08_singlecell", "donor_level", ...)

volcano <- function(file, sub, xlab) {
  d <- read.csv(DL(file), stringsAsFactors = FALSE)
  d <- d[!is.na(d$padj) & !is.na(d$log2FoldChange), ]
  d$hub <- d$gene %in% HUB6
  d$sig <- d$padj < 0.05

  ggplot(d, aes(log2FoldChange, -log10(padj))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    geom_point(aes(colour = sig), size = 0.35, alpha = 0.55) +
    geom_point(data = d[d$hub, ], shape = 21, size = 1.4, stroke = 0.35,
               colour = pal[["ink"]], fill = pal[["brick"]]) +
    ggrepel::geom_text_repel(data = d[d$hub, ], aes(label = gene), size = 7 / .pt,
                             family = FONT_BODY, fontface = "italic",
                             colour = pal[["ink"]], min.segment.length = 0,
                             segment.size = 0.2, box.padding = 0.2,
                             max.overlaps = Inf, seed = 6) +
    scale_colour_manual(values = c(`TRUE` = pal[["brick"]], `FALSE` = "#C4C9CE"),
                        labels = c(`TRUE` = "adjusted P < 0.05", `FALSE` = "not significant"),
                        name = NULL) +
    labs(x = xlab, y = "-log10 adjusted P", subtitle = sub) +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

# --- panel C ---------------------------------------------------------------

panel_c <- function() {
  d <- read.csv(DL("GSE167363_donor_level_UPP1_CEBPB.csv"), stringsAsFactors = FALSE)
  d$grp <- factor(ifelse(d$group == "HC", "Healthy", "Sepsis"),
                  levels = c("Healthy", "Sepsis"))

  ggplot(d, aes(CEBPB_pb_log2CPM, UPP1_pb_log2CPM, colour = grp)) +
    geom_point(size = 1.6) +
    ggrepel::geom_text_repel(aes(label = donor), size = 7 / .pt, family = FONT_BODY,
                             colour = pal[["ink"]], min.segment.length = 0,
                             segment.size = 0.2, box.padding = 0.2,
                             max.overlaps = Inf, seed = 8) +
    scale_colour_manual(values = c(Healthy = pal[["sage"]], Sepsis = pal[["brick"]]),
                        name = NULL) +
    labs(x = "CEBPB pseudobulk expression (log2 CPM)",
         y = "UPP1 pseudobulk expression (log2 CPM)") +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"))
}

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid,
                c(list(volcano("pseudobulk_DESeq2_GSE151263_ARDS_vs_Sepsis.csv",
                               "PBMC, ARDS versus sepsis", "log2 fold change (ARDS vs sepsis)"),
                       volcano("pseudobulk_DESeq2_GSE167363_Sepsis_vs_HC.csv",
                               "CD14+ monocytes, sepsis versus healthy",
                               "log2 fold change (sepsis vs healthy)"),
                       ncol = 2, labels = c("A", "B")), lab_args))
row2 <- do.call(cowplot::plot_grid, c(list(panel_c(), labels = "C"), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(1.15, 0.85))
save_fig(fig, "FigureS16", dir = out_dir("supplementary"), height_mm = 206)
msg("FigureS16 written")
