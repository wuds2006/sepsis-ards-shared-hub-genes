# ---------------------------------------------------------------------------
# 39_supp_S19.R -- Figure S19. Protein-level validation of VNN2
#   A  plasma concentration of VNN1 and VNN2 in the Human Protein Atlas
#   B  VNN2 protein change across ten infectious diseases
# Input  : data/frozen/10_HPA_validation
# Output : Figures/supplementary/FigureS19.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

HP <- function(...) file.path("data", "frozen", "10_HPA_validation", ...)

# --- panel A ---------------------------------------------------------------

panel_a <- function() {
  ms <- read.delim(HP("blood_ms_concentration.tsv"), stringsAsFactors = FALSE)
  names(ms) <- c("gene", "ensg", "conc")
  ms <- ms[ms$gene %in% c("VNN1", "VNN2"), ]

  ggplot(ms, aes(conc, reorder(gene, conc), fill = gene)) +
    geom_col(width = 0.55, colour = NA) +
    geom_text(aes(label = formatC(conc, format = "e", digits = 1)), hjust = -0.15,
              family = FONT_BODY, size = 7 / .pt, colour = pal[["ink"]]) +
    scale_fill_manual(values = c(VNN1 = pal[["steel"]], VNN2 = pal[["brick"]]),
                      guide = "none") +
    scale_x_log10(expand = expansion(mult = c(0, 0.28))) +
    scale_y_discrete(labels = function(v) v) +
    labs(x = "Plasma concentration (pg/L, mass spectrometry)", y = NULL) +
    theme_ajr() +
    theme(axis.text.y = element_text(face = "italic"))
}

# --- panel B ---------------------------------------------------------------

panel_b <- function() {
  d <- read.csv(HP("HPA_CD177_VNN2_infection_vs_healthy.csv"), stringsAsFactors = FALSE)
  d <- d[d$Gene == "VNN2", ]
  d$Disease <- factor(d$Disease, levels = d$Disease[order(d$logFC)])

  ggplot(d, aes(logFC, Disease)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    geom_segment(aes(x = 0, xend = logFC, yend = Disease), linewidth = 0.35,
                 colour = "#C4C9CE") +
    geom_point(aes(size = neglogp), shape = 16, colour = pal[["brick"]]) +
    geom_text(aes(label = sig), hjust = -0.55, family = FONT_BODY,
              size = 7 / .pt, colour = pal[["ink"]]) +
    scale_size_area(max_size = 2.2, breaks = c(1, 3),
                    name = "-log10 adj. P") +
    scale_y_discrete(labels = function(v) stringr::str_wrap(v, 26)) +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.16))) +
    labs(x = "VNN2 protein log2 fold change (infection vs healthy)", y = NULL) +
    theme_ajr() +
    theme(axis.text.y = element_text(lineheight = 0.85),
          legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"))
}

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid, c(list(panel_a(), labels = "A"), lab_args))
row2 <- do.call(cowplot::plot_grid, c(list(panel_b(), labels = "B"), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(0.55, 1.45))
save_fig(fig, "FigureS19", dir = out_dir("supplementary"), height_mm = 200)
msg("FigureS19 written")
