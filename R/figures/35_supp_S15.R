# ---------------------------------------------------------------------------
# 35_supp_S15.R -- Figure S15. Metabolic pathway activity in both single-cell
# cohorts
#   A  disease-versus-control change per metabolic pathway
#   B  pathway activity by subpopulation
#   C  correlation between hub-gene expression and pathway activity
# Input  : data/frozen/08_singlecell/metabolism
# Output : Figures/supplementary/FigureS15.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

MB <- function(...) file.path("data", "frozen", "08_singlecell", "metabolism", ...)

# --- panel A ---------------------------------------------------------------

panel_a <- function() {
  f <- function(cd, lab) {
    x <- read.csv(MB(sprintf("metabolism_group_stats_%s.csv", cd)), stringsAsFactors = FALSE)
    x$cohort <- lab
    x
  }
  d <- rbind(f("GSE151263", "PBMC (GSE151263)"),
             f("GSE167363", "CD14+ monocytes (GSE167363)"))
  d <- d[order(abs(d$delta)), ]
  d$pathway <- factor(d$pathway, levels = unique(d$pathway))
  d$sig <- ifelse(d$p < 0.05, "*", "")

  ggplot(d, aes(delta, pathway, fill = cohort)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = pal[["ink"]]) +
    geom_col(position = position_dodge(width = 0.68), width = 0.62, colour = NA) +
    scale_fill_manual(values = c(`PBMC (GSE151263)` = pal[["steel"]],
                                 `CD14+ monocytes (GSE167363)` = pal[["brick"]]),
                      labels = c("PBMC", "CD14+ monocytes"), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.10))) +
    labs(x = "Mean activity change (disease minus control)", y = NULL) +
    theme_ajr() +
    theme(axis.text.y = element_text(lineheight = 0.85),
          legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"))
}

# --- panel B ---------------------------------------------------------------

panel_b <- function() {
  f <- function(cd, lab) {
    x <- read.csv(MB(sprintf("metabolism_celltype_scores_%s.csv", cd)), stringsAsFactors = FALSE)
    pw <- setdiff(names(x), c("group", "celltype"))
    d <- do.call(rbind, lapply(pw, function(p) {
      data.frame(celltype = x$celltype, pathway = p, score = x[[p]])
    }))
    d$cohort <- lab
    d
  }
  d <- rbind(f("GSE151263", "PBMC"), f("GSE167363", "CD14+ monocytes"))
  d$pathway <- gsub(".", " ", d$pathway, fixed = TRUE)
  pw <- unique(d$pathway)
  d$pathway <- factor(d$pathway, levels = rev(pw))

  ggplot(d, aes(celltype, pathway, fill = score)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    facet_wrap(~ cohort, ncol = 1, scales = "free_x") +
    scale_fill_gradientn(colours = rev(ramp_seq), name = "Activity") +
    scale_x_discrete(labels = function(v) stringr::str_wrap(sub("Mono$", " Mo", v), 10)) +
    scale_y_discrete(labels = function(v) stringr::str_wrap(v, 16)) +
    labs(x = NULL, y = NULL) +
    theme_ajr() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1,
                                     lineheight = 0.85),
          axis.text.y = element_text(lineheight = 0.85),
          strip.text = element_text(size = PT_BODY, hjust = 0),
          axis.line = element_blank(), axis.ticks = element_blank(),
          legend.key.height = unit(6, "mm"), legend.key.width = unit(2.2, "mm"),
          legend.title = element_text(size = PT_BODY))
}

# --- panel C ---------------------------------------------------------------

panel_c <- function() {
  # the frozen table stores "r (p=...)" strings; recover the coefficient
  x <- read.csv(MB("metabolism_hub_correlation_GSE167363.csv"), stringsAsFactors = FALSE)
  m <- as.matrix(x[, -1])
  rho <- matrix(as.numeric(sub(" .*$", "", m)), nrow = nrow(m),
                dimnames = list(x$gene, colnames(m)))
  d <- as.data.frame(as.table(rho))
  names(d) <- c("gene", "pathway", "rho")
  d$gene <- factor(d$gene, levels = rev(HUB6))
  d$pathway <- factor(d$pathway, levels = colnames(m))

  ggplot(d, aes(pathway, gene, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_gradientn(colours = ramp_div, limits = c(-0.6, 0.6),
                         oob = scales::squish, name = "Spearman rho") +
    scale_x_discrete(labels = function(v) stringr::str_wrap(v, 14)) +
    labs(x = NULL, y = NULL) +
    theme_ajr() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1,
                                     lineheight = 0.85),
          axis.text.y = element_text(face = "italic"),
          axis.line = element_blank(), axis.ticks = element_blank(),
          legend.key.height = unit(6, "mm"), legend.key.width = unit(2.2, "mm"),
          legend.title = element_text(size = PT_BODY))
}

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid,
                c(list(panel_a(), panel_b(), ncol = 2, labels = c("A", "B")), lab_args))
row2 <- do.call(cowplot::plot_grid, c(list(panel_c(), labels = "C"), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(0.95, 0.9))
save_fig(fig, "FigureS15", dir = out_dir("supplementary"), height_mm = 204)
msg("FigureS15 written")
