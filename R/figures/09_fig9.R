# ---------------------------------------------------------------------------
# 09_fig9.R -- Figure 9, single-cell localisation of the hub genes
#   A  monocyte UMAP, coloured by cell type          (GSE167363, 12,359 cells)
#   B  cell-type composition, healthy versus sepsis  (GSE167363)
#   C  monocle3 pseudotime by cell type and group    (GSE167363)
#   D  AUCell six-gene signature score by cell type  (GSE151263, ARDS vs sepsis)
#   E  CellChat pathway strength shift               (GSE167363)
#   F  CellOracle in silico transcription-factor knockout
# Input  : data/frozen/08_singlecell
# Output : Figures/main/Figure9.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

SC <- function(...) file.path("data", "frozen", "08_singlecell", ...)

CT7 <- c("ClassicalMono", "InflammatoryMono", "IntermediateMono", "ISGMono",
         "MonoDC", "NonclassicalMono", "pDC")
CT7_LAB <- c("Classical monocyte", "Inflammatory monocyte", "Intermediate monocyte",
             "ISG monocyte", "Mono-DC", "Non-classical monocyte", "pDC")
ct_col <- setNames(c(pal[["brick"]], pal[["steel"]], pal[["sage"]], pal[["ochre"]],
                     pal[["plum"]], pal[["olive"]], pal[["slate"]]), CT7)

# --- panel A ---------------------------------------------------------------

panel_a <- function() {
  um <- read.csv(SC("GSE167363", "input", "umap.csv"), row.names = 1,
                 stringsAsFactors = FALSE)
  mt <- read.csv(SC("GSE167363", "input", "meta.csv"), row.names = 1,
                 stringsAsFactors = FALSE)
  d <- cbind(um[rownames(mt), ], celltype = mt$celltype)
  d$celltype <- factor(d$celltype, levels = CT7)
  cent <- aggregate(cbind(umap_1, umap_2) ~ celltype, data = d, FUN = median)
  cent$lab <- CT7_LAB[match(as.character(cent$celltype), CT7)]

  ggplot(d, aes(umap_1, umap_2, colour = celltype)) +
    geom_point(size = 0.16, stroke = 0) +
    scale_colour_manual(values = ct_col, labels = CT7_LAB, name = NULL) +
    # The seven population names are long, and a right-hand key is wider than
    # this half-width panel: it used to spill into panel B and clip its own
    # labels. The key therefore sits under the panel in two columns, and the
    # in-plot callouts are gone because their white boxes covered the cells.
    guides(colour = guide_legend(override.aes = list(size = 1.1), ncol = 2,
                                 byrow = TRUE)) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    coord_fixed(ratio = 1) +
    theme_ajr() +
    theme(legend.position = "bottom", legend.key.size = unit(2.4, "mm"),
          legend.spacing.y = unit(0.2, "mm"), legend.spacing.x = unit(1.2, "mm"),
          legend.text = element_text(size = PT_BODY),
          legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(1, 0, 0, 0),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_line(linewidth = 0.3))
}

# --- panel B ---------------------------------------------------------------

panel_b <- function() {
  p <- read.csv(SC("GSE167363", "GSE167363_celltype_prop.csv"),
                row.names = 1, stringsAsFactors = FALSE)
  d <- data.frame(celltype = factor(rep(CT7, 2), levels = rev(CT7)),
                  group = rep(c("Healthy", "Sepsis"), each = length(CT7)),
                  prop = c(as.numeric(p["HC", CT7]), as.numeric(p["Sepsis", CT7])))
  d$lab <- CT7_LAB[match(d$celltype, CT7)]

  ggplot(d, aes(prop, celltype, fill = group)) +
    geom_col(position = position_dodge(width = 0.66), width = 0.6,
             colour = NA) +
    scale_fill_manual(values = c(Healthy = pal[["sage"]], Sepsis = pal[["brick"]]),
                      name = NULL) +
    scale_y_discrete(labels = setNames(CT7_LAB, CT7)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Fraction of cells", y = NULL) +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0),
          legend.key.size = unit(2.6, "mm"))
}

# --- panel C ---------------------------------------------------------------

panel_c <- function() {
  d <- read.csv(SC("monocle3", "monocle3_pseudotime_GSE167363.csv"),
                stringsAsFactors = FALSE)
  d$celltype <- factor(d$celltype, levels = CT7)
  d$grp <- factor(ifelse(d$group == "HC", "Healthy", "Sepsis"),
                  levels = c("Healthy", "Sepsis"))
  pw <- read.csv(SC("monocle3", "monocle3_pseudotime_wilcox_GSE167363.csv"),
                 stringsAsFactors = FALSE)
  plab <- sprintf("P = %s", formatC(pw$p.value[1], format = "e", digits = 1))

  ggplot(d, aes(celltype, pseudotime, fill = grp)) +
    geom_violin(width = 0.78, linewidth = 0.25, colour = "#4A4A4A",
                scale = "width", position = position_dodge(0.82)) +
    geom_boxplot(width = 0.09, outlier.shape = NA, fill = "white",
                 linewidth = 0.25, position = position_dodge(0.82)) +
    annotate("text", x = 0.55, y = max(d$pseudotime) * 1.06, label = plab,
             hjust = 0, family = FONT_BODY, size = 7 / .pt, colour = pal[["ink"]]) +
    scale_fill_manual(values = c(Healthy = pal[["sage"]], Sepsis = pal[["brick"]]),
                      name = NULL) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(CT7_LAB[match(x, CT7)], 12)) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.10))) +
    labs(x = NULL, y = "monocle3 pseudotime") +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0),
          axis.text.x = element_text(angle = 32, hjust = 1, vjust = 1,
                                     lineheight = 0.86),
          legend.key.size = unit(2.6, "mm"))
}

# --- panel D ---------------------------------------------------------------

panel_d <- function() {
  d <- read.csv(SC("aucell", "GSE151263_AUCell_scores.csv"), stringsAsFactors = FALSE)
  d$ct <- factor(d$celltype,
                 levels = c("Mono", "CD4_T", "CD8_T", "NK", "B", "DC", "Platelet"))
  st <- read.csv(SC("aucell", "GSE151263_AUCell_stats.csv"), stringsAsFactors = FALSE)
  st <- st[st$GeneSet == "HubSignature", ]
  plab <- data.frame(ct = factor(st$celltype, levels = levels(d$ct)),
                     lab = ifelse(st$padj < 0.001, "***",
                                  ifelse(st$padj < 0.01, "**",
                                         ifelse(st$padj < 0.05, "*", ""))))

  ggplot(d, aes(ct, HubSignature, fill = group)) +
    geom_violin(width = 0.78, linewidth = 0.25, colour = "#4A4A4A",
                scale = "width", position = position_dodge(0.82)) +
    geom_boxplot(width = 0.09, outlier.shape = NA, fill = "white",
                 linewidth = 0.25, position = position_dodge(0.82)) +
    geom_text(data = plab, aes(x = ct, y = max(d$HubSignature) * 1.04, label = lab),
              inherit.aes = FALSE, family = FONT_BODY, size = 7 / .pt,
              colour = pal[["ink"]]) +
    scale_fill_manual(values = c(ARDS = pal[["steel"]], Sepsis = pal[["brick"]]),
                      name = NULL) +
    scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.12))) +
    labs(x = NULL, y = "AUCell score, six-gene signature") +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0),
          legend.key.size = unit(2.6, "mm"))
}

# --- panel E ---------------------------------------------------------------

panel_e <- function() {
  d <- read.csv(SC("GSE167363", "cellchat", "pathway_strength_comparison.csv"),
                stringsAsFactors = FALSE)
  d$diff <- d$Sepsis - d$HC
  d$lr <- log2((d$Sepsis + 1e-6) / (d$HC + 1e-6))
  d <- d[order(-abs(d$diff)), ][1:15, ]
  d$pw <- gsub("_", " ", d$pathway)
  d$pw <- factor(d$pw, levels = rev(d$pw[order(d$lr)]))
  d$side <- factor(ifelse(d$lr > 0, "Higher in sepsis", "Lower in sepsis"),
                   levels = c("Higher in sepsis", "Lower in sepsis"))

  ggplot(d, aes(lr, pw, fill = side)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = pal[["ink"]]) +
    geom_col(width = 0.66, colour = NA) +
    scale_fill_manual(values = c(`Higher in sepsis` = pal[["brick"]],
                                 `Lower in sepsis` = pal[["steel"]]), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
    labs(x = "CellChat pathway strength, log2 fold change (sepsis vs healthy)",
         y = NULL) +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0),
          legend.key.size = unit(2.6, "mm"))
}

# --- panel F ---------------------------------------------------------------

panel_f <- function() {
  d <- read.csv(SC("celloracle", "CellOracle_TF_KO_hub_effect.csv"),
                stringsAsFactors = FALSE)
  d$target <- factor(d$target, levels = c("UPP1", "EXOSC4"))
  ord <- d$KO_TF[d$target == "UPP1"][order(d$delta[d$target == "UPP1"])]
  d$KO_TF <- factor(d$KO_TF, levels = ord)

  ggplot(d, aes(delta, KO_TF, fill = target)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = pal[["ink"]]) +
    geom_col(width = 0.6, colour = NA) +
    facet_wrap(~ target, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = c(UPP1 = pal[["brick"]], EXOSC4 = pal[["steel"]]),
                      guide = "none") +
    scale_x_continuous(breaks = scales::pretty_breaks(3),
                       expand = expansion(mult = c(0.16, 0.16))) +
    labs(x = "Simulated expression change after knockout (delta)", y = NULL) +
    theme_ajr() +
    theme(strip.text = element_text(size = PT_BODY, hjust = 0,
                                    face = "italic"),
          panel.spacing.x = unit(4, "mm"))
}

# --- assemble --------------------------------------------------------------

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid,
                c(list(panel_a(), panel_b(), ncol = 2,
                       rel_widths = c(0.6, 0.4), labels = c("A", "B")), lab_args))
row2 <- do.call(cowplot::plot_grid,
                c(list(panel_c(), panel_d(), ncol = 2, labels = c("C", "D")), lab_args))
row3 <- do.call(cowplot::plot_grid,
                c(list(panel_e(), panel_f(), ncol = 2, labels = c("E", "F")), lab_args))

fig9 <- cowplot::plot_grid(row1, row2, row3, ncol = 1,
                           rel_heights = c(1.05, 0.95, 0.95))

save_fig(fig9, "Figure9", width_mm = FIG_W_MM, height_mm = 208)
msg("Figure9 written")
