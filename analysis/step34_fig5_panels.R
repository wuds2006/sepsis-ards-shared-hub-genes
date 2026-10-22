#!/usr/bin/env Rscript
# Figure 5 panels (ssGSEA immune infiltration, ggplot2 + ggsci NPG):
#   A,B: immune-cell score boxplots (GSE76293, GSE69528)
#   C,D: hub gene vs immune cell Spearman correlation heatmaps
# Per-sample ssGSEA scores are recomputed with the same gene sets as step2c.
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
suppressPackageStartupMessages({ library(GSVA); library(dplyr); library(tidyr) })

BASE <- PROJ
IMM <- file.path(BASE, "07_hub_immune")

immune_sets <- list(
  CD8_T = c("CD8A","CD8B","GZMA","GZMB"),
  CD4_T = c("CD4","IL7R","CD3D","CD3E"),
  Treg = c("FOXP3","IL2RA","CTLA4","TIGIT"),
  Th1 = c("TBX21","IFNG","STAT4","IL12RB2"),
  Th2 = c("GATA3","IL4","IL5","IL13"),
  Th17 = c("RORC","IL17A","IL22","IL23R"),
  NK = c("NKG7","KLRD1","NCAM1","GNLY"),
  B_cell = c("MS4A1","CD79A","CD19","CD79B"),
  Plasma = c("MZB1","XBP1","TNFRSF17","SDC1"),
  Monocyte = c("CD14","LYZ","FCN1","S100A12"),
  M1_Macro = c("IRF5","NOS2","TNF","IL1B"),
  M2_Macro = c("CD163","MRC1","ARG1","MSR1"),
  Neutrophil = c("FCGR3B","CSF3R","S100A8","CEACAM8"),
  DC = c("ITGAX","CD1C","FCER1A","CST3"),
  pDC = c("LILRA4","CLEC4C","IRF7","TCF4"),
  Mast = c("KIT","TPSAB1","CPA3","MS4A2"),
  Eosinophil = c("SIGLEC8","CCR3","RNASE2","IL5RA"),
  MDSC = c("CD33","FUT4","ITGAM","S100A9"))

star_fun <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

immune_box <- function(ds, dis, panel, ylab = "ssGSEA score") {
  expr <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_gene_expr.rds")))
  grp  <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_grp.rds")))
  sets <- lapply(immune_sets, function(g) intersect(g, rownames(expr)))
  sets <- sets[lengths(sets) >= 2]
  sc <- gsva(ssgseaParam(exprData = as.matrix(expr), geneSets = sets, minSize = 2))
  diff <- read.csv(file.path(IMM, paste0(ds, "_immune_diff.csv")))
  ord <- diff$cell[order(diff$log2FC, decreasing = TRUE)]
  ord <- ord[ord %in% rownames(sc)]
  long <- as.data.frame(t(sc[ord, , drop = FALSE]))
  long$sample <- rownames(long)
  long$group <- grp$group[match(long$sample, rownames(grp))]
  if (is.null(long$group) || anyNA(long$group)) long$group <- grp$group[match(long$sample, grp$sample)]
  long <- pivot_longer(long, cols = all_of(ord), names_to = "cell", values_to = "score")
  long$cell <- factor(long$cell, levels = ord)
  long$group <- factor(long$group, levels = c("Control", dis))
  sig <- diff %>% mutate(star = star_fun(p)) %>% select(cell, star)
  ymax <- long %>% group_by(cell) %>% summarise(y = max(score) * 1.02 + 0.005, .groups = "drop")
  sig <- left_join(sig, ymax, by = "cell") %>% filter(star != "ns")
  ggplot(long, aes(cell, score, fill = group)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.28, width = 0.68,
                 position = position_dodge(0.75)) +
    geom_text(data = sig, aes(x = cell, y = y, label = star), inherit.aes = FALSE,
              size = 2.5, vjust = 0.4) +
    scale_fill_manual(values = setNames(c(COL_CONTROL, COL_DISEASE), c("Control", dis))) +
    labs(title = paste0(panel, "  ", ds), x = NULL, y = ylab, fill = NULL) +
    theme_ajrccm() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
          legend.position = "top", legend.text = element_text(size = 6))
}

p5a <- immune_box("GSE76293", "ARDS", "A")
save_pdf(p5a, file.path(OUT_MAIN, "Fig5A_ssGSEA_boxplot_GSE76293.pdf"), 140, 72)
p5b <- immune_box("GSE69528", "Sepsis", "B")
save_pdf(p5b, file.path(OUT_MAIN, "Fig5B_ssGSEA_boxplot_GSE69528.pdf"), 140, 72)

cor_heat <- function(ds, panel) {
  cm <- as.matrix(read.csv(file.path(IMM, paste0(ds, "_hub_immune_cor.csv")), row.names = 1))
  rn <- rownames(cm); cn <- colnames(cm)
  long <- expand.grid(gene = seq_along(rn), cell = seq_along(cn))
  long$rho <- as.vector(cm[as.matrix(long)])
  long$gene <- factor(rn[long$gene], levels = rev(rn))
  long$cell <- factor(cn[long$cell], levels = cn)
  ggplot(long, aes(cell, gene, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#3C5488FF", mid = "white", high = "#E64B35FF",
                         midpoint = 0, limits = c(-1, 1), name = "Spearman rho") +
    scale_y_discrete(labels = function(x) parse(text = gene_italic(x))) +
    labs(title = paste0(panel, "  ", ds), x = NULL, y = NULL) +
    theme_ajrccm() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
          legend.position = "right", legend.key.height = unit(4, "mm"),
          legend.title = element_text(size = 6))
}

p5c <- cor_heat("GSE76293", "C")
save_pdf(p5c, file.path(OUT_MAIN, "Fig5C_corr_heatmap_GSE76293.pdf"), 120, 62)
p5d <- cor_heat("GSE69528", "D")
save_pdf(p5d, file.path(OUT_MAIN, "Fig5D_corr_heatmap_GSE69528.pdf"), 120, 62)

cat("Figure 5 panels done.\n")
