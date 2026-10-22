# ============================================================
# Step 2c —— Hub 基因确定 + 免疫浸润分析 (ssGSEA)
# Hub genes = ML关键基因 ∩ 双病种外验强基因
# 免疫浸润: ssGSEA (GSVA) 免疫细胞打分 + Hub基因-免疫细胞相关性
# ============================================================
suppressPackageStartupMessages({
  library(GSVA)
  library(ggplot2)
  library(pheatmap)
  library(patchwork)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "07_hub_immune")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- 确定 Hub genes ----
ml_genes <- readLines(file.path(BASE, "05_ML_SHAP", "ML_key_genes.txt"))
# 双病种外验 AUC > 0.9 的基因(强验证)
val <- read.csv(file.path(BASE, "06_external_validation", "validation_AUC_combined.csv"))
strong <- val$gene[val$type == "single" & val$AUC > 0.9]
strong <- unique(strong)
hub_genes <- intersect(ml_genes, strong)
cat("Hub genes (ML ∩ 双病种外验AUC>0.9):", paste(hub_genes, collapse=", "), "\n")
writeLines(hub_genes, file.path(OUT, "hub_genes.txt"))

# ---- 免疫细胞基因集(Charoentong 2017 经典标记) ----
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
  MDSC = c("CD33","FUT4","ITGAM","S100A9")
)

run_immune <- function(ds, disease) {
  expr <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_gene_expr.rds")))
  grp  <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_grp.rds")))
  # ssGSEA
  sets <- lapply(immune_sets, function(g) intersect(g, rownames(expr)))
  sets <- sets[lengths(sets) >= 2]
  gsva_res <- gsva(ssgseaParam(exprData = as.matrix(expr), geneSets = sets, minSize = 2))
  gsva_res <- as.data.frame(t(gsva_res))
  gsva_res$group <- grp$group
  # 组间差异(免疫细胞)
  diff_res <- do.call(rbind, lapply(colnames(gsva_res)[colnames(gsva_res)!="group"], function(cell) {
    x <- gsva_res[[cell]][gsva_res$group == disease]
    y <- gsva_res[[cell]][gsva_res$group == "Control"]
    p <- tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA)
    data.frame(cell = cell, disease_mean = mean(x), control_mean = mean(y),
               log2FC = mean(x) - mean(y), p = p)
  }))
  diff_res <- diff_res[order(diff_res$p), ]
  write.csv(diff_res, file.path(OUT, paste0(ds, "_immune_diff.csv")), row.names = FALSE)
  cat(sprintf("\n[%s] 显著免疫细胞(Top5):\n", ds))
  print(head(diff_res, 5)[, c("cell","disease_mean","control_mean","p")])

  # Hub 基因 - 免疫细胞相关性
  hub_avail <- intersect(hub_genes, rownames(expr))
  cor_mat <- cor(t(expr[hub_avail, , drop = FALSE]),
                 gsva_res[, colnames(gsva_res) != "group"], use = "pairwise.complete.obs")
  rownames(cor_mat) <- hub_avail
  write.csv(cor_mat, file.path(OUT, paste0(ds, "_hub_immune_cor.csv")))

  # 热图
  p <- pheatmap(cor_mat, cluster_cols = TRUE, cluster_rows = TRUE,
                color = colorRampPalette(c(PAL_DOWN, "white", PAL_UP))(100),
                main = paste(ds, "Hub-Immune Correlation"),
                fontsize = 9, silent = TRUE)
  ggsave(file.path(OUT, paste0(ds, "_hub_immune_heatmap.pdf")), p, width = 10, height = 5)
  ggsave(file.path(OUT, paste0(ds, "_hub_immune_heatmap.png")), p, width = 10, height = 5, dpi = 300)

  # 免疫细胞箱线图(显著 top)
  top_cells <- head(diff_res$cell, 4)
  plt <- lapply(top_cells, function(cell) {
    ggplot(gsva_res, aes(x = group, y = .data[[cell]], fill = group)) +
      geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, size = 0.8) +
      scale_fill_manual(values = PAL_GROUP) +
      labs(y = "ssGSEA score", title = cell) + theme_pub() + theme(legend.position = "none")
  })
  save_fig(wrap_plots(plt, ncol = 2), paste0(ds, "_immune_boxplot"), w = 9, h = 7, dir = OUT)
  list(gsva = gsva_res, diff = diff_res, cor = cor_mat)
}

r1 <- run_immune("GSE76293", "ARDS")
r2 <- run_immune("GSE69528", "Sepsis")

cat("\n===== Step2c 免疫浸润完成 =====\n")
cat("Hub genes:", paste(hub_genes, collapse=", "), "\n")
