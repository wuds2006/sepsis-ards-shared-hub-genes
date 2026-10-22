# ============================================================
# Step 3a-part2 —— 细胞注释 + 可视化 (加载已聚类 RDS)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "08_singlecell")
set.seed(42)  # AddModuleScore 对照基因抽样随机，固定种子保证注释可复现
merged <- readRDS(file.path(OUT, "GSE151263_seurat_annotated.rds"))
cat("细胞数:", ncol(merged), "\n")

cand_genes <- readLines(file.path(BASE, "05_ML_SHAP", "ML_key_genes.txt"))

# 细胞注释(经典 PBMC 标记, AddModuleScore)
markers <- list(
  `CD4_T` = c("CD3D", "CD3E", "CD4", "IL7R"),
  `CD8_T` = c("CD3D", "CD8A", "CD8B"),
  `NK` = c("NKG7", "GNLY", "NCAM1", "KLRD1"),
  `B` = c("MS4A1", "CD79A", "CD19"),
  `Mono` = c("CD14", "LYZ", "FCN1", "S100A8"),
  `DC` = c("FCER1A", "CST3", "CLEC9A"),
  `Platelet` = c("PPBP", "PF4")
)
for (ct in names(markers)) {
  merged <- AddModuleScore(merged, features = list(markers[[ct]]), name = ct)
}
score_cols <- paste0(names(markers), "1")
score_mat <- FetchData(merged, vars = score_cols)
merged$celltype <- names(markers)[apply(score_mat, 1, which.max)]
cat("\n细胞类型分布:\n"); print(table(merged$celltype))
saveRDS(merged, file.path(OUT, "GSE151263_seurat_annotated.rds"))

# UMAP(raster 省内存)
p1 <- DimPlot(merged, group.by = "celltype", label = TRUE, repel = TRUE, raster = TRUE) +
  scale_color_manual(values = unname(rep(PAL, 3))) + ggtitle("Cell Types") + theme_pub()
p2 <- DimPlot(merged, group.by = "group", raster = TRUE) +
  scale_color_manual(values = c(ARDS = "#8FA2D6", Sepsis = "#D99AAE")) + ggtitle("ARDS vs Sepsis") + theme_pub()
save_fig(p1 + p2, "GSE151263_UMAP", w = 12, h = 5.5, dir = OUT)

# 细胞类型比例(组间)
# 注: 规范口径的 GSE151263_celltype_prop.csv 由 scanpy/h5ad 管线生成(下游所有分析基准)，
#     Seurat 版单独存为 *_seurat.csv，避免覆盖规范文件
prop <- table(merged$group, merged$celltype)
prop <- prop / rowSums(prop)
write.csv(prop, file.path(OUT, "GSE151263_celltype_prop_seurat.csv"))

# 候选基因表达
avail_genes <- intersect(cand_genes, rownames(merged))
cat("单细胞可验证候选基因:", paste(avail_genes, collapse=", "), "\n")
if (length(avail_genes) > 0) {
  p4 <- DotPlot(merged, features = avail_genes, group.by = "celltype") +
    RotatedAxis() + ggtitle("Candidate Gene Expression by Cell Type") + theme_pub()
  save_fig(p4, "GSE151263_cand_gene_dotplot", w = 11, h = 5, dir = OUT)
  p5 <- FeaturePlot(merged, features = avail_genes[1:min(6, length(avail_genes))],
                    reduction = "umap", ncol = 3, raster = TRUE) & theme_pub()
  save_fig(p5, "GSE151263_cand_gene_featureplot", w = 14, h = 10, dir = OUT)
  # 组间差异(VlnPlot)
  p6 <- VlnPlot(merged, features = avail_genes, group.by = "group", pt.size = 0,
                cols = PAL_GROUP, ncol = 3) & theme_pub()
  save_fig(p6, "GSE151263_cand_gene_vln", w = 14, h = 8, dir = OUT)
}

cat("\n===== Step3a (GSE151263) 注释完成 =====\n")
