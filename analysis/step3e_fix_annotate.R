# ============================================================
# Step 3e-fix —— 修正注释列名映射 + 重出全部 GSE167363 图
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(ragg)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "08_singlecell/GSE167363")

save_fig3 <- function(p, name, w = 7, h = 6, dir = ".", line_art = FALSE) {
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf"))
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300))
  tdpi <- if (line_art) 1000 else 500
  agg_tiff(file.path(dir, paste0(name, ".tiff")), width = w, height = h,
           units = "in", res = tdpi, compression = "lzw")
  print(p); invisible(dev.off())
  invisible(p)
}

HUB <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")

merged <- readRDS(file.path(OUT, "GSE167363_seurat_annotated.rds"))
cat("细胞数:", ncol(merged), "| 聚类数:", length(unique(merged$seurat_clusters)), "\n")

# ---------- 修正注释 ----------
markers <- list(
  `Classical Mono` = c("CD14", "S100A8", "S100A9", "VCAN", "FCN1"),
  `Intermediate Mono` = c("CD14", "FCGR3A", "MS4A7", "COTL1", "LST1"),
  `Non-classical Mono` = c("FCGR3A", "CX3CR1", "CSF1R", "MS4A7", "RGS10"),
  `Mono-DC` = intersect(c("FCER1A", "CLEC10A", "CD1C"), rownames(merged)),
  `pDC` = c("LILRA4", "GZMB", "IRF7", "IL3RA"),
  `Cycling` = c("MKI67", "TOP2A")
)
markers <- markers[sapply(markers, length) >= 2]

# 清掉旧的 ms_ 列
old <- grep("^ms_", colnames(merged@meta.data), value = TRUE)
if (length(old)) merged[[old]] <- NULL

for (ct in names(markers)) {
  merged <- AddModuleScore(merged, features = list(markers[[ct]]),
                           name = paste0("ms_", gsub("[^A-Za-z]", "", ct)))
}
score_cols <- grep("^ms_", colnames(merged@meta.data), value = TRUE)
score_mat <- FetchData(merged, vars = score_cols)
# 去掉 AddModuleScore 的尾缀数字
cn <- sub("^ms_", "", colnames(score_mat))
cn <- sub("\\d+$", "", cn)
colnames(score_mat) <- cn
merged$celltype <- colnames(score_mat)[apply(score_mat, 1, which.max)]

# 各 cluster 与 celltype 交叉核对（防止单一 cluster 被错误覆盖）
cat("\ncluster x celltype:\n"); print(table(merged$seurat_clusters, merged$celltype))
cat("\n细胞类型分布:\n"); print(table(merged$celltype))

write.csv(as.data.frame(merged[[]]), file.path(OUT, "GSE167363_cell_annot.csv"))
saveRDS(merged, file.path(OUT, "GSE167363_seurat_annotated.rds"))

# ---------- 图 ----------
CT_COLS <- c(
  "ClassicalMono"      = "#6A5C95",
  "IntermediateMono"   = "#D99AAE",
  "NonclassicalMono"   = "#E7A889",
  "MonoDC"             = "#A9CBB8",
  "pDC"                = "#CFE3E6",
  "Cycling"            = "#8FA2D6"
)
CT_COLS <- CT_COLS[intersect(names(CT_COLS), unique(merged$celltype))]

p1 <- DimPlot(merged, group.by = "celltype", label = TRUE, repel = TRUE) +
  scale_color_manual(values = CT_COLS) +
  ggtitle("GSE167363 CD14+ Monocyte Subsets") + theme_pub()
p2 <- DimPlot(merged, group.by = "group") +
  scale_color_manual(values = c(HC = "#8FA2D6", Sepsis = "#D99AAE")) +
  ggtitle("HC vs Sepsis") + theme_pub()
p3 <- DimPlot(merged, group.by = "sample") + ggtitle("Samples") + theme_pub()
save_fig3(p1 + p2 + p3, "GSE167363_UMAP", w = 16, h = 5.5, dir = OUT)

prop <- table(merged$group, merged$celltype)
prop <- prop / rowSums(prop)
write.csv(prop, file.path(OUT, "GSE167363_celltype_prop.csv"))
prop_df <- as.data.frame(prop.table(table(merged$group, merged$celltype), 1))
colnames(prop_df) <- c("Group", "CellType", "Prop")
p4 <- ggplot(prop_df, aes(x = CellType, y = Prop, fill = Group)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = c(HC = "#8FA2D6", Sepsis = "#D99AAE")) +
  ylab("Cell Proportion") + xlab(NULL) +
  ggtitle("Cell Type Proportion: HC vs Sepsis") + theme_pub() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig3(p4, "GSE167363_celltype_prop", w = 8, h = 5, dir = OUT)

# 标记基因 dotplot（注释依据）
mk_genes <- unique(unlist(markers))
mk_genes <- intersect(mk_genes, rownames(merged))
pd <- DotPlot(merged, features = mk_genes, group.by = "celltype") +
  scale_color_gradientn(colors = c("#CFE3E6", "#6A5C95", "#5B608C")) +
  ggtitle("Marker Genes by Cell Type") + theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig3(pd, "GSE167363_marker_dotplot", w = 9, h = 4.5, dir = OUT)

# ---------- hub 基因 ----------
hub_avail <- intersect(HUB, rownames(merged))
cnt <- GetAssayData(merged, layer = "counts")
expr_pct <- sapply(hub_avail, function(g) mean(cnt[g, ] > 0))
cat("\nGSE167363 可检测 hub 基因:", paste(hub_avail, collapse = ", "), "\n")
cat("表达细胞占比:\n"); print(round(expr_pct, 3))
write.csv(data.frame(gene = hub_avail, pct_expressing = round(expr_pct, 4)),
          file.path(OUT, "GSE167363_hub_expression_pct.csv"), row.names = FALSE)

pd2 <- DotPlot(merged, features = hub_avail, group.by = "celltype") +
  scale_color_gradientn(colors = c("#CFE3E6", "#6A5C95", "#5B608C")) +
  ggtitle("Hub Gene Expression by Cell Type") + theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig3(pd2, "GSE167363_hub_dotplot", w = 8, h = 4.5, dir = OUT)

pf <- FeaturePlot(merged, features = hub_avail, reduction = "umap", ncol = 3, order = TRUE) &
  scale_color_gradientn(colors = c("#CFE3E6", "#8FA2D6", "#6A5C95", "#5B608C")) &
  theme_pub()
save_fig3(pf, "GSE167363_hub_featureplot", w = 14, h = 9, dir = OUT)

pv <- VlnPlot(merged, features = hub_avail, group.by = "celltype", split.by = "group",
              split.plot = TRUE, ncol = 3, pt.size = 0) &
  scale_fill_manual(values = c(HC = "#8FA2D6", Sepsis = "#D99AAE")) &
  theme_pub()
save_fig3(pv, "GSE167363_hub_violin_group", w = 14, h = 9, dir = OUT)

gv <- FetchData(merged, vars = c("celltype", "group", hub_avail))
agg <- aggregate(gv[hub_avail], by = list(celltype = gv$celltype, group = gv$group), mean)
write.csv(agg, file.path(OUT, "GSE167363_hub_celltype_mean.csv"), row.names = FALSE)

cat("\n===== Step3e-fix 完成 =====\n")
