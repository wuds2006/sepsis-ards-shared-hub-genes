# ============================================================
# Step 3e-fix3 —— 补充炎症/ISG 单核模块，修正 GSE167363 注释
# clusters 1/2/4/5 全负分 -> 脓毒症炎症/ISG-high 单核
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
if (length(Layers(merged[["RNA"]])) > 1) merged <- JoinLayers(merged)
cnt <- GetAssayData(merged, layer = "counts")

# ---------- 扩展模块 ----------
markers <- list(
  `ClassicalMono`    = c("CD14", "S100A8", "S100A9", "VCAN", "FCN1"),
  `IntermediateMono` = c("CD14", "FCGR3A", "MS4A7", "COTL1", "LST1"),
  `NonclassicalMono` = c("FCGR3A", "CX3CR1", "CSF1R", "MS4A7", "RGS10"),
  `InflammatoryMono` = c("IL1B", "CCL3", "CCL4", "NFKBIA", "SAA1", "TNF"),
  `ISGMono`          = c("ISG15", "IFIT1", "IFITM3", "IFI6", "MX1", "IFI44L"),
  `MonoDC`           = intersect(c("FCER1A", "CLEC10A", "CD1C"), rownames(merged)),
  `pDC`              = c("LILRA4", "GZMB", "IRF7", "IL3RA"),
  `Cycling`          = c("MKI67", "TOP2A")
)
markers <- markers[sapply(markers, length) >= 2]
cat("可用模块:\n"); print(sapply(markers, paste, collapse = ","))

old <- grep("^ms_", colnames(merged@meta.data), value = TRUE)
if (length(old)) merged[[old]] <- NULL
for (ct in names(markers)) {
  merged <- AddModuleScore(merged, features = list(markers[[ct]]), name = paste0("ms_", ct))
}
score_cols <- grep("^ms_", colnames(merged@meta.data), value = TRUE)
sm <- FetchData(merged, vars = c("seurat_clusters", score_cols))
agg <- aggregate(sm[, score_cols], by = list(cluster = sm$seurat_clusters), mean)
rownames(agg) <- agg$cluster
score_tab <- as.matrix(agg[, -1])
colnames(score_tab) <- sub("^ms_", "", colnames(score_tab))
cat("\ncluster 级模块得分:\n"); print(round(score_tab, 3))

cl_label <- colnames(score_tab)[apply(score_tab, 1, which.max)]
cl_label <- sub("\\d+$", "", cl_label)   # 去掉 AddModuleScore 尾缀数字
names(cl_label) <- rownames(score_tab)
cat("\ncluster -> celltype:\n"); print(cl_label)
merged$celltype <- unname(cl_label[as.character(merged$seurat_clusters)])

cat("\n细胞类型分布:\n"); print(table(merged$celltype))
cat("\ngroup x celltype:\n"); print(table(merged$group, merged$celltype))

write.csv(as.data.frame(merged[[]]), file.path(OUT, "GSE167363_cell_annot.csv"))
write.csv(score_tab, file.path(OUT, "GSE167363_cluster_marker_scores.csv"))
saveRDS(merged, file.path(OUT, "GSE167363_seurat_annotated.rds"))

# ---------- 图 ----------
CT_COLS <- c(
  "ClassicalMono"      = "#6A5C95",
  "IntermediateMono"   = "#D99AAE",
  "NonclassicalMono"   = "#E7A889",
  "InflammatoryMono"   = "#5B608C",
  "ISGMono"            = "#8FA2D6",
  "MonoDC"             = "#A9CBB8",
  "pDC"                = "#CFE3E6",
  "Cycling"            = "#A9CBB8"
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

prop <- table(merged$group, merged$celltype); prop <- prop / rowSums(prop)
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

mk_genes <- unique(unlist(markers)); mk_genes <- intersect(mk_genes, rownames(merged))
pd <- DotPlot(merged, features = mk_genes, group.by = "celltype") +
  scale_color_gradientn(colors = c("#CFE3E6", "#6A5C95", "#5B608C")) +
  ggtitle("Marker Genes by Cell Type") + theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig3(pd, "GSE167363_marker_dotplot", w = 11, h = 4.5, dir = OUT)

# ---------- hub 基因 ----------
hub_avail <- intersect(HUB, rownames(merged))
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
agg2 <- aggregate(gv[hub_avail], by = list(celltype = gv$celltype, group = gv$group), mean)
write.csv(agg2, file.path(OUT, "GSE167363_hub_celltype_mean.csv"), row.names = FALSE)

cat("\n===== Step3e-fix3 完成 =====\n")
