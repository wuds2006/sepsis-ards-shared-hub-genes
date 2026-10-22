# ============================================================
# Step 3e —— GSE167363 (Sepsis CD14+ monocyte scRNA-seq)
# QC + 整合 + 聚类 + 单核细胞亚群注释 + hub 基因表达
# 5 样本: HC1/HC2 (健康对照) + S2_T0/S2_T6/S3_T0 (脓毒症)
# 输出: PDF + PNG(300dpi) + TIFF(ragg, 500dpi)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(Matrix)
  library(harmony)
  library(ragg)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
RAW  <- file.path(BASE, "00_rawdata/GSE167363")
OUT  <- file.path(BASE, "08_singlecell/GSE167363")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

# 统一三格式输出: PDF + PNG(300) + TIFF(ragg 500; line_art=1000)
save_fig3 <- function(p, name, w = 7, h = 6, dir = ".", line_art = FALSE) {
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf"))
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300))
  tdpi <- if (line_art) 1000 else 500
  agg_tiff(file.path(dir, paste0(name, ".tiff")), width = w, height = h,
           units = "in", res = tdpi, compression = "LZW")
  print(p); invisible(dev.off())
  invisible(p)
}

HUB <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")

# ---------- 1. 读取 + QC ----------
samples <- c("GSM5102900_HC1", "GSM5102901_HC2",
             "GSM5511353_S2_T0", "GSM5511354_S2_T6", "GSM5511355_S3_T0")
groups <- c("HC", "HC", "Sepsis", "Sepsis", "Sepsis")
timepoints <- c("HC", "HC", "T0", "T6", "T0")

objs <- list()
for (i in seq_along(samples)) {
  cat("读取", samples[i], "...\n")
  m <- Read10X(file.path(RAW, samples[i]))
  o <- CreateSeuratObject(counts = m, project = samples[i], min.cells = 3, min.features = 200)
  # 每样本下采样至 2500 细胞（控制内存）
  if (ncol(o) > 2500) {
    set.seed(42); keep <- sort(sample(colnames(o), 2500))
    o <- subset(o, cells = keep)
  }
  o$sample <- samples[i]
  o$group <- groups[i]
  o$timepoint <- timepoints[i]
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
  objs[[samples[i]]] <- o
  rm(m); gc()
}

for (nm in names(objs)) {
  objs[[nm]] <- subset(objs[[nm]], subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                         nCount_RNA < 40000 & percent.mt < 20)
}
cat("\n各样本细胞数(QC 后):\n"); print(sapply(objs, ncol))

# ---------- 2. 整合 + 聚类 ----------
merged <- merge(objs[[1]], objs[-1])
rm(objs); gc()
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, nfeatures = 2000)
merged <- ScaleData(merged)
merged <- RunPCA(merged, npcs = 30)
merged <- RunHarmony(merged, group.by.vars = "sample", reduction = "pca", reduction.save = "harmony")
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:20)
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:20)
merged <- FindClusters(merged, resolution = 0.4)
saveRDS(merged, file.path(OUT, "GSE167363_seurat_annotated.rds"))

# ---------- 3. 单核细胞亚群注释 ----------
# CD14+ monocyte 富集数据：经典/中间/非经典单核 + moDC + pDC + 少量混杂
markers <- list(
  `Classical Mono` = c("CD14", "S100A8", "S100A9", "VCAN", "FCN1"),
  `Intermediate Mono` = c("CD14", "FCGR3A", "MS4A7", "COTL1", "LST1"),
  `Non-classical Mono` = c("FCGR3A", "CX3CR1", "CSF1R", "MS4A7", "RGS10"),
  `Mono-DC` = c("FCER1A", "CLEC10A", "CD1C"),
  `pDC` = c("LILRA4", "GZMB", "IRF7", "IL3RA"),
  `Cycling` = c("MKI67", "TOP2A")
)
markers <- lapply(markers, function(g) intersect(g, rownames(merged)))
for (ct in names(markers)) {
  if (length(markers[[ct]]) >= 2) {
    merged <- AddModuleScore(merged, features = list(markers[[ct]]),
                             name = paste0("ms_", gsub("[^A-Za-z]", "", ct)))
  }
}
score_cols <- grep("^ms_", colnames(merged@meta.data), value = TRUE)
score_mat <- FetchData(merged, vars = score_cols)
lab <- sub("^ms_", "", colnames(score_mat))
lab <- c("Classical Mono"="ClassicalMono","Intermediate Mono"="IntermediateMono",
         "Non-classical Mono"="NonclassicalMono","Mono-DC"="MonoDC","pDC"="pDC","Cycling"="Cycling")
colnames(score_mat) <- unname(lab[colnames(score_mat)])
merged$celltype <- colnames(score_mat)[apply(score_mat, 1, which.max)]
cat("\n细胞类型分布:\n"); print(table(merged$celltype))

write.csv(as.data.frame(merged[[]]), file.path(OUT, "GSE167363_cell_annot.csv"))
saveRDS(merged, file.path(OUT, "GSE167363_seurat_annotated.rds"))

# ---------- 4. 图 ----------
CT_COLS <- c(
  "ClassicalMono"      = "#6A5C95",
  "IntermediateMono"   = "#D99AAE",
  "NonclassicalMono"   = "#E7A889",
  "MonoDC"             = "#A9CBB8",
  "pDC"                = "#CFE3E6",
  "Cycling"            = "#8FA2D6"
)
p1 <- DimPlot(merged, group.by = "celltype", label = TRUE, repel = TRUE) +
  scale_color_manual(values = CT_COLS) +
  ggtitle("GSE167363 CD14+ Monocyte Subsets") + theme_pub()
p2 <- DimPlot(merged, group.by = "group") +
  scale_color_manual(values = c(HC = "#8FA2D6", Sepsis = "#D99AAE")) +
  ggtitle("HC vs Sepsis") + theme_pub()
p3 <- DimPlot(merged, group.by = "sample") + ggtitle("Samples") + theme_pub()
save_fig3(p1 + p2 + p3, "GSE167363_UMAP", w = 16, h = 5.5, dir = OUT)

# 细胞类型比例（组间）
prop <- table(merged$group, merged$celltype)
prop <- prop / rowSums(prop)
write.csv(prop, file.path(OUT, "GSE167363_celltype_prop.csv"))
prop_df <- as.data.frame(prop.table(table(merged$group, merged$celltype), 1))
colnames(prop_df) <- c("Group", "CellType", "Prop")
p4 <- ggplot(prop_df, aes(x = CellType, y = Prop, fill = Group)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = c(HC = "#8FA2D6", Sepsis = "#D99AAE")) +
  ylab("Cell Proportion") + xlab(NULL) +
  ggtitle("Cell Type Proportion: HC vs Sepsis") + theme_pub() + RotatedAxis()
save_fig3(p4, "GSE167363_celltype_prop", w = 8, h = 5, dir = OUT)

# ---------- 5. hub 基因表达 ----------
hub_avail <- intersect(HUB, rownames(merged))
cat("\nGSE167363 可检测 hub 基因:", paste(hub_avail, collapse = ", "), "\n")
expr_pct <- sapply(hub_avail, function(g) mean(GetAssayData(merged, layer="counts")[g, ] > 0))
cat("表达细胞占比:\n"); print(round(expr_pct, 3))
write.csv(data.frame(gene = hub_avail, pct_expressing = round(expr_pct, 4)),
          file.path(OUT, "GSE167363_hub_expression_pct.csv"), row.names = FALSE)

if (length(hub_avail) > 0) {
  pd <- DotPlot(merged, features = hub_avail, group.by = "celltype") +
    scale_color_gradientn(colors = c("#CFE3E6", "#6A5C95", "#5B608C")) +
    ggtitle("Hub Gene Expression by Cell Type") + theme_pub() + RotatedAxis()
  save_fig3(pd, "GSE167363_hub_dotplot", w = 8, h = 4.5, dir = OUT)

  pf <- FeaturePlot(merged, features = hub_avail, reduction = "umap", ncol = 3,
                    order = TRUE) &
    scale_color_gradientn(colors = c("#CFE3E6", "#8FA2D6", "#6A5C95", "#5B608C")) &
    theme_pub()
  save_fig3(pf, "GSE167363_hub_featureplot", w = 14, h = 9, dir = OUT)

  # hub 基因 分组 小提琴图（按细胞类型拆分）
  pv <- VlnPlot(merged, features = hub_avail, group.by = "celltype", split.by = "group",
                split.plot = TRUE, ncol = 3, pt.size = 0) &
    scale_fill_manual(values = c(HC = "#8FA2D6", Sepsis = "#D99AAE")) &
    theme_pub()
  save_fig3(pv, "GSE167363_hub_violin_group", w = 14, h = 9, dir = OUT)

  # 数值汇总
  gv <- FetchData(merged, vars = c("celltype", "group", hub_avail))
  agg <- aggregate(gv[hub_avail], by = list(celltype = gv$celltype, group = gv$group), mean)
  write.csv(agg, file.path(OUT, "GSE167363_hub_celltype_mean.csv"), row.names = FALSE)
}

cat("\n===== Step3e (GSE167363 QC+注释) 完成 =====\n")
cat("细胞数:", ncol(merged), "| 细胞类型:", length(unique(merged$celltype)), "\n")
