# ============================================================
# Step 3a —— 单细胞 QC + 整合 + 聚类 + 细胞注释 + 候选基因表达
# GSE151263 (ARDS+Sepsis PBMC, 7样本) + GSE167363 (Sepsis, 10x)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(Matrix)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
PAPER <- "/Users/mac/Desktop/2026年撰写的文章/2026论文撰写一Acute Respiratory Distress Syndrome and Sepsis"
OUT  <- file.path(BASE, "08_singlecell")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

# 候选基因(ML关键基因, 后续并入 hub genes)
cand_genes <- readLines(file.path(BASE, "05_ML_SHAP", "ML_key_genes.txt"))

# ---------- GSE151263 ----------
raw_dir <- file.path(PAPER, "ARDS/单细胞/GSE151263_RAW")
files <- list.files(raw_dir, pattern = "_processed_UMI.txt.gz$", full.names = TRUE)
names(files) <- sub("_processed_UMI.txt.gz", "", basename(files))
cat("GSE151263 样本文件:\n"); print(names(files))

objs <- list()
for (nm in names(files)) {
  cat("读取", nm, "...\n")
  # 用 fread 高效读取 基因x细胞 矩阵
  m <- fread(files[nm], header = TRUE, check.names = FALSE, sep = "\t", showProgress = FALSE)
  genes <- m[[1]]
  m <- as.matrix(m[, -1, with = FALSE])
  storage.mode(m) <- "double"
  # 重复基因名聚合
  if (anyDuplicated(genes)) {
    m <- rowsum(m, group = genes, reorder = FALSE)
  }
  # 转稀疏 + 下采样细胞(限2000)
  m <- Matrix(m, sparse = TRUE)
  gc()
  if (ncol(m) > 2000) {
    set.seed(42); keep <- sort(sample(1:ncol(m), 2000))
    m <- m[, keep]
  }
  group <- ifelse(grepl("ARDS", nm), "ARDS", "Sepsis")
  o <- CreateSeuratObject(counts = m, project = nm, min.cells = 3, min.features = 200)
  o$sample <- nm
  o$group <- group
  o$orig.ident <- nm
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
  objs[[nm]] <- o
  rm(m); gc()
}

# QC 过滤
for (nm in names(objs)) {
  objs[[nm]] <- subset(objs[[nm]], subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                         nCount_RNA < 30000 & percent.mt < 20)
}
cat("\n各样本细胞数(过滤后):\n")
print(sapply(objs, ncol))

# 合并 + 标准化
merged <- merge(objs[[1]], objs[-1])
rm(objs); gc()
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, nfeatures = 1500)
merged <- ScaleData(merged)
merged <- RunPCA(merged, npcs = 30)

# Harmony 批量校正(省内存)
library(harmony)
merged <- RunHarmony(merged, group.by.vars = "sample", reduction = "pca", reduction.save = "harmony")
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:20)
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:20)
merged <- FindClusters(merged, resolution = 0.5)
saveRDS(merged, file.path(OUT, "GSE151263_seurat_annotated.rds"))

# 细胞注释(经典 PBMC 标记)
markers <- list(
  `T cells` = c("CD3D", "CD3E", "IL7R"),
  `CD4 T` = c("CD4"),
  `CD8 T` = c("CD8A", "CD8B"),
  `NK cells` = c("NKG7", "GNLY", "NCAM1"),
  `B cells` = c("MS4A1", "CD79A", "CD19"),
  `Monocytes` = c("CD14", "LYZ", "FCN1"),
  `cDC` = c("FCER1A", "CST3", "CLEC9A"),
  `Neutrophils` = c("FCGR3B", "CSF3R", "S100A8"),
  `Platelets` = c("PPBP", "PF4")
)
# 用 AddModuleScore 打分注释
for (ct in names(markers)) {
  merged <- AddModuleScore(merged, features = list(markers[[ct]]), name = paste0("score_", gsub(" ", "", ct)))
}
score_cols <- paste0("score_", gsub(" ", "", names(markers)), "1")
score_mat <- FetchData(merged, vars = score_cols)
merged$celltype <- colnames(score_mat)[apply(score_mat, 1, which.max)]
merged$celltype <- gsub("score_", "", merged$celltype)
merged$celltype <- gsub("1$", "", merged$celltype)
cat("\n细胞类型分布:\n"); print(table(merged$celltype))

# UMAP 图
p1 <- DimPlot(merged, group.by = "celltype", label = TRUE, repel = TRUE) +
  scale_color_manual(values = unname(rep(PAL, 3))) + ggtitle("GSE151263 Cell Types") + theme_pub()
p2 <- DimPlot(merged, group.by = "group") +
  scale_color_manual(values = c(ARDS = "#8FA2D6", Sepsis = "#D99AAE")) + ggtitle("ARDS vs Sepsis") + theme_pub()
p3 <- DimPlot(merged, group.by = "sample") + ggtitle("Samples") + theme_pub()
save_fig(p1 + p2 + p3, "GSE151263_UMAP", w = 16, h = 5.5, dir = OUT)

# 细胞类型比例(组间)
prop <- table(merged$group, merged$celltype)
prop <- prop / rowSums(prop)
write.csv(prop, file.path(OUT, "GSE151263_celltype_prop.csv"))

# 候选基因表达
avail_genes <- intersect(cand_genes, rownames(merged))
cat("\n单细胞可验证候选基因:", paste(avail_genes, collapse=", "), "\n")
if (length(avail_genes) > 0) {
  p4 <- DotPlot(merged, features = avail_genes, group.by = "celltype") +
    RotatedAxis() + ggtitle("Candidate Gene Expression by Cell Type") + theme_pub()
  save_fig(p4, "GSE151263_cand_gene_dotplot", w = 11, h = 5, dir = OUT)
  p5 <- FeaturePlot(merged, features = avail_genes[1:min(6, length(avail_genes))],
                    reduction = "umap", ncol = 3) & theme_pub()
  save_fig(p5, "GSE151263_cand_gene_featureplot", w = 14, h = 10, dir = OUT)
}

cat("\n===== Step3a (GSE151263) 完成 =====\n")
cat("细胞数:", ncol(merged), "| 细胞类型:", length(unique(merged$celltype)), "\n")
