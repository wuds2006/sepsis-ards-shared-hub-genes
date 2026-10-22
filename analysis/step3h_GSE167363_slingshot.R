# ============================================================
# Step 3h —— GSE167363 slingshot 拟时序
# 单核细胞亚群分化轨迹: Classical -> Inflammatory/ISG -> Intermediate -> Non-classical
# 排除 pDC/MonoDC, 起点设为 ClassicalMono
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(slingshot)
  library(ggplot2)
  library(reshape2)
  library(ragg)
  library(patchwork)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "08_singlecell/GSE167363/slingshot")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

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

merged <- readRDS(file.path(BASE, "08_singlecell/GSE167363/GSE167363_seurat_annotated.rds"))
if (length(Layers(merged[["RNA"]])) > 1) merged <- JoinLayers(merged)

# 只保留单核细胞谱系
mono_types <- grep("Mono", unique(merged$celltype), value = TRUE)
mono <- subset(merged, subset = celltype %in% mono_types)
cat("单核谱系细胞:", ncol(mono), "| 亚群:", paste(sort(unique(mono$celltype)), collapse=", "), "\n")

# 重跑 PCA（子集）
DefaultAssay(mono) <- "RNA"
mono <- NormalizeData(mono, verbose = FALSE)
mono <- FindVariableFeatures(mono, nfeatures = 1500, verbose = FALSE)
mono <- ScaleData(mono, features = VariableFeatures(mono), verbose = FALSE)
mono <- RunPCA(mono, npcs = 30, verbose = FALSE)

pca <- Embeddings(mono, "pca")
umap_df <- data.frame(Embeddings(merged, "umap"))
umap_df$cell <- rownames(umap_df)
umap_df <- umap_df[colnames(mono), ]
meta <- mono[[]][, c("celltype", "group", "sample", "seurat_clusters")]

# hub 基因表达
hub_avail <- intersect(HUB, rownames(mono))
expr <- t(as.matrix(GetAssayData(mono, layer = "data")[hub_avail, ]))

write.csv(pca[, 1:20], file.path(OUT, "mono_pca.csv"))
write.csv(umap_df, file.path(OUT, "mono_umap.csv"))
write.csv(meta, file.path(OUT, "mono_meta.csv"))
write.csv(expr, file.path(OUT, "mono_gene_expr.csv"))

# ---------- slingshot ----------
rd <- pca[, 1:20]
cl <- as.character(meta$celltype)
start_cl <- "ClassicalMono"
if (!(start_cl %in% cl)) start_cl <- names(sort(table(cl), decreasing = TRUE))[1]
cat("起点亚群:", start_cl, "\n")

sds <- slingshot(rd, cl, start.clus = start_cl, approx_points = 100)
print(slingLineages(sds))

pt <- slingPseudotime(sds)
colnames(pt) <- paste0("Lineage", seq_len(ncol(pt)))
pseudo <- apply(pt, 1, min, na.rm = TRUE)

saveRDS(sds, file.path(OUT, "slingshot_sds.rds"))
write.csv(pt, file.path(OUT, "slingshot_pseudotime.csv"))
write.csv(data.frame(cell = rownames(pca), pseudotime = pseudo),
          file.path(OUT, "slingshot_main_pseudotime.csv"), row.names = FALSE)

# ---------- 可视化 ----------
CT_COLS <- c(
  "ClassicalMono"      = "#6A5C95",
  "IntermediateMono"   = "#D99AAE",
  "NonclassicalMono"   = "#E7A889",
  "InflammatoryMono"   = "#5B608C",
  "ISGMono"            = "#8FA2D6"
)
CT_COLS <- CT_COLS[intersect(names(CT_COLS), unique(cl))]

df <- data.frame(UMAP1 = umap_df[, 1], UMAP2 = umap_df[, 2],
                 celltype = cl, group = meta$group, pseudotime = pseudo)
rownames(df) <- rownames(pca)

# 轨迹曲线（在 UMAP 上）
lin <- slingCurves(sds)
# 把曲线从 PCA 空间映射到 UMAP（用线性插值近似：直接用 PCA 曲线细胞的最近邻）
# 更直接：在 PCA 空间画曲线 + UMAP 着色细胞
p1 <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = celltype)) +
  geom_point(size = 0.4, stroke = 0) +
  scale_color_manual(values = CT_COLS) +
  ggtitle("GSE167363 Monocyte Subsets (UMAP)") + theme_pub()
p2 <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = pseudotime)) +
  geom_point(size = 0.4, stroke = 0) +
  scale_color_gradientn(colors = c("#CFE3E6", "#8FA2D6", "#6A5C95", "#5B608C")) +
  ggtitle("Pseudotime (slingshot)") + theme_pub()
p3 <- ggplot(df, aes(x = pseudotime, fill = celltype)) +
  geom_density(alpha = 0.6, color = NA, position = "stack") +
  scale_fill_manual(values = CT_COLS) +
  ggtitle("Celltype Composition along Pseudotime") + theme_pub()
save_fig3(p1 + p2 + p3, "GSE167363_slingshot_overview", w = 16, h = 5, dir = OUT)

# PCA 空间轨迹（真实曲线）
pca_df <- data.frame(PC1 = pca[, 1], PC2 = pca[, 2], celltype = cl, pseudotime = pseudo)
base <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = celltype), size = 0.3, stroke = 0) +
  scale_color_manual(values = CT_COLS) + theme_pub()
for (i in seq_along(lin)) {
  lc <- as.data.frame(lin[[i]]$s[lin[[i]]$ord, , drop = FALSE])
  names(lc) <- colnames(rd)   # 与 slingshot 输入空间一致(20 维)
  base <- base + geom_path(data = lc, aes(x = PC1, y = PC2),
                           color = "#5B608C", linewidth = 1.1, inherit.aes = FALSE)
}
base <- base + ggtitle(paste("Trajectory in PCA space (", length(lin), "lineages )"))
save_fig3(base, "GSE167363_slingshot_trajectory_PCA", w = 7, h = 6, dir = OUT)

# ---------- hub 基因沿伪时间 ----------
cor_res <- data.frame()
for (g in hub_avail) {
  ct <- suppressWarnings(cor.test(pseudo, expr[, g], method = "spearman"))
  cor_res <- rbind(cor_res, data.frame(gene = g, rho = ct$estimate, p = ct$p.value))
}
cor_res$fdr <- p.adjust(cor_res$p, method = "BH")
cat("\n=== hub 基因与伪时间 Spearman 相关 ===\n"); print(cor_res)
write.csv(cor_res, file.path(OUT, "slingshot_gene_trajectory_cor.csv"), row.names = FALSE)

# 趋势（5 分箱）
pseudo_bin <- cut(pseudo, breaks = 5, labels = c("T1","T2","T3","T4","T5"))
trend <- data.frame()
for (g in hub_avail) {
  m <- tapply(expr[, g], pseudo_bin, mean, na.rm = TRUE)
  trend <- rbind(trend, data.frame(gene = g, t(m)))
}
write.csv(trend, file.path(OUT, "slingshot_gene_trend.csv"), row.names = FALSE)
cat("\n=== hub 基因沿伪时间表达趋势(5 分箱均值) ===\n"); print(trend, row.names = FALSE)

# 基因-伪时间平滑曲线图
gdf <- data.frame(pseudotime = rep(pseudo, length(hub_avail)),
                  expr = as.vector(expr),
                  gene = rep(hub_avail, each = length(pseudo)))
gdf$expr[log1p(abs(gdf$expr)) < 0] <- NA
pg <- ggplot(na.omit(gdf), aes(x = pseudotime, y = expr)) +
  geom_point(size = 0.2, color = "#CFE3E6", alpha = 0.4) +
  geom_smooth(aes(color = gene), method = "loess", se = FALSE, linewidth = 0.9,
              method.args = list(span = 0.3)) +
  scale_color_manual(values = c("#6A5C95", "#D99AAE", "#E7A889", "#A9CBB8", "#CFE3E6", "#8FA2D6")[seq_along(hub_avail)]) +
  facet_wrap(~gene, ncol = 3, scales = "free_y") +
  ggtitle("Hub Gene Expression along Pseudotime") + theme_pub()
save_fig3(pg, "GSE167363_slingshot_gene_trend", w = 12, h = 7, dir = OUT, line_art = TRUE)

cat("\n===== Step3h slingshot (GSE167363) 完成 =====\n")
