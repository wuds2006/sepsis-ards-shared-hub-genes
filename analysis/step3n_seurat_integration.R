# ============================================================
# Step 3n —— 模块2补做：Seurat integration (RPCA) 批次校正
# 对 GSE151263 (PBMC) 与 GSE167363 (CD14+ Mono) 各自按样本做
# Seurat RPCA anchor 整合，与 Harmony 并排对比，并量化批次混合度。
# 输出: integration/ 目录下 UMAP 对比图 + 批次混合指标 CSV。
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(ggplot2); library(patchwork); library(Matrix); library(ragg)
})
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/common.R")
set.seed(42)

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"

save_fig3 <- function(p, name, w = 7, h = 6, dir = ".", line_art = FALSE) {
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf"))
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300))
  tdpi <- if (line_art) 1000 else 500
  agg_tiff(file.path(dir, paste0(name, ".tiff")), width = w, height = h,
           units = "in", res = tdpi, compression = "lzw")
  print(p); invisible(dev.off())
  invisible(p)
}

# 批次混合度：局部逆 Simpson 指数（每细胞 15 近邻的样本多样性）
batch_mixing <- function(obj, reduction, sample_col = "sample", k = 15) {
  emb <- Embeddings(obj, reduction = reduction)
  n <- nrow(emb)
  k <- min(k, n - 1)
  nn <- FNN::get.knn(emb, k = k)$nn.index
  samples <- as.character(obj[[sample_col]][[1]])
  lisi <- vapply(seq_len(n), function(i) {
    tb <- table(samples[nn[i, ]])
    p <- as.numeric(tb) / sum(tb)
    1 / sum(p^2)   # 1 = 完全不混合, n.batch = 完全混合
  }, numeric(1))
  lisi
}

PAL7 <- unname(PAL)

run_integration <- function(rds_path, dataset, group_levels) {
  cat("\n##########", dataset, "##########\n")
  obj <- readRDS(rds_path)
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

  # ---- Harmony（既有基准，重算用于同图对比）----
  obj <- RunHarmony(obj, group.by.vars = "sample", verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:20, verbose = FALSE,
                 reduction.name = "umap.harmony", reduction.key = "umapH_")

  # ---- Seurat v5 reciprocal PCA integration：IntegrateLayers ----
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- IntegrateLayers(obj, method = RPCAIntegration, orig.reduction = "pca",
                         new.reduction = "rpca.integrated", dims = 1:20, verbose = FALSE)
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  obj <- RunUMAP(obj, reduction = "rpca.integrated", dims = 1:20, verbose = FALSE,
                 reduction.name = "umap.integrated", reduction.key = "umapI_")

  # ---- 批次混合量化 ----
  m_harmony <- batch_mixing(obj, "umap.harmony")
  m_integ   <- batch_mixing(obj, "umap.integrated")
  n_batch   <- length(unique(obj$sample))
  metric <- data.frame(
    dataset = dataset,
    method = c("Harmony", "Seurat_RPCA"),
    mean_local_inv_simpson = c(mean(m_harmony), mean(m_integ)),
    max_possible = n_batch,
    mixing_ratio = c(mean(m_harmony) / n_batch, mean(m_integ) / n_batch)
  )
  print(metric)

  # ---- 图：2x2 对比（Harmony vs Integration；按样本着色 / 按细胞类型着色）----
  p1 <- DimPlot(obj, reduction = "umap.harmony", group.by = "sample",
                cols = PAL7, pt.size = 0.3) +
    ggtitle("Harmony (by sample)") + theme_pub(10) + NoLegend() +
    theme(plot.title = element_text(size = 11))
  p2 <- DimPlot(obj, reduction = "umap.integrated", group.by = "sample",
                cols = PAL7, pt.size = 0.3) +
    ggtitle("Seurat RP integration (by sample)") + theme_pub(10) +
    theme(plot.title = element_text(size = 11))
  p3 <- DimPlot(obj, reduction = "umap.harmony", group.by = "celltype",
                pt.size = 0.3) +
    ggtitle("Harmony (by cell type)") + theme_pub(10) +
    theme(plot.title = element_text(size = 11))
  p4 <- DimPlot(obj, reduction = "umap.integrated", group.by = "celltype",
                pt.size = 0.3) +
    ggtitle("Seurat RP integration (by cell type)") + theme_pub(10) +
    theme(plot.title = element_text(size = 11))
  comp <- (p1 | p2) / (p3 | p4) + plot_annotation(
    title = paste0(dataset, ": Harmony vs Seurat integration batch correction"),
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5)))
  save_fig3(comp, paste0(dataset, "_Harmony_vs_seurat_integration"), w = 12, h = 10, dir = OUT)

  # ---- 图：批次混合度柱状图 ----
  bm <- data.frame(method = rep(c("Harmony", "Seurat_RPCA"), each = length(m_harmony)),
                   lisi = c(m_harmony, m_integ))
  bm$method <- factor(bm$method, levels = c("Harmony", "Seurat_RPCA"))
  pb <- ggplot(bm, aes(x = method, y = lisi, fill = method)) +
    geom_violin(color = "grey30", linewidth = 0.3, width = 0.8) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "grey20") +
    scale_fill_manual(values = c("#8FA2D6", "#D99AAE")) +
    geom_hline(yintercept = n_batch, linetype = "dashed", color = "#5B608C") +
    annotate("text", x = 0.55, y = n_batch, label = paste0("complete mixing = ", n_batch),
             hjust = 0, vjust = -0.3, size = 3, color = "#5B608C") +
    ggtitle(paste0(dataset, ": local inverse Simpson index (n=15 NN)")) +
    ylab("Batch mixing (local inverse Simpson)") + xlab(NULL) +
    theme_pub() + theme(legend.position = "none")
  save_fig3(pb, paste0(dataset, "_batch_mixing_metrics"), w = 5.5, h = 5, dir = OUT)

  metric
}

OUT <- file.path(BASE, "08_singlecell/integration")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# GSE151263 PBMC：3 ARDS + 4 Sepsis 样本
m1 <- tryCatch(
  run_integration(file.path(BASE, "08_singlecell/GSE151263_seurat_annotated.rds"),
                  "GSE151263", c("ARDS", "Sepsis")),
  error = function(e) { cat("GSE151263 integration 失败:", conditionMessage(e), "\n"); NULL })

# GSE167363 CD14+ Mono：2 HC + 3 Sepsis 样本
m2 <- tryCatch(
  run_integration(file.path(BASE, "08_singlecell/GSE167363/GSE167363_seurat_annotated.rds"),
                  "GSE167363", c("HC", "Sepsis")),
  error = function(e) { cat("GSE167363 integration 失败:", conditionMessage(e), "\n"); NULL })

metric_all <- do.call(rbind, Filter(Negate(is.null), list(m1, m2)))
write.csv(metric_all, file.path(OUT, "integration_batch_mixing_metrics.csv"), row.names = FALSE)
cat("\n===== Step3n Seurat integration 完成 =====\n")
print(metric_all)
