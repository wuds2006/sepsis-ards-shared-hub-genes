# ============================================================
# Step3k: GSE151263 (PBMC) 缺失图件补齐 — R 部分
#   （Seurat rds 已损坏，hub 图件改由 Python 从 scanpy h5ad 绘制）
#   1) CellChat 差异可视化补格式（diff_LR_bubble / rankNet / diff_heatmap）
#   2) slingshot 轨迹图（三格式）+ hub 基因拟时序趋势线图（三格式）
# 规范: PDF(矢量) + PNG(300dpi) + TIFF(ragg agg_tiff, lzw, 彩色500)
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ComplexHeatmap)
  library(CellChat)
  library(slingshot)
  library(ragg)
})

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT <- file.path(BASE, "08_singlecell")

PAL <- c("#6A5C95", "#D99AAE", "#E7A889", "#A9CBB8", "#CFE3E6", "#8FA2D6", "#5B608C")

save_fig3 <- function(p, name, w = 7, h = 6, dir = ".", line_art = FALSE) {
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf"))
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300))
  tdpi <- if (line_art) 1000 else 500
  agg_tiff(file.path(dir, paste0(name, ".tiff")), width = w, height = h,
           units = "in", res = tdpi, compression = "lzw")
  print(p); invisible(dev.off())
  invisible(p)
}

save_base3 <- function(name, plot_expr, w = 7, h = 6, dir = ".", line_art = FALSE) {
  tdpi <- if (line_art) 1000 else 500
  pdf(file.path(dir, paste0(name, ".pdf")), width = w, height = h)
  eval(plot_expr); invisible(dev.off())
  png(file.path(dir, paste0(name, ".png")), width = w, height = h, units = "in", res = 300)
  eval(plot_expr); invisible(dev.off())
  agg_tiff(file.path(dir, paste0(name, ".tiff")), width = w, height = h,
           units = "in", res = tdpi, compression = "lzw")
  eval(plot_expr); invisible(dev.off())
}

theme_pub <- function() {
  theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          axis.text = element_text(color = "black"),
          legend.position = "right")
}

# ============================================================
# 1) CellChat 差异可视化（补齐三格式）
# ============================================================
cat("[1] CellChat 差异可视化补格式 ...\n")
CC <- file.path(OUT, "cellchat")
ccm <- readRDS(file.path(CC, "CellChat_merged.rds"))

# 1a. 差异 L-R bubble
tryCatch({
  p <- netVisual_bubble(ccm, comparison = c(1, 2), angle.x = 45,
                        remove.isolate = TRUE, return.data = FALSE)
  save_fig3(p, "diff_LR_bubble", w = 13, h = 11, dir = CC)
  cat("    diff_LR_bubble 完成\n")
}, error = function(e) cat("    bubble 失败:", conditionMessage(e), "\n"))

# 1b. rankNet
tryCatch({
  gg2 <- rankNet(ccm, mode = "comparison", stacked = TRUE, do.stat = TRUE)
  save_fig3(gg2, "rankNet_pathways", w = 8, h = 7, dir = CC)
  cat("    rankNet 完成\n")
}, error = function(e) cat("    rankNet 失败:", conditionMessage(e), "\n"))

# 1c. 差异交互强度 heatmap（ComplexHeatmap）
tryCatch({
  gg1 <- netVisual_heatmap(ccm)
  pdf(file.path(CC, "diff_heatmap.pdf"), width = 8, height = 6)
  draw(gg1); dev.off()
  png(file.path(CC, "diff_heatmap.png"), width = 8, height = 6, units = "in", res = 300)
  draw(gg1); dev.off()
  agg_tiff(file.path(CC, "diff_heatmap.tiff"), width = 8, height = 6,
           units = "in", res = 500, compression = "lzw")
  draw(gg1); dev.off()
  cat("    diff_heatmap 完成\n")
}, error = function(e) cat("    heatmap 失败:", conditionMessage(e), "\n"))

# ============================================================
# 2) slingshot 轨迹图 + 基因趋势（三格式）
# ============================================================
cat("\n[2] slingshot 轨迹图补格式 ...\n")
SL <- file.path(OUT, "slingshot")
umap <- read.csv(file.path(SL, "mono_umap.csv"), row.names = 1)
meta <- read.csv(file.path(SL, "mono_meta.csv"), row.names = 1)
sds <- readRDS(file.path(SL, "slingshot_sds.rds"))

grp_cols <- c(ARDS = "#D99AAE", Sepsis = "#8FA2D6")
g <- factor(meta$group, levels = c("ARDS", "Sepsis"))
plot_traj <- function() {
  par(mar = c(4, 4, 2, 1))
  plot(umap$UMAP1, umap$UMAP2, col = grp_cols[as.character(g)], pch = 16, cex = 0.5,
       xlab = "UMAP1", ylab = "UMAP2",
       main = "Monocyte differentiation trajectory (slingshot, GSE151263 PBMC)")
  lines(SlingshotDataSet(sds), lwd = 2.5, col = "black")
  legend("topright", legend = levels(g), col = grp_cols[levels(g)], pch = 16, cex = 0.9, bty = "n")
}
save_base3("slingshot_trajectory_UMAP", quote(plot_traj()), w = 7.5, h = 6, dir = SL)
cat("    轨迹图完成（覆盖旧 PNG）\n")

# 基因沿伪时序趋势线图（z-score 归一化便于同轴展示）
trend <- read.csv(file.path(SL, "slingshot_gene_trend.csv"), row.names = 1)
genes_avail <- rownames(trend)
zn <- t(scale(t(as.matrix(trend))))
zdf <- as.data.frame(zn) %>% tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "bin", values_to = "z")
zdf$gene <- factor(zdf$gene, levels = genes_avail)
p_trend <- ggplot(zdf, aes(x = bin, y = z, color = gene, group = gene)) +
  geom_point(size = 2.2) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = rep(PAL, 3)) +
  labs(x = "Pseudotime bin (T1 early -> T5 late)", y = "Expression (z-score)",
       color = "Gene", title = "Hub gene dynamics along pseudotime (GSE151263 PBMC)") +
  theme_pub() +
  theme(axis.text.x = element_text(face = "plain"))
save_fig3(p_trend, "slingshot_gene_trend", w = 8, h = 5, dir = SL)
cat("    基因趋势图完成\n")

cat("\n===== Step3k (R 部分) 完成 =====\n")
