# step2j_go_bp_cnetplot_fc.R
# 6 core hub genes GO BP enrichment + cnetplot（foldChange 着色 + colorEdge + kk 布局）
# 基因颜色：低 #CFE3E6 -> 高 #D99AAE；圆形节点；边按 category 带颜色；layout = kk

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT <- file.path(BASE, "AJRCCM_submission/07_figure_exports/GO_BP_cnetplot")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

hub <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")
eg <- bitr(hub, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

ego <- enrichGO(
  eg$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)
res <- as.data.frame(ego)
cat("GO BP 显著 term 数:", nrow(res), "\n")

# ---- 6 基因 foldChange（GSE76293 ARDS 发现队列 logFC）----
fc <- c(
  CDK5RAP2 = 3.187,
  CD177    = 3.883,
  GRB10    = 2.611,
  UPP1     = 1.257,
  EXOSC4   = 2.051,
  VNN1     = 2.245
)

# ---- cnetplot：kk 布局 + foldChange 着色 + 边带颜色 ----
showCat <- 15  # 覆盖全部 6 基因（CD177 的 top term 排第 15）

p <- cnetplot(
  ego,
  showCategory = showCat,
  foldChange = fc,
  layout = "kk",
  color_edge = "category",
  node_label = "all"
) +
  # 基因节点渐变：低表达 #CFE3E6 -> 高表达 #D99AAE
  scale_color_gradient(low = "#CFE3E6", high = "#D99AAE", name = "log2 fold change") +
  theme(legend.position = "right")

# ---- 导出 ----
w_in <- 9.0
h_in <- 7.5

# 1) PDF（矢量）
pdf(file.path(OUT, "GO_BP_cnetplot_fc_6hub.pdf"), width = w_in, height = h_in, family = "Helvetica")
print(p)
dev.off()

# 2) TIFF（line art >=1000 dpi，LZW）
ragg::agg_tiff(file.path(OUT, "GO_BP_cnetplot_fc_6hub.tiff"),
               width = w_in, height = h_in, units = "in", res = 1000, compression = "lzw")
print(p)
dev.off()

# 3) PNG（彩色 >=300 dpi）
ragg::agg_png(file.path(OUT, "GO_BP_cnetplot_fc_6hub.png"),
              width = w_in, height = h_in, units = "in", res = 300)
print(p)
dev.off()

cat("\n完成。输出目录:", OUT, "\n")
cat("PDF / TIFF(1000dpi) / PNG(300dpi) 已生成\n")
