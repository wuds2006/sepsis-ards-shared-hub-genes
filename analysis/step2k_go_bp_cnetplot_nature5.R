# step2k_go_bp_cnetplot_nature5.R
# 6 core hub genes GO BP cnetplot —— 按参考图样式：
#   kk 布局 + gene 按 logFC 渐变 + 边按 category 着色 + 3 个图例
# 配色：Nature 5 模块色 (#6A5C95 / #D99AAE / #E7A889 / #A9CBB8 / #8FA2D6)

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
  eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
)

# ---- foldChange（GSE76293 ARDS 发现队列 logFC）----
fc <- c(
  CDK5RAP2 = 3.187, CD177 = 3.883, GRB10 = 2.611,
  UPP1 = 1.257, EXOSC4 = 2.051, VNN1 = 2.245
)

# ---- Nature 5 模块色（5 功能模块 -> 颜色）----
module_color <- c(
  "pyrimidine metabolism" = "#6A5C95",
  "glycogen metabolism"   = "#D99AAE",
  "centrosome regulation" = "#E7A889",
  "RNA processing"        = "#A9CBB8",
  "neutrophil immunity"   = "#8FA2D6"
)

showCat <- 15
top_desc <- ego@result$Description[1:showCat]
module_map <- c(
  "nucleoside phosphate catabolic process"                       = "pyrimidine metabolism",
  "organophosphate catabolic process"                            = "pyrimidine metabolism",
  "negative regulation of glycogen biosynthetic process"         = "glycogen metabolism",
  "CMP metabolic process"                                        = "pyrimidine metabolism",
  "UMP biosynthetic process"                                     = "pyrimidine metabolism",
  "pyrimidine ribonucleotide catabolic process"                  = "pyrimidine metabolism",
  "negative regulation of centrosome duplication"                = "centrosome regulation",
  "dTMP metabolic process"                                       = "pyrimidine metabolism",
  "negative regulation of centrosome cycle"                      = "centrosome regulation",
  "nuclear mRNA surveillance"                                    = "RNA processing",
  "pyrimidine ribonucleoside monophosphate biosynthetic process" = "pyrimidine metabolism",
  "sno(s)RNA 3'-end processing"                                  = "RNA processing",
  "sno(s)RNA processing"                                         = "RNA processing",
  "negative regulation of glycogen metabolic process"            = "glycogen metabolism",
  "neutrophil extravasation"                                     = "neutrophil immunity"
)

# category 节点颜色（长度 15，按 top_desc 顺序对应 5 模块色）
cat_cols <- unname(module_color[module_map[top_desc]])

# 边按 category 着色：让边也用相同的 5 模块色（与 category 节点一致）
edge_color_map <- setNames(cat_cols, top_desc)

# ---- cnetplot：kk 布局 + foldChange + 边按 category + 5 模块配色 ----
p <- cnetplot(
  ego,
  showCategory   = showCat,
  foldChange     = fc,
  layout         = "kk",
  color_edge     = "category",
  color_category = cat_cols,
  node_label     = "all"
) +
  # gene 节点：logFC 3 色渐变（Nature 5 色：蓝 -> 紫 -> 桃）
  scale_color_gradientn(
    colors = c("#8FA2D6", "#6A5C95", "#D99AAE"),
    name   = expression(log[2]~"(Fold change)")
  ) +
  # 边的离散 scale：覆盖 colour_ggnewscale_1，让边用 5 模块色（与节点一致）
  scale_color_manual(
    values      = edge_color_map,
    aesthetics = "colour_ggnewscale_1",
    name        = "category"
  ) +
  theme(legend.position = "right")

# ---- 导出 ----
w_in <- 9.5
h_in <- 7.5

pdf(file.path(OUT, "GO_BP_cnetplot_nature5_6hub.pdf"), width = w_in, height = h_in, family = "Helvetica")
print(p)
dev.off()

ragg::agg_tiff(file.path(OUT, "GO_BP_cnetplot_nature5_6hub.tiff"),
               width = w_in, height = h_in, units = "in", res = 1000, compression = "lzw")
print(p)
dev.off()

ragg::agg_png(file.path(OUT, "GO_BP_cnetplot_nature5_6hub.png"),
              width = w_in, height = h_in, units = "in", res = 300)
print(p)
dev.off()

cat("\n完成。输出目录:", OUT, "\n")
cat("PDF / TIFF(1000dpi) / PNG(300dpi) 已生成\n")
