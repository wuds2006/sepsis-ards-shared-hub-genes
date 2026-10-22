# step2i_go_bp_cnetplot.R
# 6 core hub genes GO Biological Process enrichment + circular cnetplot
# Nature color palette; PDF (vector) / TIFF (1000 dpi) / PNG (300 dpi)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(patchwork)
  library(igraph)
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
write.csv(res, file.path(OUT, "GO_BP_enrichment_6hub.csv"), row.names = FALSE)
cat("GO BP 显著 term 数:", nrow(res), "\n")

# ---- showCategory = 15，覆盖全部 6 基因（CD177 的 top term 排第 15）----
showCat <- 15
top_desc <- res$Description[1:showCat]

# 每个 term 归属的功能模块（5 大模块，用于统一配色）
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

module_color <- c(
  "pyrimidine metabolism" = "#6A5C95",
  "glycogen metabolism"   = "#D99AAE",
  "centrosome regulation" = "#E7A889",
  "RNA processing"        = "#A9CBB8",
  "neutrophil immunity"   = "#8FA2D6"
)

# category 节点颜色（长度 = showCat，按 top_desc 顺序对应功能模块）
cat_cols <- unname(module_color[module_map[top_desc]])

# gene 节点统一深色
gene_col <- "#5B608C"

# ---- 主图：cnetplot 环状布局 ----
p <- cnetplot(
  ego,
  showCategory = showCat,
  color_category = cat_cols,
  color_item = gene_col,
  node_label = "all",
  layout = igraph::layout_in_circle,
  color_edge = "grey70"
) +
  theme(legend.position = "none")

# ---- category 功能模块图例（独立 ggplot，patchwork 拼在右侧）----
legend_df <- data.frame(
  module = names(module_color),
  y = seq_along(module_color),
  stringsAsFactors = FALSE
)
legend_p <- ggplot(legend_df, aes(x = 0, y = y)) +
  geom_point(aes(color = module), size = 4) +
  geom_text(aes(x = 0.35, label = module), hjust = 0, size = 3.2) +
  scale_color_manual(values = unname(module_color)) +
  scale_y_continuous(limits = c(0.5, 5.5), breaks = NULL) +
  scale_x_continuous(limits = c(-0.3, 3.8), breaks = NULL) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(0, 0, 0, 6))

# ---- 拼接：主图 + 右侧图例 ----
p_final <- p + legend_p + plot_layout(widths = c(5, 1.6))

# ---- 导出 ----
w_in <- 9.5
h_in <- 8.5

# 1) PDF（矢量）
pdf(file.path(OUT, "GO_BP_cnetplot_6hub.pdf"), width = w_in, height = h_in, family = "Helvetica")
print(p_final)
dev.off()

# 2) TIFF（line art >=1000 dpi，LZW）
ragg::agg_tiff(file.path(OUT, "GO_BP_cnetplot_6hub.tiff"),
               width = w_in, height = h_in, units = "in", res = 1000, compression = "lzw")
print(p_final)
dev.off()

# 3) PNG（彩色 >=300 dpi）
ragg::agg_png(file.path(OUT, "GO_BP_cnetplot_6hub.png"),
              width = w_in, height = h_in, units = "in", res = 300)
print(p_final)
dev.off()

cat("\n完成。输出目录:", OUT, "\n")
cat("PDF / TIFF(1000dpi) / PNG(300dpi) 已生成\n")
