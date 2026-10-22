# ============================================================
# common.R —— 共享配色、主题、工具函数（AJRCCM 论文统一视觉规范）
# 配色: #6A5C95 #D99AAE #E7A889 #A9CBB8 #CFE3E6 #8FA2D6 #5B608C
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(RColorBrewer)
})

# 统一配色
PAL <- c(
  purple    = "#6A5C95",
  pink      = "#D99AAE",
  peach     = "#E7A889",
  green     = "#A9CBB8",
  lightblue = "#CFE3E6",
  periwinkle= "#8FA2D6",
  darkblue  = "#5B608C"
)

# 语义配色（差异表达等）
PAL_UP   <- "#D99AAE"   # 上调（红粉色系）
PAL_DOWN <- "#8FA2D6"   # 下调（蓝紫色系）
PAL_NS   <- "#C9CCD3"   # 不显著（灰）
PAL_GROUP <- c(Control = "#8FA2D6", Disease = "#D99AAE")

# 全局 ggplot 主题
theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 1),
      axis.text = element_text(color = "black"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

# 保存图：PDF(矢量, line art) + PNG(300dpi 预览)
save_fig <- function(p, name, w = 7, h = 6, dir = ".", dpi = 300) {
  # PDF 用默认 pdf 设备(不依赖 cairo/X11)
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf"))
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = dpi))
  invisible(p)
}

# 读取 GEO 系列矩阵清洗后的表达文件
read_expr <- function(path) {
  expr <- read.delim(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE,
                     row.names = 1, sep = "\t")
  expr <- as.matrix(expr)
  storage.mode(expr) <- "double"
  expr
}

# 探针 -> 基因(取平均表达最高的探针)
collapse_by_gene <- function(expr, gmap) {
  # gmap: named vector probe -> gene
  common <- intersect(rownames(expr), names(gmap))
  expr <- expr[common, , drop = FALSE]
  gene <- gmap[common]
  gmean <- rowMeans(expr, na.rm = TRUE)
  df <- data.frame(probe = common, gene = gene, mean = gmean, stringsAsFactors = FALSE)
  df <- df[order(df$mean, decreasing = TRUE), ]
  df <- df[!duplicated(df$gene), ]
  out <- expr[df$probe, , drop = FALSE]
  rownames(out) <- df$gene
  out
}

cat("common.R loaded: 配色与主题就绪\n")
