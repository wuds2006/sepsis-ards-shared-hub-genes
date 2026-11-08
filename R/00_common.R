# ============================================================
# 00_common.R —— 分析阶段共享：路径、配色、主题、工具函数
# 由早期 scripts/common.R 与 scripts/ajrccm_theme.R 合并而来，
# 合并原因：两份文件定义了同名的 palette / theme，散落引用易冲突。
#
# 路径约定：所有脚本假定工作目录 = 仓库根（Code/）。
#   Rscript R/10_deg.R        在 Code/ 下执行
# 输出目录沿用原分析阶段的编号（01_DEG/、05_ML_SHAP/ …），
# 已在 .gitignore 中排除，重跑后不会污染仓库。
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(RColorBrewer)
  library(ggsci)
})

# ---- 路径 ----------------------------------------------------
# 允许用环境变量指到别处（例如把大矩阵放在外置盘）
ROOT <- normalizePath(Sys.getenv("ARDS_ROOT", "."), mustWork = FALSE)
PROJ <- ROOT          # 旧脚本里的变量名，保留以兼容
setwd(ROOT)

# ---- 配色（低饱和，AJRCCM 印刷友好） -------------------------
PAL <- c(
  purple     = "#6A5C95",
  pink       = "#D99AAE",
  peach      = "#E7A889",
  green      = "#A9CBB8",
  lightblue  = "#CFE3E6",
  periwinkle = "#8FA2D6",
  darkblue   = "#5B608C"
)

PAL_UP    <- "#D99AAE"   # 上调
PAL_DOWN  <- "#8FA2D6"   # 下调
PAL_NS    <- "#C9CCD3"   # 不显著
PAL_GROUP <- c(Control = "#8FA2D6", Disease = "#D99AAE")

COL_DISEASE <- "#E64B35FF"
COL_CONTROL <- "#4DBBD5FF"
COL_UP      <- "#E64B35FF"
COL_DOWN    <- "#3C5488FF"
COL_NS      <- "grey75"

# NPG 分类色（面板超过 7 类时按需扩展）
pal_npg <- function(n) {
  cols <- ggsci::pal_npg("nrc")(10)
  if (n <= 10) cols[seq_len(n)] else colorRampPalette(cols)(n)
}

# ---- 版面（AJRCCM 双栏 178 mm，正文 7 pt） -------------------
AJRCCM_W_SINGLE  <- 89
AJRCCM_W_ONEHALF <- 140
AJRCCM_W_DOUBLE  <- 178
AJRCCM_FONT      <- "Helvetica"

HUB6 <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")

# ---- 主题 ----------------------------------------------------
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

theme_ajrccm <- function(base_size = 7, base_family = AJRCCM_FONT) {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid       = element_blank(),
      panel.border     = element_rect(color = "black", linewidth = 0.4),
      axis.line        = element_blank(),
      axis.ticks       = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(1.2, "mm"),
      axis.text        = element_text(color = "black", size = base_size),
      axis.title       = element_text(color = "black", size = base_size),
      plot.title       = element_text(size = base_size, face = "bold", hjust = 0),
      plot.subtitle    = element_text(size = base_size * 0.9),
      legend.text      = element_text(size = base_size * 0.9),
      legend.title     = element_text(size = base_size * 0.9),
      legend.key.size  = unit(3, "mm"),
      legend.margin    = margin(0, 0, 0, 0),
      plot.margin      = margin(2, 2, 2, 2),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(size = base_size, face = "bold",
                                      margin = margin(1, 1, 1, 1))
    )
}

# ---- 输出 ----------------------------------------------------
# PDF 用基础 pdf 设备 + Helvetica（PDF 核心字体）：文字保持可编辑，
# 不用 cairo_pdf —— 本机 Arial 的 CID 查找失败。
save_pdf <- function(plot, filename, width_mm, height_mm, ...) {
  grDevices::pdf(filename, width = width_mm / 25.4, height = height_mm / 25.4,
                 family = "Helvetica", useDingbats = FALSE, ...)
  print(plot)
  dev.off()
  invisible(filename)
}

save_fig <- function(p, name, w = 7, h = 6, dir = ".", dpi = 300) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".pdf")), p,
                          width = w, height = h, device = "pdf"))
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".png")), p,
                          width = w, height = h, dpi = dpi))
  invisible(p)
}

# ---- 数据工具 -------------------------------------------------
read_expr <- function(path) {
  expr <- read.delim(path, header = TRUE, check.names = FALSE,
                     stringsAsFactors = FALSE, row.names = 1, sep = "\t")
  expr <- as.matrix(expr)
  storage.mode(expr) <- "double"
  expr
}

# 探针 -> 基因：同名基因保留平均表达最高的探针
collapse_by_gene <- function(expr, gmap) {
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

# ---- 小工具 ---------------------------------------------------
sig_stars <- function(p) {
  ifelse(p < 1e-4, "****", ifelse(p < 1e-3, "***",
         ifelse(p < 1e-2, "**", ifelse(p < 0.05, "*", "ns"))))
}

gene_italic <- function(g) paste0("italic('", g, "')")

cat("00_common.R loaded; root =", ROOT, "\n")
