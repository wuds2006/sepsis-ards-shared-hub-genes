# ============================================================
# 重做 Graphical Abstract (三区域圆角矩形 + 中心枢纽)
# 参考 Desktop/Graphical abstract.pdf 样式
# ============================================================
library(ggplot2)
library(patchwork)
library(ggrepel)
library(dplyr)

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "AJRCCM_submission/04_graphical_abstract")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# 配色 (Nature 风格)
COL_GREEN <- "#A9C8A8"  # 左上 识别
COL_ORANGE <- "#E8B17B" # 右上 机制
COL_BLUE <- "#A3BCD3"   # 下 药物
COL_CENTER <- "#E8C8C8" # 中心

# ---- 画布 ----
# 用 ggplot 拼 3 个 panel + 中心 panel (patchwork)
# 简单方案：用 ggplot facet + 自定义 annotation，模仿参考图

# 方案：分别画三个 panel + 一个 center panel，patchwork 拼接 + annotation 加连接线

p_green <- ggplot() + theme_void() +
  # 框
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
           fill=COL_GREEN, alpha=0.5, color=COL_GREEN, linewidth=0) +
  annotate("text", x=5, y=9.5, label="Identification of potential core hub genes",
           fontface="bold", size=4, color="grey15") +
  annotate("text", x=2.5, y=8.4, label="GSE76293\nARDS\n(n=12 vs 12)",
           size=3, fontface="italic") +
  annotate("text", x=7.5, y=8.4, label="GSE69528\nSepsis\n(n=83 vs 28)",
           size=3, fontface="italic") +
  annotate("text", x=2.5, y=7.2, label="557 / 1,458\nDEGs", size=3) +
  annotate("text", x=7.5, y=7.2, label="13 / 11\nmodules", size=3) +
  annotate("text", x=5, y=6.0, label="179 shared candidate genes",
           fontface="bold", size=3.5, color="grey25") +
  annotate("text", x=5, y=5.0, label="11 ML classifiers + SHAP", size=3) +
  annotate("text", x=5, y=4.0, label="(AUC 0.75-1.00)",
           size=2.5, color="grey35") +
  annotate("text", x=5, y=2.5, label="AUC", fontface="bold",
           size=3, color="grey15") +
  # 简化柱状+折线
  annotate("segment", x=2, y=1.0, xend=8, yend=1.0, linetype="dashed", color="grey60") +
  annotate("segment", x=2, y=0.5, xend=8, yend=0.5, linetype="dashed", color="grey60") +
  annotate("segment", x=2.0, y=1.5, xend=2.5, yend=1.8, color="#6A5C95", linewidth=1.5) +
  annotate("segment", x=2.5, y=1.8, xend=3.0, yend=1.7, color="#6A5C95", linewidth=1.5) +
  annotate("segment", x=3.0, y=1.7, xend=3.5, yend=1.9, color="#6A5C95", linewidth=1.5) +
  annotate("segment", x=3.5, y=1.9, xend=4.0, yend=2.0, color="#6A5C95", linewidth=1.5) +
  annotate("segment", x=4.0, y=2.0, xend=4.5, yend=2.0, color="#6A5C95", linewidth=1.5) +
  annotate("segment", x=4.5, y=2.0, xend=5.0, yend=2.0, color="#6A5C95", linewidth=1.5) +
  annotate("segment", x=5.0, y=2.0, xend=5.5, yend=2.0, color="#6A5C95", linewidth=1.5) +
  annotate("segment", x=5.5, y=2.0, xend=6.0, yend=2.0, color="#6A5C95", linewidth=1.5) +
  annotate("segment", x=6.0, y=2.0, xend=6.5, yend=2.0, color="#6A5C95", linewidth=1.5) +
  annotate("segment", x=6.5, y=2.0, xend=7.0, yend=2.0, color="#6A5C95", linewidth=1.5) +
  annotate("segment", x=7.0, y=2.0, xend=7.5, yend=2.0, color="#6A5C95", linewidth=1.5) +
  annotate("segment", x=7.5, y=2.0, xend=8.0, yend=2.0, color="#6A5C95", linewidth=1.5) +
  xlim(0, 10) + ylim(0, 10)

p_orange <- ggplot() + theme_void() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
           fill=COL_ORANGE, alpha=0.5, color=COL_ORANGE, linewidth=0) +
  annotate("text", x=5, y=9.5, label="Exploration of molecular mechanisms",
           fontface="bold", size=4, color="grey15") +
  annotate("text", x=5, y=8.4, label="Single-cell transcriptomics",
           fontface="bold", size=3.2) +
  annotate("text", x=2.5, y=7.4, label="CellChat\n(MIF signaling)",
           size=2.8, fontface="italic") +
  annotate("text", x=7.5, y=7.4, label="Slingshot\n(monocyte-to-M2)",
           size=2.8, fontface="italic") +
  annotate("text", x=2.5, y=6.0, label="scTenifoldKnk", size=2.8, fontface="italic") +
  annotate("text", x=7.5, y=6.0, label="CellOracle\n(CEBPB KO)", size=2.8, fontface="italic") +
  annotate("text", x=5, y=4.6, label="Pathway enrichment", fontface="bold", size=3.2) +
  annotate("text", x=2.5, y=3.6, label="GO: complement\nactivation",
           size=2.8, fontface="italic") +
  annotate("text", x=7.5, y=3.6, label="M2 macrophage\npolarization",
           size=2.8, fontface="italic", color="#D99AAE") +
  annotate("text", x=5, y=2.2, label="shared immune microenvironment", fontface="bold",
           size=3.5, color="#6A5C95") +
  xlim(0, 10) + ylim(0, 10)

p_blue <- ggplot() + theme_void() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
           fill=COL_BLUE, alpha=0.5, color=COL_BLUE, linewidth=0) +
  annotate("text", x=5, y=8.5, label="Drug prediction & docking",
           fontface="bold", size=4, color="grey15") +
  annotate("text", x=2.5, y=7.3, label="CMap\ndrug repositioning",
           size=3, fontface="italic") +
  annotate("text", x=7.5, y=7.3, label="AutoDock Vina\nmolecular docking",
           size=3, fontface="italic") +
  annotate("text", x=5, y=5.8, label="Glucocorticoids", fontface="bold",
           size=3.5, color="#6A5C95") +
  annotate("text", x=5, y=4.6, label="(prednisolone, dexamethasone,\n hydrocortisone)",
           size=2.6, color="grey35") +
  annotate("text", x=5, y=3.0, label="VNN1 — flavonoid natural products",
           fontface="bold", size=3.5, color="#D99AAE") +
  annotate("text", x=5, y=1.6, label="VNN1–luteolin −10.01 kcal/mol\nGRB10–EGCG −9.61 kcal/mol",
           size=2.6, color="grey35") +
  xlim(0, 10) + ylim(0, 10)

p_center <- ggplot() + theme_void() +
  annotate("point", x=5, y=5, size=4, shape=21, fill=COL_CENTER, color="#6A5C95") +
  annotate("text", x=5, y=5.5, label="Six shared",
           fontface="bold", size=2.8, color="grey15") +
  annotate("text", x=5, y=4.7, label="hub genes",
           fontface="bold", size=2.8, color="grey15") +
  annotate("text", x=5, y=4.0, label="CDK5RAP2  CD177", size=2, color="grey25", fontface="italic") +
  annotate("text", x=5, y=3.5, label="GRB10  UPP1", size=2, color="grey25", fontface="italic") +
  annotate("text", x=5, y=3.0, label="EXOSC4  VNN1", size=2, color="grey25", fontface="italic") +
  xlim(0, 10) + ylim(0, 10)

# 拼接
# 上排: p_green, p_orange (左、右)
# 下排: p_blue (居中)
# 中心: p_center (覆盖在中间)

p_combo <- (p_green | p_orange) / p_blue +
  plot_annotation(
    title = "Shared Hub Genes and Immune Signatures of\nSepsis-Associated Acute Respiratory Distress Syndrome",
    theme = theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, color = "#2A2A3A"))
  )

# 保存
ggsave(file.path(OUT, "graphical_abstract.pdf"), p_combo, width = 16, height = 11)
ggsave(file.path(OUT, "graphical_abstract.png"), p_combo, width = 16, height = 11, dpi = 300)
cat("Graphical abstract 已保存\n")