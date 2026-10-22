# ============================================================
# Step28: Figure 7 (GSVA top20) 矢量重绘
#   基于已保存的 GSVA scores/top20/分组数据，用 pheatmap 重新生成
#   矢量 PDF（Helvetica 文字，可编辑），复刻 step2c3 的绘图样式：
#     - 行 = top20 免疫 GO-BP 通路（行聚类，scale="row"）
#     - 列 = 样本（Control 在前、疾病在后），列不聚类
#     - 配色 #2166AC-white-#B2182B；分组注释条 Control #8FA2D6 / 疾病 #D99AAE
#   输出到 07_hub_immune/gsva_top20/vector/（原始全尺寸，后续再等比
#   缩放到宽 120mm 并导出 PDF/PNG/TIFF）。
# ============================================================
suppressPackageStartupMessages({
  library(pheatmap)
  library(grid)
})

ROOT <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT <- file.path(ROOT, "07_hub_immune/gsva_top20/vector")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

build_panel <- function(ds, disease) {
  scores <- read.csv(file.path(ROOT, sprintf("07_hub_immune/gsva_top20/%s_GSVA_scores.csv", ds)),
                     row.names = 1, check.names = FALSE)
  top <- read.csv(file.path(ROOT, sprintf("07_hub_immune/gsva_top20/%s_GSVA_top20.csv", ds)))
  grp_raw <- read.csv(file.path(ROOT, sprintf("01_DEG/%s_groups.csv", ds)),
                      stringsAsFactors = FALSE)
  # 仅保留 疾病 vs 对照（与 step2c3 口径一致）
  grp_raw <- grp_raw[grp_raw$group %in% c("Control", disease), , drop = FALSE]
  grp_raw$group <- factor(grp_raw$group, levels = c("Control", disease))
  common <- intersect(grp_raw$sample, colnames(scores))
  grp_raw <- grp_raw[match(common, grp_raw$sample), , drop = FALSE]
  scores <- scores[, common, drop = FALSE]
  ord <- order(grp_raw$group)
  mat <- scores[top$pathway, ord, drop = FALSE]
  rownames(mat) <- sub("^GOBP_", "", rownames(mat))
  grp_ord <- droplevels(grp_raw$group[ord])
  ann <- data.frame(Group = grp_ord, row.names = colnames(mat))
  lv <- levels(grp_ord)
  ann_colors <- list(Group = c(Control = "#8FA2D6", ARDS = "#D99AAE", Sepsis = "#D99AAE")[lv])
  palette <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  p <- pheatmap(mat, annotation_col = ann, annotation_colors = ann_colors,
                color = palette, scale = "row",
                cluster_cols = FALSE, show_colnames = FALSE,
                fontsize = 8, fontsize_row = 8, fontsize_col = 8,
                main = sprintf("%s: Immune GO-BP GSVA (top 20)", disease),
                silent = TRUE)
  cat(sprintf("[OK] %s (%s): %d 通路 x %d 样本\n", ds, disease, nrow(mat), ncol(mat)))
  p
}

pA <- build_panel("GSE76293", "ARDS")
pB <- build_panel("GSE69528", "Sepsis")

# ---- 独立面板 A / B（全尺寸矢量，11 x 9 in，同原始 PNG 设备尺寸）----
save_one <- function(p, fn) {
  fp <- file.path(OUT, fn)
  pdf(fp, width = 11, height = 9, family = "Helvetica", useDingbats = FALSE)
  grid.newpage(); grid.draw(p$gtable)
  dev.off()
  cat("saved:", fp, "\n")
}
save_one(pA, "A_GSVA_top20_ARDS_full.pdf")
save_one(pB, "B_GSVA_top20_Sepsis_full.pdf")

# ---- 组合图（A/B 并排 + 字母标签独立于面板标题，17 x 8 in 矢量）----
fp <- file.path(OUT, "Figure5_composite_full.pdf")
pdf(fp, width = 17, height = 8, family = "Helvetica", useDingbats = FALSE)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
# 面板上沿下移 8 mm，为字母标签留出空间
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1,
                      y = unit(1, "npc") - unit(8, "mm"), just = c("centre", "top"),
                      height = unit(1, "npc") - unit(10, "mm")))
grid.draw(pA$gtable)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2,
                      y = unit(1, "npc") - unit(8, "mm"), just = c("centre", "top"),
                      height = unit(1, "npc") - unit(10, "mm")))
grid.draw(pB$gtable)
upViewport()
upViewport()
grid.text("A", x = unit(1.5, "mm"), y = unit(1, "npc") - unit(0.5, "mm"),
          just = c("left", "top"),
          gp = gpar(fontsize = 22, fontface = "bold", fontfamily = "Helvetica"))
grid.text("B", x = unit(0.5, "npc") + unit(1.5, "mm"),
          y = unit(1, "npc") - unit(0.5, "mm"),
          just = c("left", "top"),
          gp = gpar(fontsize = 22, fontface = "bold", fontfamily = "Helvetica"))
dev.off()
cat("saved:", fp, "\n")
