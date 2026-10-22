# ============================================================
# Step 3f —— GSE167363 CellChat 细胞通讯（HC vs Sepsis）
# 从注释好的 Seurat 对象导出 counts + meta，分组跑 CellChat
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ragg)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "08_singlecell/GSE167363/cellchat")
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

merged <- readRDS(file.path(BASE, "08_singlecell/GSE167363/GSE167363_seurat_annotated.rds"))
cat("细胞数:", ncol(merged), "\n"); print(table(merged$group, merged$celltype))

# 导出 raw counts + meta（供 CellChat / scTenifoldKnk / CellOracle 复用）
cnt <- GetAssayData(merged, layer = "counts")
meta <- merged[[]][, c("sample", "group", "timepoint", "celltype", "seurat_clusters")]
INP <- file.path(BASE, "08_singlecell/GSE167363/input")
dir.create(INP, showWarnings = FALSE)
writeMM(cnt, file.path(INP, "counts.mtx"))
writeLines(rownames(cnt), file.path(INP, "genes.tsv"))
write.csv(meta, file.path(INP, "meta.csv"))
cat("counts 导出:", dim(cnt), "\n")

# ---------- CellChat ----------
run_cellchat <- function(grp) {
  idx <- which(meta$group == grp)
  data <- cnt[, idx, drop = FALSE]
  m <- meta[idx, , drop = FALSE]
  m$labels <- droplevels(factor(m$celltype))
  cat(sprintf("=== %s: %d 细胞 ===\n", grp, length(idx)))
  cellchat <- createCellChat(object = data, meta = m, group.by = "labels")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat
}

cat("\n运行 HC CellChat...\n")
cc_HC <- run_cellchat("HC")
saveRDS(cc_HC, file.path(OUT, "CellChat_HC.rds"))
cat("\n运行 Sepsis CellChat...\n")
cc_Sep <- run_cellchat("Sepsis")
saveRDS(cc_Sep, file.path(OUT, "CellChat_Sepsis.rds"))

for (g in c("HC", "Sepsis")) {
  cc <- if (g == "HC") cc_HC else cc_Sep
  df.net <- subsetCommunication(cc)
  write.csv(df.net, file.path(OUT, paste0(g, "_interactions.csv")), row.names = FALSE)
  cat(sprintf("%s 显著配受体对: %d\n", g, nrow(df.net)))
}

# ---------- 可视化（单组）----------
groupSize <- table(cc_HC@idents)
net_cols <- c("#6A5C95", "#D99AAE", "#E7A889", "#A9CBB8", "#CFE3E6", "#8FA2D6")
par_save_net <- function(cc, grp, slots = c("count", "weight")) {
  for (s in slots) {
    png(file.path(OUT, paste0("net_", s, "_", grp, ".png")), width = 2000, height = 1800, res = 300)
    par(mar = c(4, 4, 2, 1))
    netVisual_circle(cc@net[[s]], vertex.weight = groupSize, weight.scale = TRUE,
                     label.edge = FALSE, title.name = paste(s, "(", grp, ")"),
                     vertex.label.cex = 0.7, color.use = net_cols[seq_along(groupSize)])
    dev.off()
  }
}
tryCatch(net_save <- netVisual_circle, error = function(e) NULL)
# 圈图
groupSize_HC <- table(cc_HC@idents)
groupSize_Sep <- table(cc_Sep@idents)
for (s in c("count", "weight")) {
  for (grp in c("HC", "Sepsis")) {
    cc <- if (grp == "HC") cc_HC else cc_Sep
    gs <- if (grp == "HC") groupSize_HC else groupSize_Sep
    png(file.path(OUT, paste0("net_", s, "_", grp, ".png")), width = 2000, height = 1800, res = 300)
    par(mar = c(2, 2, 2, 1))
    netVisual_circle(cc@net[[s]], vertex.weight = gs, weight.scale = TRUE,
                     label.edge = FALSE, title.name = paste(toupper(s), "-", grp),
                     vertex.label.cex = 0.8, color.use = net_cols[seq_along(gs)])
    dev.off()
  }
}

# ---------- 合并比较 ----------
cat("\n合并比较 HC vs Sepsis...\n")
object.list <- list(HC = cc_HC, Sepsis = cc_Sep)
cellchat_merged <- mergeCellChat(object.list, add.names = c("HC", "Sepsis"), cell.prefix = TRUE)
saveRDS(cellchat_merged, file.path(OUT, "CellChat_merged.rds"))

tryCatch({
  gg1 <- netVisual_heatmap(cellchat_merged)
  save_fig3(gg1, "diff_heatmap", w = 8, h = 6, dir = OUT)
}, error = function(e) cat("diff heatmap 失败:", conditionMessage(e), "\n"))

tryCatch({
  gg2 <- rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
  save_fig3(gg2, "rankNet_pathways", w = 7, h = 7, dir = OUT)
}, error = function(e) cat("rankNet 失败:", conditionMessage(e), "\n"))

tryCatch({
  pb <- netVisual_bubble(cellchat_merged, comparison = c(1, 2), angle.x = 45)
  save_fig3(pb, "diff_LR_bubble", w = 12, h = 10, dir = OUT)
}, error = function(e) cat("bubble 图失败:", conditionMessage(e), "\n"))

# 关键通路比较：MIF、CCL、CXCL、IL1 等炎症通路
tryCatch({
  for (pw in intersect(c("MIF", "CCL", "CXCL", "IL1", "IL6", "TNF", "SPP1"),
                       cellchat_merged@varnames$interaction)) {
    # 不在此处逐通路画, 后续汇总
  }
}, error = function(e) NULL)

cat("\n===== Step3f CellChat (GSE167363) 完成 =====\n")
cat("输出目录:", OUT, "\n")
