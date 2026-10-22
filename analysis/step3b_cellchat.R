# Step3b: CellChat 细胞通讯分析 (ARDS vs Sepsis)
suppressPackageStartupMessages({
  library(CellChat)
  library(Matrix)
  library(dplyr)
})

setwd("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis")
OUT <- "08_singlecell/cellchat"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# 读数据
cat("读取 count 矩阵...\n")
counts <- readMM("08_singlecell/cellchat_input/counts.mtx")  # 基因 x 细胞
genes <- readLines("08_singlecell/cellchat_input/genes.tsv")
meta <- read.csv("08_singlecell/cellchat_input/meta.csv", row.names = 1)
rownames(counts) <- genes
colnames(counts) <- rownames(meta)
cat("counts:", dim(counts), "| meta:", dim(meta), "\n")

run_cellchat <- function(grp) {
  idx <- which(meta$group == grp)
  data <- counts[, idx, drop = FALSE]
  m <- meta[idx, , drop = FALSE]
  m$labels <- m$celltype
  cat(sprintf("=== %s: %d 细胞 ===\n", grp, length(idx)))

  cellchat <- createCellChat(object = data, meta = m, group.by = "labels")
  cellchat@DB <- CellChatDB.human
  # 用所有基因(已QC过滤)
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # 通讯概率(非空间数据, raw.use=TRUE 用 data.signaling)
  cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat
}

cat("\n运行 ARDS CellChat...\n")
cc_ARDS <- run_cellchat("ARDS")
saveRDS(cc_ARDS, file.path(OUT, "CellChat_ARDS.rds"))

cat("\n运行 Sepsis CellChat...\n")
cc_Sepsis <- run_cellchat("Sepsis")
saveRDS(cc_Sepsis, file.path(OUT, "CellChat_Sepsis.rds"))

# 交互强度
for (g in c("ARDS", "Sepsis")) {
  cc <- if (g == "ARDS") cc_ARDS else cc_Sepsis
  df.net <- subsetCommunication(cc)
  write.csv(df.net, file.path(OUT, paste0(g, "_interactions.csv")), row.names = FALSE)
  cat(sprintf("%s 显著配受体对: %d\n", g, nrow(df.net)))
}

# 合并比较
cat("\n合并比较 ARDS vs Sepsis...\n")
object.list <- list(ARDS = cc_ARDS, Sepsis = cc_Sepsis)
cellchat_merged <- mergeCellChat(object.list, add.names = c("ARDS", "Sepsis"), cell.prefix = TRUE)

# 差异交互强度
net_diff <- netComparison_diff <- NULL
tryCatch({
  df.net_diff <- netAnalysis_computeCentrality <- NULL
  # 总交互强度差异
  gg1 <- netVisual_heatmap(cellchat_merged)
  ggsave(file.path(OUT, "diff_heatmap.png"), gg1, width = 8, height = 6, dpi = 300)
  # 差异信号通路排名
  gg2 <- rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
  ggsave(file.path(OUT, "rankNet_pathways.png"), gg2, width = 7, height = 6, dpi = 300)
  # 差异配受体对
  df_diff <- netAnalysis_compareSignaling <- NULL
  tryCatch({
    p <- netVisual_bubble(cellchat_merged, comparison = c(1, 2), angle.x = 45)
    ggsave(file.path(OUT, "diff_LR_bubble.png"), p, width = 12, height = 10, dpi = 300)
  }, error = function(e) cat("bubble 图失败:", conditionMessage(e), "\n"))
}, error = function(e) cat("合并比较失败:", conditionMessage(e), "\n"))

saveRDS(cellchat_merged, file.path(OUT, "CellChat_merged.rds"))
cat("\n===== Step3b CellChat 完成 =====\n")
cat("输出目录:", OUT, "\n")
