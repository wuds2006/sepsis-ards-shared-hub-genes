# Step3b 补充: CellChat 差异可视化 (ARDS vs Sepsis)
suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(ComplexHeatmap)
})

setwd("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis")
OUT <- "08_singlecell/cellchat"

cc_ARDS <- readRDS(file.path(OUT, "CellChat_ARDS.rds"))
cc_Sepsis <- readRDS(file.path(OUT, "CellChat_Sepsis.rds"))

object.list <- list(ARDS = cc_ARDS, Sepsis = cc_Sepsis)
cellchat_merged <- mergeCellChat(object.list, add.names = c("ARDS", "Sepsis"), cell.prefix = TRUE)

# 1) 差异配受体对 bubble (ggplot)
tryCatch({
  p <- netVisual_bubble(cellchat_merged, comparison = c(1, 2), angle.x = 45,
                        remove.isolate = TRUE)
  ggsave(file.path(OUT, "diff_LR_bubble.png"), p, width = 13, height = 11, dpi = 300)
  cat("bubble 图已保存\n")
}, error = function(e) cat("bubble 失败:", conditionMessage(e), "\n"))

# 2) 差异信号通路排名 (ggplot)
tryCatch({
  gg2 <- rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
  ggsave(file.path(OUT, "rankNet_pathways.png"), gg2, width = 8, height = 7, dpi = 300)
  cat("rankNet 图已保存\n")
}, error = function(e) cat("rankNet 失败:", conditionMessage(e), "\n"))

# 3) 总交互强度差异 heatmap (ComplexHeatmap)
tryCatch({
  gg1 <- netVisual_heatmap(cellchat_merged)
  pdf(file.path(OUT, "diff_heatmap.pdf"), width = 8, height = 6)
  draw(gg1)
  dev.off()
  png(file.path(OUT, "diff_heatmap.png"), width = 2000, height = 1500, res = 300)
  draw(gg1)
  dev.off()
  cat("heatmap 已保存\n")
}, error = function(e) cat("heatmap 失败:", conditionMessage(e), "\n"))

# 4) 导出差异配受体对表 (分别导出的 interactions 已含, 这里跳过 merged 合并)
cat("差异交互表: 见 ARDS_interactions.csv 与 Sepsis_interactions.csv\n")

cat("\n===== CellChat 差异可视化完成 =====\n")
