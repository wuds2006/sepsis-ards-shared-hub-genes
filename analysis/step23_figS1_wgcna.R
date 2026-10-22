# ============================================================
# Step 23 —— Figure S1: WGCNA 软阈值选择 + 模块聚类树图
# 两队列 (GSE76293 ARDS / GSE69528 Sepsis)，面板输出到
# 02_WGCNA/figS1_parts/ 供 Python 拼装 composite
# 复现 step1b_WGCNA.R 参数: top3000 MAD, unsigned, minModuleSize=30,
# mergeCutHeight=0.25, power 取原运行值 (02_WGCNA/*_wgcna.rds)
# ============================================================
suppressPackageStartupMessages({
  library(WGCNA); library(ggplot2)
})
allowWGCNAThreads(nThreads = 2)
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
PARTS <- file.path(BASE, "02_WGCNA", "figS1_parts")
dir.create(PARTS, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

datasets <- list(list(ds = "GSE76293", disease = "ARDS"),
                list(ds = "GSE69528", disease = "Sepsis"))

for (k in seq_along(datasets)) {
  ds <- datasets[[k]]$ds; disease <- datasets[[k]]$disease
  cat("\n##########", ds, "##########\n")

  expr <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_gene_expr.rds")))
  grp  <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_grp.rds")))

  mads <- apply(expr, 1, mad, na.rm = TRUE)
  expr <- expr[order(mads, decreasing = TRUE), , drop = FALSE]
  expr <- expr[1:min(3000, nrow(expr)), , drop = FALSE]
  datExpr <- t(expr)
  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) datExpr <- datExpr[, gsg$goodGenes, drop = FALSE]

  # 原运行使用的软阈值
  power_used <- readRDS(file.path(BASE, "02_WGCNA", paste0(ds, "_wgcna.rds")))$power
  cat("原 soft power =", power_used, "\n")

  # ---- 软阈值扫描 ----
  sft <- pickSoftThreshold(datExpr, powerVector = 1:20, verbose = 0,
                           networkType = "unsigned")
  fi <- sft$fitIndices
  fi$SFT <- -sign(fi$slope) * fi$SFT.R.sq

  p_r2 <- ggplot(fi, aes(x = Power, y = SFT)) +
    geom_hline(yintercept = 0.8, color = "#b00020", linetype = "dashed", linewidth = 0.45) +
    geom_line(color = "#8FA2D6", linewidth = 0.5) +
    geom_point(size = 1.7, color = "#5B608C") +
    geom_point(data = subset(fi, Power == power_used), size = 3.2,
               shape = 21, fill = "#6A5C95", color = "#2b2b2b", stroke = 0.5) +
    annotate("text", x = power_used + 0.5, y = fi$SFT[fi$Power == power_used] - 0.06,
             label = paste0("power = ", power_used),
             size = 2.6, color = "#6A5C95", fontface = "bold", hjust = 0) +
    scale_x_continuous(breaks = seq(0, 20, 4)) +
    labs(x = "Soft-thresholding power", y = "Scale-free topology fit (R²)",
         title = paste0(ds, " (", disease, "): scale-free topology fit")) +
    theme_pub() +
    theme(plot.title = element_text(size = 9, face = "bold"),
          axis.title = element_text(size = 7.5), axis.text = element_text(size = 6.5))
  save_fig(p_r2, paste0(letters[3*k-2], "_sft_", ds), w = 4.4, h = 3.2, dir = PARTS)

  p_mc <- ggplot(fi, aes(x = Power, y = mean.k.)) +
    geom_line(color = "#D99AAE", linewidth = 0.5) +
    geom_point(size = 1.7, color = "#A66E8C") +
    geom_point(data = subset(fi, Power == power_used), size = 3.2,
               shape = 21, fill = "#D99AAE", color = "#2b2b2b", stroke = 0.5) +
    scale_x_continuous(breaks = seq(0, 20, 4)) +
    labs(x = "Soft-thresholding power", y = "Mean connectivity",
         title = paste0(ds, ": mean connectivity")) +
    theme_pub() +
    theme(plot.title = element_text(size = 9, face = "bold"),
          axis.title = element_text(size = 7.5), axis.text = element_text(size = 6.5))
  save_fig(p_mc, paste0(letters[3*k-1], "_conn_", ds), w = 4.4, h = 3.2, dir = PARTS)

  # ---- 模块聚类树 ----
  net <- blockwiseModules(datExpr, power = power_used, networkType = "unsigned",
                          TOMType = "unsigned", minModuleSize = 30,
                          reassignThreshold = 0, mergeCutHeight = 0.25,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = FALSE, verbose = 0, maxBlockSize = 20000)
  moduleColors <- labels2colors(net$colors)
  cat("模块数:", length(unique(moduleColors)), "\n")

  pdf(file.path(PARTS, paste0(letters[3*k], "_dendro_", ds, ".pdf")),
      width = 7.2, height = 3.4)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module", dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = paste0(ds, " (", disease,
                                    "): gene clustering dendrogram and module colors"))
  dev.off()
  rm(net); gc()
}

cat("\n===== Step 23 Figure S1 面板完成 =====\n")
cat("OK\n")
