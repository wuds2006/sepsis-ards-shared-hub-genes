# ============================================================
# Step 1b —— WGCNA 加权基因共表达网络分析
# GSE76293(ARDS) 与 GSE69528(Sepsis) 分别构建，筛选疾病相关模块基因
# ============================================================
suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
})
allowWGCNAThreads(nThreads = 2)
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "02_WGCNA")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# 命令行参数: 数据集名
args <- commandArgs(trailingOnly = TRUE)
DATASET <- if (length(args) >= 1) args[1] else "GSE76293"
DISEASE <- if (DATASET == "GSE76293") "ARDS" else "Sepsis"

run_wgcna <- function(ds, disease, n_top = 3000, power = NULL) {
  cat("\n########## WGCNA:", ds, "##########\n")
  expr <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_gene_expr.rds")))
  grp  <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_grp.rds")))

  # 选择高变异基因
  mads <- apply(expr, 1, mad, na.rm = TRUE)
  expr <- expr[order(mads, decreasing = TRUE), , drop = FALSE]
  expr <- expr[1:min(n_top, nrow(expr)), , drop = FALSE]
  cat(sprintf("筛选 %d 个高变异基因, 样本 %d\n", nrow(expr), ncol(expr)))

  datExpr <- t(expr)
  # 样本聚类检测离群
  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) {
    datExpr <- datExpr[, gsg$goodGenes, drop = FALSE]
  }

  # 表型
  trait <- data.frame(disease = as.numeric(grp$group == disease))
  rownames(trait) <- grp$sample

  # 软阈值
  if (is.null(power)) {
    sft <- pickSoftThreshold(datExpr, powerVector = 1:20, verbose = 0,
                             networkType = "unsigned")
    power <- sft$powerEstimate
    if (is.na(power)) power <- 6
  }
  cat(sprintf("软阈值 power = %d\n", power))

  # 一步法网络构建 + 模块检测
  net <- blockwiseModules(datExpr, power = power, networkType = "unsigned",
                          TOMType = "unsigned", minModuleSize = 30,
                          reassignThreshold = 0, mergeCutHeight = 0.25,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = FALSE, verbose = 0, maxBlockSize = 20000)
  moduleColors <- labels2colors(net$colors)
  rm(net); gc()
  cat("模块数:", length(unique(moduleColors)), "\n")
  cat("模块大小:", paste(names(table(moduleColors)), table(moduleColors), collapse = " "), "\n")

  # 模块-表型相关性
  MElist <- moduleEigengenes(datExpr, moduleColors)$eigengenes
  MElist <- orderMEs(MElist)
  modTrait <- cor(MElist, trait$disease, use = "p")
  modTraitP <- corPvalueStudent(modTrait, nrow(datExpr))
  mt <- data.frame(module = sub("^ME", "", rownames(modTrait)),
                   cor = modTrait[, 1], p = modTraitP[, 1])
  mt <- mt[order(mt$p), ]
  cat("模块-疾病相关性(Top):\n"); print(head(mt, 10))

  # 显著相关模块 (p < 0.05)
  sig_modules <- mt$module[mt$p < 0.05]
  cat("显著相关模块:", paste(sig_modules, collapse = ", "), "\n")

  # 模块基因分配
  gene_module <- data.frame(gene = colnames(datExpr), module = moduleColors,
                            stringsAsFactors = FALSE)
  write.csv(gene_module, file.path(OUT, paste0(ds, "_gene_module.csv")), row.names = FALSE)

  # 疾病相关模块基因
  disease_genes <- gene_module$gene[gene_module$module %in% sig_modules]
  writeLines(disease_genes, file.path(OUT, paste0(ds, "_disease_module_genes.txt")))
  cat("疾病相关模块基因数:", length(disease_genes), "\n")

  # 保存 eigengene + trait 供后续
  saveRDS(list(gene_module = gene_module, sig_modules = sig_modules,
               disease_genes = disease_genes, mt = mt, power = power),
          file.path(OUT, paste0(ds, "_wgcna.rds")))

  # 模块-表型热图
  df <- data.frame(module = mt$module, cor = mt$cor, p = mt$p)
  df$module <- factor(df$module, levels = rev(df$module))
  p <- ggplot(df, aes(x = "ARDS/Sepsis", y = module, fill = cor)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f\n(p=%.2g)", cor, p)), size = 2.6) +
    scale_fill_gradient2(low = PAL_DOWN, mid = "white", high = PAL_UP, midpoint = 0) +
    labs(x = "", y = "Module", title = paste(ds, "模块-疾病相关性")) +
    theme_pub()
  save_fig(p, paste0(ds, "_module_trait"), w = 5, h = 7, dir = OUT)

  invisible(list(expr = expr, datExpr = datExpr, trait = trait, moduleColors = moduleColors,
                 mt = mt, disease_genes = disease_genes))
}

r1 <- run_wgcna(DATASET, DISEASE)

cat("\n===== Step1b WGCNA 完成:", DATASET, "=====\n")
cat("疾病相关模块基因数:", length(r1$disease_genes), "\n")
cat("OK\n")
