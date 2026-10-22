#!/usr/bin/env Rscript
# Step27 — 审稿修复分析（R 部分）
# 修复 II-1: 发现队列内部跨疾病候选基因验证（不使用外部验证队列的筛选规则）
# 修复 II-6: 六 hub 基因的 WGCNA kME（模块隶属度）与共表达度中心性
# 输出:
#   05_ML_SHAP/discovery_crossdisease_AUC.csv
#   02_WGCNA/hub_gene_kME.csv

suppressMessages({
  library(WGCNA)
})
options(stringsAsFactors = FALSE)
ROOT <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"

## ---------- 工具：AUC (Mann-Whitney) ----------
auc_mw <- function(x, g) {
  # g: 逻辑向量，TRUE = 疾病
  x1 <- x[g]; x0 <- x[!g]
  if (length(x1) == 0 || length(x0) == 0) return(c(NA, NA))
  wi <- wilcox.test(x1, x0, exact = TRUE)
  n1 <- length(x1); n0 <- length(x0)
  auc <- wi$statistic / (n1 * n0)
  # R exact wilcoxon CDF 在极小 P 时存在浮点负值；截断至机器精度下限
  p <- pmax(wi$p.value, .Machine$double.eps)
  c(auc = as.numeric(auc), p = p)
}

## ---------- Part A: 发现队列内部跨疾病候选验证 ----------
ARD_CAND <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1", "DYRK2", "TCTN1", "MUC1")
SEP_CAND <- c("METTL9", "MTHFD2", "TCTN1", "BATF", "PLP2", "ZSCAN18")
HUB6 <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")

expr76 <- readRDS(file.path(ROOT, "01_DEG/GSE76293_gene_expr.rds"))
expr95 <- readRDS(file.path(ROOT, "01_DEG/GSE69528_gene_expr.rds"))
g76 <- read.csv(file.path(ROOT, "01_DEG/GSE76293_groups.csv"))
g95 <- read.csv(file.path(ROOT, "01_DEG/GSE69528_groups.csv"))

# 对齐样本（groups CSV 含 EXCLUDE 行，expr 矩阵仅含策划后样本）
g76 <- g76[g76$sample %in% colnames(expr76), ]
g95 <- g95[g95$sample %in% colnames(expr95), ]
expr76 <- expr76[, g76$sample]; expr95 <- expr95[, g95$sample]
d76 <- g76$group == "ARDS"; d95 <- g95$group == "Sepsis" | g95$group == "Disease"
if (!any(d95)) d95 <- !grepl("Control|Healthy", g95$group, ignore.case = TRUE)
cat("GSE76293 groups:", table(g76$group), "\n")
cat("GSE69528 groups:", table(g95$group), "\n")

# ARDS 候选基因 -> 在 sepsis 发现队列中检验（跨疾病、发现队列内部）
res <- data.frame()
for (gn in unique(c(ARD_CAND, SEP_CAND, HUB6))) {
  if (gn %in% rownames(expr95)) {
    a <- auc_mw(as.numeric(expr95[gn, ]), d95)
  } else a <- c(NA, NA)
  if (gn %in% rownames(expr76)) {
    b <- auc_mw(as.numeric(expr76[gn, ]), d76)
  } else b <- c(NA, NA)
  res <- rbind(res, data.frame(
    gene = gn,
    list = ifelse(gn %in% ARD_CAND, "ARDS_candidate", ""),
    AUC_GSE76293_ARDSdiscovery = round(b[1], 3),
    P_GSE76293 = signif(b[2], 3),
    AUC_GSE69528_SepsisDiscovery = round(a[1], 3),
    P_GSE69528 = signif(a[2], 3)
  ))
}
res$list[res$gene %in% SEP_CAND] <- paste(res$list[res$gene %in% SEP_CAND], "Sepsis_candidate")
write.csv(res, file.path(ROOT, "05_ML_SHAP/discovery_crossdisease_AUC.csv"), row.names = FALSE)
print(res)

## ---------- Part B: hub 基因 kME（两个发现队列） ----------
kme_list <- list()
for (ds in c("GSE76293", "GSE69528")) {
  w <- readRDS(file.path(ROOT, sprintf("02_WGCNA/%s_wgcna.rds", ds)))
  ex <- if (ds == "GSE76293") expr76 else expr95
  gm <- w$gene_module
  genes <- intersect(gm$gene, rownames(ex))
  ex <- ex[genes, , drop = FALSE]
  # WGCNA 输入要行=样本
  dat <- t(ex)
  me <- moduleEigengenes(dat, colors = gm$module[match(genes, gm$gene)])
  MEs <- orderMEs(me$eigengenes)
  kme_all <- signedKME(dat, MEs, outputColumnName = "kME.")
  # 基因所在模块 + 该模块的 kME
  hub_in <- intersect(HUB6, genes)
  for (gn in hub_in) {
    mod <- gm$module[gm$gene == gn]
    me_col <- paste0("kME.", mod)
    if (me_col %in% colnames(kme_all)) {
      k <- kme_all[gn, me_col]
    } else k <- NA
    # 度中心性：与 179 共享候选基因中显著相关(|r|>0.5 / >0.7)的基因数
    cg <- readLines(file.path(ROOT, "03_intersection/Shared_CGs.txt"))
    cg <- intersect(cg, genes); cg <- setdiff(cg, gn)
    if (length(cg) > 1) {
      cm <- cor(t(ex[c(gn, cg), , drop = FALSE]))
      deg5 <- sum(abs(cm[gn, -1]) > 0.5)
      deg7 <- sum(abs(cm[gn, -1]) > 0.7)
    } else { deg5 <- NA; deg7 <- NA }
    kme_list[[length(kme_list) + 1]] <- data.frame(
      cohort = ds, gene = gn, module = mod,
      kME = round(as.numeric(k), 3),
      module_cor_disease = w$mt$cor[w$mt$module == mod],
      module_p = w$mt$p[w$mt$module == mod],
      degree_r0.5_in_179CG = deg5,
      degree_r0.7_in_179CG = deg7,
      n_179CG = length(cg)
    )
  }
  cat(ds, "module size of hub genes done\n")
}
kme <- do.call(rbind, kme_list)
write.csv(kme, file.path(ROOT, "02_WGCNA/hub_gene_kME.csv"), row.names = FALSE)
print(kme)
cat("\n[done] step27 R part\n")
