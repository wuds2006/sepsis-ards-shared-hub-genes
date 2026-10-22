# ============================================================
# Step 49: CD177 消融分析 (M3)
# ARDS 模型去掉 CD177 后 5 基因 coefficient-transfer 在 GSE243066 (PBMC) 的 AUC
# 流程与 step2e_nomogram.R 完全一致: lrm(penalty=1) + 训练集 z-score 参数
# ============================================================
suppressPackageStartupMessages({
  library(rms)
  library(pROC)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "07_hub_immune")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

HUB6 <- c("CDK5RAP2","CD177","GRB10","UPP1","EXOSC4","VNN1")
HUB5 <- setdiff(HUB6, "CD177")

# ---- 数据加载（与 step47 loader 一致） ----
d76293 <- list(expr = readRDS(file.path(BASE, "01_DEG", "GSE76293_gene_expr.rds")),
               grp  = readRDS(file.path(BASE, "01_DEG", "GSE76293_grp.rds"))$group)
# GSE243066 ARDS 验证（与 step2e_nomogram.R 完全同款的 loader）
dd <- file.path(BASE, "00_rawdata/validation/GSE243066")
files <- list.files(dd, pattern = "raw_count.txt.gz$", full.names = TRUE)
mat <- NULL
for (f in files) {
  x <- read.delim(gzfile(f), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(x)[2] <- basename(f)
  mat <- if (is.null(mat)) x else merge(mat, x, by = "Gene_symbol", all = TRUE)
}
rownames(mat) <- mat$Gene_symbol; mat <- as.matrix(mat[, -1, drop = FALSE])
mat[is.na(mat)] <- 0; storage.mode(mat) <- "double"
sn <- colnames(mat)
g243 <- ifelse(grepl("_A[0-9]", sn), "ARDS", ifelse(grepl("_H[0-9]", sn), "Control", "Exclude"))
cpm <- t(t(mat) / colSums(mat) * 1e6); lc <- log2(cpm + 1)
keep <- g243 %in% c("ARDS", "Control"); lc <- lc[, keep, drop = FALSE]; g243 <- g243[keep]
d243066 <- list(expr = lc, grp = g243)

fit_transfer <- function(train, valid, genes, tag) {
  tr <- t(train$expr[genes, , drop = FALSE])
  y_tr <- as.numeric(train$grp != "Control")
  m_tr <- colMeans(tr); s_tr <- apply(tr, 2, sd)
  tr_z <- sweep(sweep(tr, 2, m_tr), 2, s_tr, "/")
  colnames(tr_z) <- genes
  d_tr <- data.frame(y = y_tr, tr_z)
  f <- lrm(y ~ ., data = d_tr, x = TRUE, y = TRUE, penalty = 1)

  va <- t(valid$expr[genes, , drop = FALSE])
  y_va <- as.numeric(valid$grp != "Control")
  va_z <- sweep(sweep(va, 2, m_tr), 2, s_tr, "/"); colnames(va_z) <- genes
  pred_va <- plogis(predict(f, data.frame(va_z)))
  pred_tr <- plogis(predict(f, d_tr))
  roc_tr <- roc(y_tr, pred_tr, quiet = TRUE)
  roc_va <- roc(y_va, pred_va, quiet = TRUE)
  list(tag = tag, n_genes = length(genes),
       auc_tr = as.numeric(auc(roc_tr)), auc_va = as.numeric(auc(roc_va)),
       coef = coef(f))
}

cat("组水平: GSE76293 =", paste(table(d76293$grp), collapse="/"),
    "| GSE243066 =", paste(table(d243066$grp), collapse="/"), "\n")

r6 <- fit_transfer(d76293, d243066, HUB6, "six-gene (with CD177)")
r5 <- fit_transfer(d76293, d243066, HUB5, "five-gene (CD177 removed)")

cat("六基因: 训练 AUC =", sprintf("%.3f", r6$auc_tr), "| 转移 AUC =", sprintf("%.3f", r6$auc_va), "\n")
cat("五基因(去CD177): 训练 AUC =", sprintf("%.3f", r5$auc_tr), "| 转移 AUC =", sprintf("%.3f", r5$auc_va), "\n")

res <- data.frame(
  model = c(r6$tag, r5$tag),
  n_genes = c(r6$n_genes, r5$n_genes),
  train_AUC = round(c(r6$auc_tr, r5$auc_tr), 4),
  transfer_AUC_GSE243066 = round(c(r6$auc_va, r5$auc_va), 4)
)
write.csv(res, file.path(OUT, "cd177_ablation_transfer_AUC.csv"), row.names = FALSE)
cat("结果已写出:", file.path(OUT, "cd177_ablation_transfer_AUC.csv"), "\n")
