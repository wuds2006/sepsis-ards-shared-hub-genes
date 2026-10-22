# ============================================================
# Step 2b —— GSE243066 (ARDS 外验, RNA-seq 计数)
# ============================================================
suppressPackageStartupMessages({
  library(pROC); library(ggplot2)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "06_external_validation")
key_genes <- readLines(file.path(BASE, "05_ML_SHAP", "ML_key_genes.txt"))

d <- file.path(BASE, "00_rawdata/validation/GSE243066")
files <- list.files(d, pattern = "raw_count.txt.gz$", full.names = TRUE)
cat("样本数:", length(files), "\n")

# 读取合并
mat <- NULL
for (f in files) {
  x <- read.delim(gzfile(f), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(x)[2] <- basename(f)
  if (is.null(mat)) mat <- x else mat <- merge(mat, x, by = "Gene_symbol", all = TRUE)
}
rownames(mat) <- mat$Gene_symbol
mat <- as.matrix(mat[, -1, drop = FALSE])
mat[is.na(mat)] <- 0
storage.mode(mat) <- "double"

# 分组: A* = ARDS, H* = Healthy
sample_names <- colnames(mat)
grp <- ifelse(grepl("_A[0-9]", sample_names), "ARDS",
       ifelse(grepl("_H[0-9]", sample_names), "Control", "Exclude"))
cat("分组: ARDS=", sum(grp=="ARDS"), " Control=", sum(grp=="Control"), "\n")

# log2 标准化 (CPM + 1)
libsize <- colSums(mat)
cpm <- t(t(mat) / libsize * 1e6)
logcpm <- log2(cpm + 1)

keep <- grp %in% c("ARDS", "Control")
logcpm <- logcpm[, keep, drop = FALSE]; grp <- grp[keep]

# 单基因 + 联合 ROC
avail <- intersect(key_genes, rownames(logcpm))
cat("可验证关键基因:", length(avail), "/", length(key_genes), "\n")
res <- list()
for (g in avail) {
  ro <- roc(grp, logcpm[g, ], levels = c("Control","ARDS"), direction = "<", quiet = TRUE)
  res[[g]] <- data.frame(gene = g, type = "single", AUC = as.numeric(ro$auc), cohort = "GSE243066_ARDS")
}
if (length(avail) >= 3) {
  X <- t(logcpm[avail, , drop = FALSE])
  yy <- as.numeric(grp == "ARDS")
  m <- glm(yy ~ X, family = binomial)
  score <- predict(m, type = "response")
  ro <- roc(yy, score, quiet = TRUE)
  res[["signature"]] <- data.frame(gene = "Combined_signature", type = "signature",
                                   AUC = as.numeric(ro$auc), cohort = "GSE243066_ARDS")
  saveRDS(list(y = yy, score = score), file.path(OUT, "GSE243066_sig.rds"))
}
r2 <- do.call(rbind, res); r2 <- r2[order(r2$AUC, decreasing = TRUE), ]
write.csv(r2, file.path(OUT, "GSE243066_validation_AUC.csv"), row.names = FALSE)
print(r2)

# 合并两个队列结果 (GSE100159 Sepsis 候选基因 + 联合签名)
r1c <- read.csv(file.path(OUT, "GSE100159_AUC_CI_candidates.csv"))
r1s <- read.csv(file.path(OUT, "GSE100159_AUC_CI.csv"))
r1 <- data.frame(
  gene   = c("Combined_signature", r1c$gene),
  type   = c("signature", rep("single", nrow(r1c))),
  AUC    = c(r1s$AUC[r1s$gene == "Combined_signature"], r1c$AUC),
  cohort = "GSE100159_Sepsis",
  stringsAsFactors = FALSE)
r1 <- r1[order(r1$AUC, decreasing = TRUE), ]
all_res <- rbind(r1, r2)
write.csv(all_res, file.path(OUT, "validation_AUC_combined.csv"), row.names = FALSE)

# ROC 曲线
plt <- list()
for (g in head(r2$gene[r2$type=="single"], 6)) {
  ro <- roc(grp, logcpm[g, ], levels=c("Control","ARDS"), direction="<", quiet=TRUE)
  plt[[g]] <- data.frame(FPR=1-ro$specificities, TPR=ro$sensitivities,
                         gene=paste0(g," (AUC=",sprintf("%.2f",ro$auc),")"))
}
plt <- do.call(rbind, plt)
p <- ggplot(plt, aes(FPR, TPR, color=gene)) + geom_line(size=0.9) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
  labs(x="1 - Specificity", y="Sensitivity", title="GSE243066 ARDS: Key Gene ROC") +
  theme_pub()
save_fig(p, "GSE243066_ROC", w = 8, h = 6, dir = OUT)

cat("\n===== Step2b (GSE243066) 完成 =====\n")
