# ============================================================
# Step 2b 补充3: GSE243066 (ARDS 外验证) 对 Sepsis ML 候选基因
# 计算单基因 AUC —— 用于「双队列候选合并 → 跨病种外验证」口径核验
# 输出: 06_external_validation/GSE243066_AUC_sepsis_candidates.csv
# ============================================================
suppressPackageStartupMessages({ library(pROC) })

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
ml_ards <- readLines(file.path(BASE, "05_ML_SHAP", "ML_key_genes.txt"))
ml_seps <- readLines(file.path(BASE, "05_ML_SHAP", "sepsis", "ML_key_genes_sepsis.txt"))
genes <- setdiff(unique(c(ml_ards, ml_seps)), ml_ards)  # 只补 ARDS 候选之外的基因

d <- file.path(BASE, "00_rawdata/validation/GSE243066")
files <- list.files(d, pattern = "raw_count.txt.gz$", full.names = TRUE)
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

sample_names <- colnames(mat)
grp <- ifelse(grepl("_A[0-9]", sample_names), "ARDS",
       ifelse(grepl("_H[0-9]", sample_names), "Control", "Exclude"))

libsize <- colSums(mat)
cpm <- t(t(mat) / libsize * 1e6)
logcpm <- log2(cpm + 1)
keep <- grp %in% c("ARDS", "Control")
logcpm <- logcpm[, keep, drop = FALSE]; grp <- grp[keep]

cat("待验证基因:", paste(genes, collapse=", "), "\n")
avail <- intersect(genes, rownames(logcpm))
res <- list()
for (g in avail) {
  ro <- roc(grp, logcpm[g, ], levels = c("Control","ARDS"), direction = "<", quiet = TRUE)
  res[[g]] <- data.frame(gene = g, source = "Sepsis_ML", AUC = as.numeric(ro$auc), cohort = "GSE243066_ARDS")
}
r <- do.call(rbind, res)
r <- r[order(-r$AUC), ]
print(r, row.names = FALSE)
write.csv(r, file.path(BASE, "06_external_validation", "GSE243066_AUC_sepsis_candidates.csv"), row.names = FALSE)

# ---- 汇总: 双队列候选合并 → 跨病种外验证口径核验 ----
# 规则: AUC > 0.9 于任一外验证队列 且 AUC > 0.7 于另一队列
g100159 <- read.csv(file.path(BASE, "06_external_validation", "GSE100159_AUC_CI_candidates.csv"))
g243066 <- rbind(
  data.frame(gene = c("CDK5RAP2","CD177","GRB10","UPP1","EXOSC4","VNN1","MUC1","TCTN1","DYRK2"),
             AUC = c(0.782352941176471, 0.998039215686275, 0.970588235294118,
                     0.935294117647059, 0.831372549019608, 0.982352941176471,
                     0.574509803921569, 0.119607843137255, 0.00392156862745098)),
  data.frame(gene = r$gene, AUC = r$AUC))
a1 <- setNames(g100159$AUC, g100159$gene)  # GSE100159 (Sepsis)
a2 <- setNames(g243066$AUC, g243066$gene)  # GSE243066 (ARDS)
allc <- unique(c(ml_ards, ml_seps))
core <- c()
for (g in allc) {
  if (is.na(a1[g]) || is.na(a2[g])) next
  if ((a1[g] > 0.9 && a2[g] > 0.7) || (a2[g] > 0.9 && a1[g] > 0.7)) core <- c(core, g)
}
cat("\n合并候选 (", length(allc), "):", paste(allc, collapse=", "), "\n")
cat("跨病种外验证口径核心基因 (", length(core), "):", paste(sort(core), collapse=", "), "\n")
for (g in allc) cat(sprintf("  %-9s GSE100159=%.3f  GSE243066=%.3f\n", g, a1[[g]], a2[[g]]))
