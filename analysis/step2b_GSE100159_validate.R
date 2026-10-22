# ============================================================
# GSE100159 可行性验证: 6 hub 基因 ROC/AUC (脓毒症 vs 健康)
# 病例 = 18 melioidosis + 15 other infection = 33 (排除 2 Recovery)
# 对照 = 12 Control
# 平台 GPL6884 (Illumina HumanWG-6 v3.0)
# ============================================================
suppressPackageStartupMessages({ library(pROC) })

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
DD <- file.path(BASE, "00_rawdata/validation/GSE100159")
hub <- c("CDK5RAP2","CD177","GRB10","UPP1","EXOSC4","VNN1")

# ---- 读平台注释(探针->基因) ----
lines <- readLines(file.path(DD, "GPL6884.annot"))
b <- which(lines == "!platform_table_begin")
hdr <- strsplit(lines[b+1], "\t")[[1]]
ann <- read.delim(file.path(DD, "GPL6884.annot"), header=FALSE, stringsAsFactors=FALSE,
                  skip=b, quote="", comment.char="", check.names=FALSE)
colnames(ann) <- hdr
ann <- ann[ann[["Gene symbol"]] != "", c("ID", "Gene symbol")]
colnames(ann) <- c("probe", "gene")

# ---- 读 series matrix 表达矩阵 ----
sm <- readLines(file.path(DD, "GSE100159_series_matrix.txt"))
tb <- which(sm == "!series_matrix_table_begin")
te <- which(sm == "!series_matrix_table_end")
hdr2 <- gsub('"', '', strsplit(sm[tb+1], "\t")[[1]])
mat <- read.delim(textConnection(sm[(tb+2):(te-1)]), header=FALSE, sep="\t",
                  stringsAsFactors=FALSE, check.names=FALSE)
colnames(mat) <- hdr2
rownames(mat) <- mat[, 1]
expr <- as.matrix(mat[, -1, drop=FALSE])
storage.mode(expr) <- "double"
cat("表达矩阵: 探针", nrow(expr), "x 样本", ncol(expr), "\n")

# ---- 样本分类 ----
i <- which(grepl("sample class:", sm, fixed=TRUE))[1]
cls <- strsplit(gsub("!Sample_characteristics_ch1\t|\"", "", sm[i]), "\t")[[1]]
cls <- cls[cls != ""]
samples <- colnames(expr)
cls <- cls[1:length(samples)]
names(cls) <- samples
cat("样本分类统计:\n"); print(table(cls))

# ---- 构建分组: 病例(melioidosis+other infection) vs Control ----
case <- cls %in% c("sample class: melioidosis", "sample class: other infection")
ctrl <- cls == "sample class: Control"
grp <- rep(NA, length(cls)); grp[case] <- "Sepsis"; grp[ctrl] <- "Control"
names(grp) <- samples
cat(sprintf("病例: %d, 对照: %d, 排除(Recovery等): %d\n",
            sum(grp=="Sepsis", na.rm=TRUE), sum(grp=="Control", na.rm=TRUE),
            sum(is.na(grp))))

# ---- 探针->基因(取平均表达最高探针) ----
keep <- grp %in% c("Sepsis","Control")
expr2 <- expr[, keep, drop=FALSE]; grp2 <- grp[keep]
common <- intersect(rownames(expr2), ann$probe)
expr2 <- expr2[common, , drop=FALSE]
gmap <- setNames(ann$gene[match(common, ann$probe)], common)
gmean <- rowMeans(expr2)
df <- data.frame(probe=common, gene=gmap, mean=gmean, stringsAsFactors=FALSE)
df <- df[order(df$mean, decreasing=TRUE), ]
df <- df[!duplicated(df$gene), ]
expr2 <- expr2[df$probe, , drop=FALSE]
rownames(expr2) <- df$gene[match(df$probe, df$probe)]
cat("基因数(去重后):", nrow(expr2), "\n")

# ---- 6 hub 基因 ROC/AUC ----
avail <- intersect(hub, rownames(expr2))
cat("\n可验证 hub 基因:", length(avail), "/", length(hub), "\n")
res <- list()
for(g in avail){
  ro <- roc(grp2, expr2[g, ], levels=c("Control","Sepsis"), direction="<", quiet=TRUE)
  ci <- ci.auc(ro)
  res[[g]] <- data.frame(gene=g, AUC=as.numeric(ro$auc),
                         CI_lower=ci[1], CI_upper=ci[3])
}
# 联合签名(逻辑回归)
X <- t(expr2[avail, , drop=FALSE])
yy <- as.numeric(grp2 == "Sepsis")
m <- glm(yy ~ X, family=binomial)
score <- predict(m, type="response")
ro <- roc(yy, score, quiet=TRUE); ci <- ci.auc(ro)
res[["signature"]] <- data.frame(gene="Combined_signature", AUC=as.numeric(ro$auc),
                                 CI_lower=ci[1], CI_upper=ci[3])
r <- do.call(rbind, res)
r <- r[order(-r$AUC), ]
cat("\n===== GSE100159 AUC 结果 =====\n")
print(r, row.names=FALSE)
write.csv(r, file.path(BASE, "06_external_validation", "GSE100159_AUC_CI.csv"), row.names=FALSE)
cat("\n已保存: 06_external_validation/GSE100159_AUC_CI.csv\n")
