# ============================================================
# Step 2b 补充: 计算外验证 AUC 的 95% 置信区间 (DeLong 法)
# 外验队列: GSE100159 (Sepsis) + GSE243066 (ARDS)
# 输出: 06_external_validation/AUC_CI.csv
# ============================================================
suppressPackageStartupMessages({ library(pROC) })
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
DD   <- file.path(BASE, "00_rawdata/validation/GSE100159")
OUT  <- file.path(BASE, "06_external_validation")

hub_genes <- readLines(file.path(BASE, "07_hub_immune", "hub_genes.txt"))
cat("Hub genes:", paste(hub_genes, collapse=", "), "\n\n")

# ---- GSE100159 (Sepsis; GPL6884, 33 病例 + 12 对照) ----
lines <- readLines(file.path(DD, "GPL6884.annot"))
b <- which(lines == "!platform_table_begin")
hdr <- strsplit(lines[b+1], "\t")[[1]]
ann <- read.delim(file.path(DD, "GPL6884.annot"), header=FALSE, stringsAsFactors=FALSE,
                  skip=b, quote="", comment.char="", check.names=FALSE)
colnames(ann) <- hdr
ann <- ann[ann[["Gene symbol"]] != "", c("ID", "Gene symbol")]
colnames(ann) <- c("probe", "gene")

sm <- readLines(file.path(DD, "GSE100159_series_matrix.txt"))
tb <- which(sm == "!series_matrix_table_begin")
te <- which(sm == "!series_matrix_table_end")
hdr2 <- gsub('"', '', strsplit(sm[tb+1], "\t")[[1]])
mmat <- read.delim(textConnection(sm[(tb+2):(te-1)]), header=FALSE, sep="\t",
                   stringsAsFactors=FALSE, check.names=FALSE)
colnames(mmat) <- hdr2
rownames(mmat) <- mmat[, 1]
expr <- as.matrix(mmat[, -1, drop=FALSE])
storage.mode(expr) <- "double"

i <- which(grepl("sample class:", sm, fixed=TRUE))[1]
cls <- strsplit(gsub("!Sample_characteristics_ch1\t|\"", "", sm[i]), "\t")[[1]]
cls <- cls[cls != ""]
samples <- colnames(expr)
cls <- cls[1:length(samples)]
names(cls) <- samples
case <- cls %in% c("sample class: melioidosis", "sample class: other infection")
ctrl <- cls == "sample class: Control"
grp <- rep(NA, length(cls)); grp[case] <- "Sepsis"; grp[ctrl] <- "Control"
names(grp) <- samples

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
rownames(expr2) <- df$gene

avail <- intersect(hub_genes, rownames(expr2))
res1 <- list()
for (g in avail) {
  ro <- roc(grp2, expr2[g, ], levels=c("Control","Sepsis"), direction="<", quiet=TRUE)
  ci <- ci.auc(ro)
  res1[[g]] <- data.frame(gene=g, cohort="GSE100159_Sepsis",
                          AUC=as.numeric(ro$auc), CI_lower=ci[1], CI_upper=ci[3])
}
X <- t(expr2[avail, , drop=FALSE]); yy <- as.numeric(grp2=="Sepsis")
m <- glm(yy ~ X, family=binomial); score <- predict(m, type="response")
ro <- roc(yy, score, quiet=TRUE); ci <- ci.auc(ro)
res1[["signature"]] <- data.frame(gene="Combined_signature", cohort="GSE100159_Sepsis",
                                  AUC=as.numeric(ro$auc), CI_lower=ci[1], CI_upper=ci[3])

# ---- GSE243066 (ARDS) ----
d <- file.path(BASE, "00_rawdata/validation/GSE243066")
files <- list.files(d, pattern="raw_count.txt.gz$", full.names=TRUE)
mat <- NULL
for (f in files) {
  x <- read.delim(gzfile(f), header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
  colnames(x)[2] <- basename(f)
  mat <- if (is.null(mat)) x else merge(mat, x, by="Gene_symbol", all=TRUE)
}
rownames(mat) <- mat$Gene_symbol; mat <- as.matrix(mat[, -1, drop=FALSE])
mat[is.na(mat)] <- 0; storage.mode(mat) <- "double"
sample_names <- colnames(mat)
grp2 <- ifelse(grepl("_A[0-9]", sample_names), "ARDS",
        ifelse(grepl("_H[0-9]", sample_names), "Control", "Exclude"))
cpm <- t(t(mat) / colSums(mat) * 1e6); logcpm <- log2(cpm + 1)
keep <- grp2 %in% c("ARDS","Control"); logcpm <- logcpm[, keep, drop=FALSE]; grp2 <- grp2[keep]

avail2 <- intersect(hub_genes, rownames(logcpm))
res2 <- list()
for (g in avail2) {
  ro <- roc(grp2, logcpm[g, ], levels=c("Control","ARDS"), direction="<", quiet=TRUE)
  ci <- ci.auc(ro)
  res2[[g]] <- data.frame(gene=g, cohort="GSE243066_ARDS",
                          AUC=as.numeric(ro$auc), CI_lower=ci[1], CI_upper=ci[3])
}
X <- t(logcpm[avail2, , drop=FALSE]); yy2 <- as.numeric(grp2=="ARDS")
m <- glm(yy2 ~ X, family=binomial); score <- predict(m, type="response")
ro <- roc(yy2, score, quiet=TRUE); ci <- ci.auc(ro)
res2[["signature"]] <- data.frame(gene="Combined_signature", cohort="GSE243066_ARDS",
                                  AUC=as.numeric(ro$auc), CI_lower=ci[1], CI_upper=ci[3])

all_ci <- rbind(do.call(rbind, res1), do.call(rbind, res2))
all_ci <- all_ci[order(all_ci$cohort, -all_ci$AUC), ]
write.csv(all_ci, file.path(OUT, "AUC_CI.csv"), row.names=FALSE)
cat("\n===== AUC 95% CI 计算完成 =====\n")
print(all_ci, row.names=FALSE)
