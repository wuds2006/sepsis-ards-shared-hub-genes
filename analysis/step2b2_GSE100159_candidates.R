# ============================================================
# Step 2b 补充2: GSE100159 (Sepsis 外验证) 对 ARDS ML 全部候选
# 基因 (9 个 key genes) + Sepsis ML 候选基因计算单基因 AUC
# 用途: 核心基因筛选口径 (ARDS 候选 ∩ 外验证 AUC>0.9) 的可溯源性闭环
# 输出: 06_external_validation/GSE100159_AUC_CI_candidates.csv
# ============================================================
suppressPackageStartupMessages({ library(pROC) })

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
DD <- file.path(BASE, "00_rawdata/validation/GSE100159")

ml_ards  <- readLines(file.path(BASE, "05_ML_SHAP", "ML_key_genes.txt"))
ml_seps  <- readLines(file.path(BASE, "05_ML_SHAP", "sepsis", "ML_key_genes_sepsis.txt"))
genes <- unique(c(ml_ards, ml_seps))

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

# ---- 样本分类: 病例(melioidosis+other infection) vs Control ----
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

# ---- 探针->基因(取平均表达最高探针) ----
keep <- grp %in% c("Sepsis", "Control")
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

cat("候选基因:", paste(genes, collapse=", "), "\n")
avail <- intersect(genes, rownames(expr2))
cat("可验证:", length(avail), "/", length(genes), "\n")

res <- list()
for (g in avail) {
  ro <- roc(grp2, expr2[g, ], levels=c("Control", "Sepsis"), direction="<", quiet=TRUE)
  ci <- ci.auc(ro)
  src <- ifelse(g %in% ml_ards, "ARDS_ML", ifelse(g %in% ml_seps, "Sepsis_ML", "both"))
  res[[g]] <- data.frame(gene=g, source=src, AUC=as.numeric(ro$auc),
                         CI_lower=ci[1], CI_upper=ci[3])
}
r <- do.call(rbind, res)
r <- r[order(-r$AUC), ]
print(r, row.names=FALSE)
write.csv(r, file.path(BASE, "06_external_validation", "GSE100159_AUC_CI_candidates.csv"), row.names=FALSE)

# ---- 核心基因判定: ARDS 候选 ∩ (GSE100159 或 GSE243066 AUC>0.9) ----
val2 <- read.csv(file.path(BASE, "06_external_validation", "GSE243066_validation_AUC.csv"))
strong <- unique(c(r$gene[r$source != "Sepsis_ML" & r$AUC > 0.9],
                   as.character(val2$gene[!is.na(val2$AUC) & val2$AUC > 0.9])))
core <- intersect(ml_ards, strong)
cat("\n外验证 AUC>0.9 基因(合并两队列):", paste(sort(strong), collapse=", "), "\n")
cat("核心基因 (ARDS ML 候选 ∩ 外验证>0.9):", paste(core, collapse=", "), "\n")
