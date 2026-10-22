# ============================================================
# Step 2b 补充4: GSE100159 (Sepsis 外验证) 6 hub 基因 ROC 曲线图
# 替换原 GSE28750_ROC（数据已清除，外验口径统一为 GSE100159 + GSE243066）
# 输出: 06_external_validation/GSE100159_ROC.pdf / .png
# ============================================================
suppressPackageStartupMessages({ library(pROC); library(ggplot2) })

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
DD <- file.path(BASE, "00_rawdata/validation/GSE100159")
OUT <- file.path(BASE, "06_external_validation")
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
cat(sprintf("Sepsis=%d Control=%d Exclude=%d\n",
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
rownames(expr2) <- df$gene

# ---- 6 hub 基因 ROC 曲线 ----
avail <- intersect(hub, rownames(expr2))
cat("可验证 hub 基因:", length(avail), "/", length(hub), "\n")
plt <- list()
for (g in avail) {
  ro <- roc(grp2, expr2[g, ], levels=c("Control","Sepsis"), direction="<", quiet=TRUE)
  plt[[g]] <- data.frame(FPR=1-ro$specificities, TPR=ro$sensitivities,
                         gene=paste0(g," (AUC=",sprintf("%.2f",as.numeric(ro$auc)),")"))
}
plt <- do.call(rbind, plt)
p <- ggplot(plt, aes(FPR, TPR, color=gene)) + geom_line(linewidth=0.9) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
  labs(x="1 - Specificity", y="Sensitivity", title="GSE100159 Sepsis: Key Gene ROC") +
  theme_classic(base_size=12) +
  theme(legend.title=element_blank(), legend.position=c(0.72, 0.22),
        legend.background=element_rect(fill="white", color=NA))
ggsave(file.path(OUT, "GSE100159_ROC.pdf"), p, width=6.5, height=5.5)
ggsave(file.path(OUT, "GSE100159_ROC.png"), p, width=6.5, height=5.5, dpi=300)
cat("已保存: GSE100159_ROC.pdf / .png\n")
