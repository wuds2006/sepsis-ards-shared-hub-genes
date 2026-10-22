# ============================================================
# Step 75: M1 敏感性分析 —— 七基因(含 MUC1) vs 六基因模型
# 目的: 量化 MUC1 剔除(验证集筛选步骤)对验证队列性能估计的影响
# 输出: 06_external_validation/muc1_sensitivity.csv
# ============================================================
suppressPackageStartupMessages({ library(rms); library(pROC) })
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "06_external_validation")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

hub <- c("CDK5RAP2","CD177","GRB10","UPP1","EXOSC4","VNN1")

load_expr <- function(ds) {
  expr <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_gene_expr.rds")))
  grp  <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_grp.rds")))
  list(expr=expr, grp=grp$group)
}
d76293 <- load_expr("GSE76293")
d69528 <- load_expr("GSE69528")

dd <- file.path(BASE, "00_rawdata/validation/GSE243066")
files <- list.files(dd, pattern="raw_count.txt.gz$", full.names=TRUE)
mat <- NULL
for (f in files) {
  x <- read.delim(gzfile(f), header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
  colnames(x)[2] <- basename(f)
  mat <- if (is.null(mat)) x else merge(mat, x, by="Gene_symbol", all=TRUE)
}
rownames(mat) <- mat$Gene_symbol; mat <- as.matrix(mat[, -1, drop=FALSE])
mat[is.na(mat)] <- 0; storage.mode(mat) <- "double"
sn <- colnames(mat)
g243 <- ifelse(grepl("_A[0-9]", sn), "ARDS", ifelse(grepl("_H[0-9]", sn), "Control", "Exclude"))
cpm <- t(t(mat) / colSums(mat) * 1e6); lc <- log2(cpm + 1)
keep <- g243 %in% c("ARDS","Control"); lc <- lc[, keep, drop=FALSE]; g243 <- g243[keep]
d243066 <- list(expr=lc, grp=g243)

DD <- file.path(BASE, "00_rawdata/validation/GSE100159")
lines <- readLines(file.path(DD, "GPL6884.annot"))
b <- which(lines == "!platform_table_begin")
hdr <- strsplit(lines[b+1], "\t")[[1]]
ann <- read.delim(file.path(DD, "GPL6884.annot"), header=FALSE, stringsAsFactors=FALSE,
                  skip=b, quote="", comment.char="", check.names=FALSE)
colnames(ann) <- hdr
ann <- ann[ann[["Gene symbol"]] != "", c("ID","Gene symbol")]; colnames(ann) <- c("probe","gene")
sm <- readLines(file.path(DD, "GSE100159_series_matrix.txt"))
tb <- which(sm == "!series_matrix_table_begin"); te <- which(sm == "!series_matrix_table_end")
hdr2 <- gsub('"','',strsplit(sm[tb+1], "\t")[[1]])
m4 <- read.delim(textConnection(sm[(tb+2):(te-1)]), header=FALSE, sep="\t",
                 stringsAsFactors=FALSE, check.names=FALSE)
colnames(m4) <- hdr2; rownames(m4) <- m4[,1]
e4 <- as.matrix(m4[, -1, drop=FALSE]); storage.mode(e4) <- "double"
i <- which(grepl("sample class:", sm, fixed=TRUE))[1]
cls <- strsplit(gsub("!Sample_characteristics_ch1\t|\"","",sm[i]),"\t")[[1]]
cls <- cls[cls != ""][1:ncol(e4)]
g100 <- ifelse(cls %in% c("sample class: melioidosis","sample class: other infection"), "Sepsis",
        ifelse(cls == "sample class: Control", "Control", "Exclude"))
keep <- g100 %in% c("Sepsis","Control"); e4 <- e4[, keep, drop=FALSE]; g100 <- g100[keep]
common <- intersect(rownames(e4), ann$probe); e4 <- e4[common,,drop=FALSE]
gmap <- setNames(ann$gene[match(common, ann$probe)], common)
gmean <- rowMeans(e4)
df <- data.frame(probe=common, gene=gmap, mean=gmean, stringsAsFactors=FALSE)
df <- df[order(df$mean, decreasing=TRUE), ]; df <- df[!duplicated(df$gene), ]
e4 <- e4[df$probe,,drop=FALSE]; rownames(e4) <- df$gene[match(df$probe, df$probe)]
d100159 <- list(expr=e4, grp=g100)

fit_eval <- function(train, valid, genes, label) {
  g <- intersect(genes, rownames(train$expr))
  g <- intersect(g, rownames(valid$expr))
  tr <- t(train$expr[g, , drop=FALSE]); y_tr <- as.numeric(train$grp != "Control")
  m_tr <- colMeans(tr); s_tr <- apply(tr, 2, sd)
  tr_z <- sweep(sweep(tr, 2, m_tr), 2, s_tr, "/"); colnames(tr_z) <- g
  d_tr <- data.frame(y = y_tr, tr_z)
  dd <- datadist(d_tr); assign("dd", dd, envir = .GlobalEnv); options(datadist = "dd")
  f <- lrm(y ~ ., data = d_tr, x=TRUE, y=TRUE, penalty = 1)
  pred_tr <- plogis(predict(f, d_tr))
  va <- t(valid$expr[g, , drop=FALSE]); y_va <- as.numeric(valid$grp != "Control")
  va_z <- sweep(sweep(va, 2, m_tr), 2, s_tr, "/"); colnames(va_z) <- g
  pred_va <- plogis(predict(f, data.frame(va_z)))
  ci_tr <- ci.auc(roc(y_tr, pred_tr, quiet=TRUE))  # pROC 返回 (lower, estimate, upper)
  ci_va <- ci.auc(roc(y_va, pred_va, quiet=TRUE))
  data.frame(model=label, genes=length(g),
             train_auc=round(as.numeric(ci_tr[2]),3),
             valid_auc=round(as.numeric(ci_va[2]),3),
             valid_lo=round(as.numeric(ci_va[1]),3),
             valid_hi=round(as.numeric(ci_va[3]),3))
}

res <- rbind(
  fit_eval(d76293,  d243066, hub,               "ARDS 6-gene (final)"),
  fit_eval(d76293,  d243066, c(hub, "MUC1"),    "ARDS 7-gene (incl MUC1)"),
  fit_eval(d69528,  d100159, hub,               "Sepsis 6-gene (final)"),
  fit_eval(d69528,  d100159, c(hub, "MUC1"),    "Sepsis 7-gene (incl MUC1)")
)
write.csv(res, file.path(OUT, "muc1_sensitivity.csv"), row.names=FALSE)
print(res)
cat("\nMUC1 presence: train ARDS:", "MUC1" %in% rownames(d76293$expr),
    "| valid ARDS:", "MUC1" %in% rownames(d243066$expr),
    "| train sepsis:", "MUC1" %in% rownames(d69528$expr),
    "| valid sepsis:", "MUC1" %in% rownames(d100159$expr), "\n")
