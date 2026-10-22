# ============================================================
# Step 2e: 疾病预测模型构建 (lrm + nomogram + ROC + calibration + DCA)
# ARDS 模型: GSE76293(训练) -> GSE243066(验证)
# Sepsis 模型: GSE69528(训练) -> GSE100159(验证)
# ============================================================
suppressPackageStartupMessages({
  library(rms); library(pROC); library(ggplot2); library(patchwork); library(ragg)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "AJRCCM_submission", "07_figure_exports", "Figure8")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

hub <- c("CDK5RAP2","CD177","GRB10","UPP1","EXOSC4","VNN1")

# ---- 数据读取(复用 step2d 逻辑) ----
load_expr <- function(ds) {
  expr <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_gene_expr.rds")))
  grp  <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_grp.rds")))
  list(expr=expr, grp=grp$group)
}
d76293 <- load_expr("GSE76293")   # ARDS 训练
d69528 <- load_expr("GSE69528")   # Sepsis 训练

# GSE243066 ARDS 验证
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

# GSE100159 Sepsis 验证
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
hdr2 <- gsub('"','',strsplit(sm[tb+1],"\t")[[1]])
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

# ---- 手写 DCA 净收益 ----
dca_nb <- function(pred, y, thr=seq(0.01, 0.99, 0.01)) {
  n <- length(y); res <- data.frame(threshold=thr, model=NA, all=NA, none=0)
  for (k in seq_along(thr)) {
    pt <- thr[k]
    tp <- sum(pred >= pt & y==1); fp <- sum(pred >= pt & y==0)
    res$model[k] <- (tp - fp*pt/(1-pt))/n
    res$all[k]    <- (sum(y==1) - sum(y==0)*pt/(1-pt))/n
  }
  res
}

# ---- 构建模型函数 ----
build_model <- function(train, valid, disease_label, tag) {
  # 提取 6 基因，z-score 标准化(训练集 mean/sd)
  g <- intersect(hub, rownames(train$expr))
  tr <- t(train$expr[g, , drop=FALSE]); y_tr <- as.numeric(train$grp != "Control")
  m_tr <- colMeans(tr); s_tr <- apply(tr, 2, sd)
  tr_z <- sweep(sweep(tr, 2, m_tr), 2, s_tr, "/")
  colnames(tr_z) <- g
  d_tr <- data.frame(y = y_tr, tr_z)
  
  dd <- datadist(d_tr)
  assign("dd", dd, envir = .GlobalEnv)
  options(datadist = "dd")
  # penalty=1 (ridge) 避免完全分离时系数爆炸
  f <- lrm(y ~ ., data = d_tr, x=TRUE, y=TRUE, penalty = 1)
  
  # nomogram
  nom <- nomogram(f, fun=plogis, funlabel="Disease Probability",
                  lp=FALSE, fun.at=c(0.05,0.1,0.3,0.5,0.7,0.9,0.95))
  pdf(file.path(OUT, paste0(tag, "_nomogram.pdf")), width=9, height=6.5)
  plot(nom)
  dev.off()
  png(file.path(OUT, paste0(tag, "_nomogram.png")), width=9, height=6.5, units="in", res=300)
  plot(nom)
  dev.off()
  agg_tiff(file.path(OUT, paste0(tag, "_nomogram.tiff")), width=2700, height=1950,
           units="px", res=300, compression="lzw", background="white")
  plot(nom)
  dev.off()
  
  # 训练集预测
  pred_tr <- plogis(predict(f, d_tr))
  # 验证集预测(z-score 用训练集参数)
  va <- t(valid$expr[g, , drop=FALSE]); y_va <- as.numeric(valid$grp != "Control")
  va_z <- sweep(sweep(va, 2, m_tr), 2, s_tr, "/"); colnames(va_z) <- g
  d_va <- data.frame(va_z)
  pred_va <- plogis(predict(f, d_va))
  
  # ROC
  roc_tr <- roc(y_tr, pred_tr, quiet=TRUE)
  roc_va <- roc(y_va, pred_va, quiet=TRUE)
  
  # calibration (训练 bootstrap)
  cal_tr <- calibrate(f, method="boot", B=200)
  # 验证外部校准 val.prob
  # DCA
  dca_tr <- dca_nb(pred_tr, y_tr)
  dca_va <- dca_nb(pred_va, y_va)
  
  list(f=f, g=g, nom=nom, pred_tr=pred_tr, pred_va=pred_va,
       y_tr=y_tr, y_va=y_va, roc_tr=roc_tr, roc_va=roc_va,
       cal_tr=cal_tr, dca_tr=dca_tr, dca_va=dca_va)
}

m_ards  <- build_model(d76293,  d243066, "ARDS",  "ARDS")
m_sepsis<- build_model(d69528,  d100159, "Sepsis","Sepsis")

# ---- ROC 组合图 ----
roc_plot <- function(m, tag) {
  r1 <- data.frame(FPR=1-m$roc_tr$specificities, TPR=m$roc_tr$sensitivities,
                   set=paste0("Training (AUC=", sprintf("%.3f", m$roc_tr$auc), ")"))
  r2 <- data.frame(FPR=1-m$roc_va$specificities, TPR=m$roc_va$sensitivities,
                   set=paste0("Validation (AUC=", sprintf("%.3f", m$roc_va$auc), ")"))
  d <- rbind(r1, r2)
  ggplot(d, aes(FPR, TPR, color=set)) + geom_line(lwd=0.8) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
    labs(x="1 - Specificity", y="Sensitivity", title=paste(tag, "Model"), color=NULL) +
    scale_color_manual(values=c("#6A5C95","#D99AAE")) +
    theme_classic(base_size=9) + theme(legend.position=c(0.7,0.25))
}

# ---- Calibration 组合图 ----
cal_plot <- function(m, tag) {
  # 训练 bootstrap 校准
  cald <- data.frame(predicted=m$cal_tr[,1], actual=m$cal_tr[,2], set="Training")
  # 验证: 用验证预测概率做校准(按预测概率分箱)
  va <- data.frame(pred=m$pred_va, y=m$y_va)
  va <- va[order(va$pred), ]
  nb <- 5; va$bin <- cut(va$pred, breaks=nb, include.lowest=TRUE, labels=FALSE)
  agg <- aggregate(cbind(pred=pred, y=y) ~ bin, data=va, FUN=mean)
  vd <- data.frame(predicted=agg$pred, actual=agg$y, set="Validation")
  d <- rbind(cald, vd)
  ggplot(d, aes(predicted, actual, color=set)) + geom_line(lwd=0.8) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
    labs(x="Predicted Probability", y="Observed Probability", title=paste(tag,"Calibration"), color=NULL) +
    scale_color_manual(values=c("#6A5C95","#D99AAE")) +
    coord_equal() + theme_classic(base_size=9) + theme(legend.position="bottom")
}

# ---- DCA 组合图 ----
dca_plot <- function(m, tag) {
  d1 <- data.frame(threshold=m$dca_tr$threshold, nb=m$dca_tr$model, set="Model (Training)")
  d2 <- data.frame(threshold=m$dca_va$threshold, nb=m$dca_va$model, set="Model (Validation)")
  d3 <- data.frame(threshold=m$dca_tr$threshold, nb=m$dca_tr$all, set="Treat All")
  d4 <- data.frame(threshold=m$dca_tr$threshold, nb=m$dca_tr$none, set="Treat None")
  d <- rbind(d1, d2, d3, d4)
  ggplot(d, aes(threshold, nb, color=set)) + geom_line(lwd=0.8) +
    labs(x="Threshold Probability", y="Net Benefit", title=paste(tag,"Decision Curve"), color=NULL) +
    scale_color_manual(values=c("#6A5C95","#D99AAE","grey40","grey70")) +
    coord_cartesian(xlim=c(0,1)) +
    theme_classic(base_size=9) + theme(legend.position="bottom")
}

# 保存组合图(每个疾病 3 面板: ROC + Calibration + DCA)
save_combo <- function(m, tag) {
  p <- (roc_plot(m, tag) | cal_plot(m, tag) | dca_plot(m, tag))
  ggsave(file.path(OUT, paste0(tag, "_ROC_cal_DCA.pdf")), p, width=15, height=5)
  ggsave(file.path(OUT, paste0(tag, "_ROC_cal_DCA.png")), p, width=15, height=5, dpi=300)
  ggsave(file.path(OUT, paste0(tag, "_ROC_cal_DCA.tiff")), p, width=15, height=5, dpi=300,
         device=ragg::agg_tiff, compression="lzw", bg="white")
}
save_combo(m_ards, "ARDS")
save_combo(m_sepsis, "Sepsis")

# ---- 模型汇总 ----
cat("\n===== 模型汇总 =====\n")
for (m in list(ARDS=m_ards, Sepsis=m_sepsis)) {
  cat("\n---", names(m), "---\n")
  print(m$f)
  cat(sprintf("训练 AUC=%.3f, 验证 AUC=%.3f\n", m$roc_tr$auc, m$roc_va$auc))
}
cat("\n完成，输出到:", OUT, "\n")
