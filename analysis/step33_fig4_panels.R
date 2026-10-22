#!/usr/bin/env Rscript
# Figure 4 panels (rms nomogram + ggplot2 NPG):
#   A,B: nomograms (ARDS, sepsis) - rms plot, vector PDF, Helvetica
#   C,D: ROC + calibration + DCA per disease (3 sub-panels each)
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
suppressPackageStartupMessages({ library(rms); library(pROC); library(dplyr) })

BASE <- PROJ
hub <- HUB6

## data loaders (same as step2e)
load_expr <- function(ds) {
  expr <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_gene_expr.rds")))
  grp  <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_grp.rds")))
  list(expr = expr, grp = grp$group)
}
d76293 <- load_expr("GSE76293")
d69528 <- load_expr("GSE69528")

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

DD <- file.path(BASE, "00_rawdata/validation/GSE100159")
lines <- readLines(file.path(DD, "GPL6884.annot"))
b <- which(lines == "!platform_table_begin")
hdr <- strsplit(lines[b + 1], "\t")[[1]]
ann <- read.delim(file.path(DD, "GPL6884.annot"), header = FALSE, stringsAsFactors = FALSE,
                  skip = b, quote = "", comment.char = "", check.names = FALSE)
colnames(ann) <- hdr
ann <- ann[ann[["Gene symbol"]] != "", c("ID", "Gene symbol")]; colnames(ann) <- c("probe", "gene")
sm <- readLines(file.path(DD, "GSE100159_series_matrix.txt"))
tb <- which(sm == "!series_matrix_table_begin"); te <- which(sm == "!series_matrix_table_end")
hdr2 <- gsub('"', '', strsplit(sm[tb + 1], "\t")[[1]])
m4 <- read.delim(textConnection(sm[(tb + 2):(te - 1)]), header = FALSE, sep = "\t",
                 stringsAsFactors = FALSE, check.names = FALSE)
colnames(m4) <- hdr2; rownames(m4) <- m4[, 1]
e4 <- as.matrix(m4[, -1, drop = FALSE]); storage.mode(e4) <- "double"
i <- which(grepl("sample class:", sm, fixed = TRUE))[1]
cls <- strsplit(gsub("!Sample_characteristics_ch1\t|\"", "", sm[i]), "\t")[[1]]
cls <- cls[cls != ""][1:ncol(e4)]
g100 <- ifelse(cls %in% c("sample class: melioidosis", "sample class: other infection"), "Sepsis",
        ifelse(cls == "sample class: Control", "Control", "Exclude"))
keep <- g100 %in% c("Sepsis", "Control"); e4 <- e4[, keep, drop = FALSE]; g100 <- g100[keep]
common <- intersect(rownames(e4), ann$probe); e4 <- e4[common, , drop = FALSE]
gmap <- setNames(ann$gene[match(common, ann$probe)], common)
gmean <- rowMeans(e4)
dfp <- data.frame(probe = common, gene = gmap, mean = gmean, stringsAsFactors = FALSE)
dfp <- dfp[order(dfp$mean, decreasing = TRUE), ]; dfp <- dfp[!duplicated(dfp$gene), ]
e4 <- e4[dfp$probe, , drop = FALSE]; rownames(e4) <- dfp$gene
d100159 <- list(expr = e4, grp = g100)

dca_nb <- function(pred, y, thr = seq(0.01, 0.99, 0.01)) {
  n <- length(y); res <- data.frame(threshold = thr, model = NA, all = NA, none = 0)
  for (k in seq_along(thr)) {
    pt <- thr[k]
    tp <- sum(pred >= pt & y == 1); fp <- sum(pred >= pt & y == 0)
    res$model[k] <- (tp - fp * pt / (1 - pt)) / n
    res$all[k]   <- (sum(y == 1) - sum(y == 0) * pt / (1 - pt)) / n
  }
  res
}

build_model <- function(train, valid, tag, panel_nom) {
  g <- intersect(hub, rownames(train$expr))
  tr <- t(train$expr[g, , drop = FALSE]); y_tr <- as.numeric(train$grp != "Control")
  m_tr <- colMeans(tr); s_tr <- apply(tr, 2, sd)
  tr_z <- sweep(sweep(tr, 2, m_tr), 2, s_tr, "/"); colnames(tr_z) <- g
  d_tr <- data.frame(y = y_tr, tr_z)
  dd <- datadist(d_tr); assign("dd", dd, envir = .GlobalEnv); options(datadist = "dd")
  f <- lrm(y ~ ., data = d_tr, x = TRUE, y = TRUE, penalty = 1)
  nom <- nomogram(f, fun = plogis, funlabel = "Disease Probability",
                  lp = FALSE, fun.at = c(0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95))
  grDevices::pdf(file.path(OUT_MAIN, paste0("Fig4", panel_nom, "_nomogram_", tag, ".pdf")),
                 width = 120/25.4, height = 85/25.4, family = "Helvetica", useDingbats = FALSE)
  par(cex = 0.62)
  plot(nom)
  dev.off()

  pred_tr <- plogis(predict(f, d_tr))
  va <- t(valid$expr[g, , drop = FALSE]); y_va <- as.numeric(valid$grp != "Control")
  va_z <- sweep(sweep(va, 2, m_tr), 2, s_tr, "/"); colnames(va_z) <- g
  pred_va <- plogis(predict(f, data.frame(va_z)))
  roc_tr <- roc(y_tr, pred_tr, quiet = TRUE)
  roc_va <- roc(y_va, pred_va, quiet = TRUE)
  cal_tr <- calibrate(f, method = "boot", B = 200)
  list(f = f, pred_tr = pred_tr, pred_va = pred_va, y_tr = y_tr, y_va = y_va,
       roc_tr = roc_tr, roc_va = roc_va, cal_tr = cal_tr,
       dca_tr = dca_nb(pred_tr, y_tr), dca_va = dca_nb(pred_va, y_va))
}

cat("building ARDS model...\n")
m_ards <- build_model(d76293, d243066, "ARDS", "A")
cat("building sepsis model...\n")
m_sepsis <- build_model(d69528, d100159, "Sepsis", "B")
cat(sprintf("ARDS  train AUC=%.3f  valid AUC=%.3f\n", m_ards$roc_tr$auc, m_ards$roc_va$auc))
cat(sprintf("Sepsis train AUC=%.3f  valid AUC=%.3f\n", m_sepsis$roc_tr$auc, m_sepsis$roc_va$auc))

roc_plot <- function(m, ttl) {
  d <- rbind(
    data.frame(FPR = 1 - m$roc_tr$specificities, TPR = m$roc_tr$sensitivities,
               set = paste0("Train (AUC = ", sprintf("%.3f", m$roc_tr$auc), ")")),
    data.frame(FPR = 1 - m$roc_va$specificities, TPR = m$roc_va$sensitivities,
               set = paste0("Valid (AUC = ", sprintf("%.3f", m$roc_va$auc), ")")))
  ggplot(d, aes(FPR, TPR, color = set)) + geom_line(linewidth = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    scale_color_npg() + coord_equal() +
    labs(x = "1 - Specificity", y = "Sensitivity", title = ttl, color = NULL) +
    theme_ajrccm() + theme(legend.position = c(0.62, 0.18), legend.text = element_text(size = 5.5),
                           legend.background = element_blank())
}
cal_plot <- function(m, ttl) {
  cald <- data.frame(predicted = m$cal_tr[, 1], actual = m$cal_tr[, 2], set = "Train (boot)")
  va <- data.frame(pred = m$pred_va, y = m$y_va); va <- va[order(va$pred), ]
  va$bin <- cut(va$pred, breaks = 5, include.lowest = TRUE, labels = FALSE)
  agg <- aggregate(cbind(pred = pred, y = y) ~ bin, data = va, FUN = mean)
  vd <- data.frame(predicted = agg$pred, actual = agg$y, set = "Valid (bin)")
  d <- rbind(cald, vd)
  ggplot(d, aes(predicted, actual, color = set)) +
    geom_line(linewidth = 0.6) + geom_point(size = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    scale_color_npg() + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Predicted probability", y = "Observed probability", title = ttl, color = NULL) +
    theme_ajrccm() + theme(legend.position = "bottom", legend.text = element_text(size = 5.5))
}
dca_plot <- function(m, ttl) {
  d <- rbind(
    data.frame(threshold = m$dca_tr$threshold, nb = m$dca_tr$model, set = "Model (train)"),
    data.frame(threshold = m$dca_va$threshold, nb = m$dca_va$model, set = "Model (valid)"),
    data.frame(threshold = m$dca_tr$threshold, nb = m$dca_tr$all, set = "Treat all"),
    data.frame(threshold = m$dca_tr$threshold, nb = m$dca_tr$none, set = "Treat none"))
  d$set <- factor(d$set, levels = c("Model (train)", "Model (valid)", "Treat all", "Treat none"))
  ggplot(d, aes(threshold, nb, color = set)) + geom_line(linewidth = 0.6) +
    scale_color_manual(values = c(pal_npg(2), "grey40", "grey75")) +
    coord_cartesian(xlim = c(0, 1), ylim = c(-0.1, NA)) +
    labs(x = "Threshold probability", y = "Net benefit", title = ttl, color = NULL) +
    theme_ajrccm() + theme(legend.position = "bottom", legend.text = element_text(size = 5.5))
}

combo <- function(m, tag, letter) {
  p <- roc_plot(m, paste0(letter, "1  ", tag, " ROC")) +
       cal_plot(m, paste0(letter, "2  ", tag, " calibration")) +
       dca_plot(m, paste0(letter, "3  ", tag, " decision curve")) +
       plot_layout(nrow = 1, guides = "collect") &
       guides(color = guide_legend(nrow = 2, byrow = TRUE))
  save_pdf(p, file.path(OUT_MAIN, paste0("Fig4", letter, "_ROC_cal_DCA_", tag, ".pdf")), 178, 78)
}
combo(m_ards, "ARDS", "C")
combo(m_sepsis, "sepsis", "D")
cat("Figure 4 panels done.\n")
