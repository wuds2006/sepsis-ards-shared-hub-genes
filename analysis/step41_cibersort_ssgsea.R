#!/usr/bin/env Rscript
# step41: CIBERSORT immune infiltration (4 bulk datasets) + 6-hub-gene ssGSEA
# Classic CIBERSORT algorithm (Newman et al. Nat Methods 2015), ported from the
# canonical reference implementation (e1071 nu-SVR + preprocessCore QN).
# Outputs (07_hub_immune/cibersort/):
#   <ds>_CIBERSORT_results.csv  sample x 22 fractions + P-value/Correlation/RMSE
#   <ds>_celltype_wilcox.csv    per-celltype Wilcoxon disease vs control
#   <ds>_spearman_cor.csv       22 x 22 Spearman correlation matrix
# Outputs (07_hub_immune/):
#   <ds>_hub6_ssGSEA.csv        per-sample ssGSEA score of 6-gene signature
#   hub6_ssGSEA_stats.csv       Wilcoxon summary across datasets
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
suppressPackageStartupMessages({
  library(e1071); library(preprocessCore); library(parallel); library(GSVA)
})

BASE <- PROJ
CBD  <- file.path(BASE, "07_hub_immune", "cibersort")
IMM  <- file.path(BASE, "07_hub_immune")

## ---------------- classic CIBERSORT core ----------------
CoreAlg <- function(X, y) {
  nus <- c(0.25, 0.5, 0.75)
  out <- lapply(nus, function(nu)
    e1071::svm(X, y, type = "nu-regression", kernel = "linear", nu = nu, scale = FALSE))
  nusvm <- rep(0, 3); corrv <- rep(0, 3)
  for (t in 1:3) {
    weights <- t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights < 0)] <- 0
    w <- if (sum(weights) > 0) weights / sum(weights) else weights * 0
    u <- sweep(X, MARGIN = 2, w, "*")
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt(mean((k - y)^2))
    corrv[t] <- suppressWarnings(cor(k, y))
  }
  mn <- which.min(nusvm)
  q <- t(out[[mn]]$coefs) %*% out[[mn]]$SV
  q[which(q < 0)] <- 0
  w <- if (sum(q) > 0) q / sum(q) else q * 0
  list(w = as.numeric(w), mix_rmse = nusvm[mn], mix_r = corrv[mn])
}

doPerm <- function(perm, X, Y, cores = 4) {
  set.seed(42)  # reproducible permutation null
  Ylist <- as.numeric(Y)
  one_perm <- function(i) {
    yr <- Ylist[sample(length(Ylist), nrow(X))]
    yr <- (yr - mean(yr)) / sd(yr)
    CoreAlg(X, yr)$mix_r
  }
  unlist(mclapply(1:perm, one_perm, mc.cores = cores))
}

cibersort_run <- function(X, Y, perm = 1000, QN = TRUE, antilog_auto = TRUE) {
  X <- X[order(rownames(X)), ]
  Y <- Y[order(rownames(Y)), ]
  # classic anti-log heuristic for log-scale mixtures
  if (antilog_auto && max(Y) < 50) Y <- 2^Y
  if (QN) {
    tmpr <- rownames(Y); tmpc <- colnames(Y)
    Y <- preprocessCore::normalize.quantiles(Y)
    rownames(Y) <- tmpr; colnames(Y) <- tmpc
  }
  # intersect genes
  Y <- Y[rownames(Y) %in% rownames(X), , drop = FALSE]
  X <- X[rownames(X) %in% rownames(Y), , drop = FALSE]
  n_int <- nrow(X)
  # standardize signature
  X <- (X - mean(X)) / sd(as.vector(X))
  # permutation null
  nulldist <- sort(doPerm(perm, X, Y))
  res <- matrix(NA_real_, nrow = ncol(Y), ncol = ncol(X) + 3,
                dimnames = list(colnames(Y), c(colnames(X), "P-value", "Correlation", "RMSE")))
  for (i in seq_len(ncol(Y))) {
    y <- Y[, i]
    y <- (y - mean(y)) / sd(y)
    r <- CoreAlg(X, y)
    res[i, seq_len(ncol(X))] <- r$w
    res[i, "P-value"] <- 1 - (which.min(abs(nulldist - r$mix_r)) / length(nulldist))
    res[i, "Correlation"] <- r$mix_r
    res[i, "RMSE"] <- r$mix_rmse
  }
  list(result = res, n_int_genes = n_int)
}

## ---------------- dataset loaders (ported from step31) ----------------
LM22 <- read.delim(file.path(CBD, "LM22.txt"), check.names = FALSE)
rownames(LM22) <- LM22[[1]]; LM22 <- as.matrix(LM22[, -1, drop = FALSE])
storage.mode(LM22) <- "double"

load_disc <- function(ds) {
  list(expr = readRDS(file.path(BASE, "01_DEG", paste0(ds, "_gene_expr.rds"))),
       grp  = readRDS(file.path(BASE, "01_DEG", paste0(ds, "_grp.rds"))))
}

load_GSE243066 <- function() {
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
  g <- ifelse(grepl("_A[0-9]", sn), "ARDS", ifelse(grepl("_H[0-9]", sn), "Control", "Exclude"))
  cpm <- t(t(mat) / colSums(mat) * 1e6)          # linear CPM for CIBERSORT
  logcpm <- log2(cpm + 1)                        # log scale for ssGSEA
  keep <- g %in% c("ARDS", "Control")
  list(expr = logcpm, expr_lin = cpm[, keep, drop = FALSE], grp = data.frame(group = g[keep]))
}

load_GSE100159 <- function() {
  DD <- file.path(BASE, "00_rawdata/validation/GSE100159")
  lines <- readLines(file.path(DD, "GPL6884.annot"))
  b <- which(lines == "!platform_table_begin")
  hdr <- strsplit(lines[b + 1], "\t")[[1]]
  ann <- read.delim(file.path(DD, "GPL6884.annot"), header = FALSE, stringsAsFactors = FALSE,
                    skip = b, quote = "", comment.char = "", check.names = FALSE)
  colnames(ann) <- hdr
  ann <- ann[ann[["Gene symbol"]] != "", c("ID", "Gene symbol")]
  colnames(ann) <- c("probe", "gene")
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
  g4 <- ifelse(cls %in% c("sample class: melioidosis", "sample class: other infection"), "Sepsis",
        ifelse(cls == "sample class: Control", "Control", "Exclude"))
  keep <- g4 %in% c("Sepsis", "Control"); e4 <- e4[, keep, drop = FALSE]; g4 <- g4[keep]
  common <- intersect(rownames(e4), ann$probe); e4 <- e4[common, , drop = FALSE]
  gmap <- setNames(ann$gene[match(common, ann$probe)], common)
  gmean <- rowMeans(e4)
  dfp <- data.frame(probe = common, gene = gmap, mean = gmean, stringsAsFactors = FALSE)
  dfp <- dfp[order(dfp$mean, decreasing = TRUE), ]; dfp <- dfp[!duplicated(dfp$gene), ]
  e4 <- e4[dfp$probe, , drop = FALSE]; rownames(e4) <- dfp$gene
  list(expr = e4, grp = data.frame(group = g4))
}

d1 <- load_disc("GSE76293")        # log2 microarray, PMN fractions
d2 <- load_disc("GSE69528")        # log2 microarray, whole blood
d3 <- load_GSE243066()             # log2 CPM + linear CPM
d4 <- load_GSE100159()             # linear Illumina

datasets <- list(
  list(id = "GSE76293",  dis = "ARDS",   lin = d1$expr,  grp = d1$grp$group, QN = TRUE),
  list(id = "GSE69528",  dis = "Sepsis", lin = d2$expr,  grp = d2$grp$group, QN = TRUE),
  list(id = "GSE243066", dis = "ARDS",   lin = d3$expr_lin, grp = d3$grp$group, QN = FALSE),
  list(id = "GSE100159", dis = "Sepsis", lin = d4$expr,  grp = d4$grp$group, QN = TRUE)
)

## ---------------- run CIBERSORT per dataset ----------------
celltypes <- colnames(LM22)
summ <- data.frame()
for (ds in datasets) {
  id <- ds$id
  expr_lin <- as.matrix(ds$lin)
  grp <- ds$grp
  names(grp) <- colnames(expr_lin)
  stopifnot(!anyNA(grp))
  cat(sprintf("== %s: %d samples (%d %s / %d Control) ==\n",
              id, length(grp), sum(grp == ds$dis), ds$dis, sum(grp == "Control")))
  run <- cibersort_run(LM22, expr_lin, perm = 1000, QN = ds$QN)
  res <- run$result
  cat(sprintf("   LM22 intersected genes: %d / 547\n", run$n_int_genes))
  outdf <- data.frame(sample = rownames(res), group = grp[rownames(res)], res, check.names = FALSE)
  write.csv(outdf, file.path(CBD, paste0(id, "_CIBERSORT_results.csv")), row.names = FALSE)

  # Wilcoxon per cell type (two-sided, exact=FALSE to avoid large-sample bug)
  dif <- do.call(rbind, lapply(celltypes, function(ct) {
    x <- res[grp == ds$dis, ct]; y0 <- res[grp == "Control", ct]
    p <- tryCatch(wilcox.test(x, y0, exact = FALSE, correct = FALSE)$p.value,
                  error = function(e) NA)
    data.frame(cell = ct, median_disease = median(x), median_control = median(y0),
               mean_disease = mean(x), mean_control = mean(y0),
               log2FC = ifelse(mean(x) > 0 & mean(y0) > 0,
                               log2((mean(x) + 1e-4) / (mean(y0) + 1e-4)), NA),
               P = p)
  }))
  dif$padj <- p.adjust(dif$P, method = "BH")
  write.csv(dif, file.path(CBD, paste0(id, "_celltype_wilcox.csv")), row.names = FALSE)
  sig <- dif[!is.na(dif$P) & dif$P < 0.05, ]
  cat(sprintf("   significant cell types (P<0.05): %d\n", nrow(sig)))
  if (nrow(sig) > 0) print(sig[, c("cell", "log2FC", "P")], row.names = FALSE)

  # Spearman correlation among 22 cell types
  ok <- rownames(res)  # all
  cor22 <- suppressWarnings(cor(res[, celltypes], method = "spearman"))
  write.csv(cor22, file.path(CBD, paste0(id, "_spearman_cor.csv")))

  summ <- rbind(summ, data.frame(dataset = id, disease = ds$dis,
                                  n_disease = sum(grp == ds$dis), n_control = sum(grp == "Control"),
                                  lm22_genes = run$n_int_genes,
                                  median_P = median(res[, "P-value"]),
                                  median_corr = median(res[, "Correlation"]),
                                  sig_cells = nrow(sig)))
}
write.csv(summ, file.path(CBD, "run_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE)

## ---------------- 6-hub-gene ssGSEA per dataset ----------------
HUB6 <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")
log_expr <- list(GSE76293 = d1$expr, GSE69528 = d2$expr,
                 GSE243066 = d3$expr, GSE100159 = d4$expr)
stats <- data.frame()
for (id in names(log_expr)) {
  expr <- as.matrix(log_expr[[id]])
  ds <- datasets[[which(vapply(datasets, function(z) z$id == id, logical(1)))]]
  grp <- ds$grp; names(grp) <- colnames(expr)
  genes <- intersect(HUB6, rownames(expr))
  sc <- gsva(ssgseaParam(exprData = expr, geneSets = list(hub6 = genes), minSize = 2))
  score <- as.numeric(sc[1, ])
  df <- data.frame(sample = colnames(expr), group = grp, score = score)
  write.csv(df, file.path(IMM, paste0(id, "_hub6_ssGSEA.csv")), row.names = FALSE)
  x <- score[grp == ds$dis]; y0 <- score[grp == "Control"]
  p <- wilcox.test(x, y0, exact = FALSE, correct = FALSE)$p.value
  stats <- rbind(stats, data.frame(dataset = id, disease = ds$dis, genes_used = paste(genes, collapse = ","),
                                   n_disease = length(x), n_control = length(y0),
                                   median_disease = median(x), median_control = median(y0),
                                   P = p))
  cat(sprintf("%s hub6 ssGSEA: median %.3f vs %.3f, P = %.3g\n", id, median(x), median(y0), p))
}
write.csv(stats, file.path(IMM, "hub6_ssGSEA_stats.csv"), row.names = FALSE)
cat("step41 done.\n")
