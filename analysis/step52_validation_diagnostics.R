# ============================================================
# Step 52 —— 验证队列诊断性能系统分析（审稿 W5 / 问题3）
# 1) 单基因 AUC + DeLong 95% CI（两验证队列，统一重算）
# 2) 每基因 训练-验证 表达分布重叠系数 OVL（核密度 min 积分）+ SMD
# 3) MetaDE 式合并：六基因单基因效应（Hedges g）跨两验证队列 FE+RE
# 4) ComBat 跨平台批次校正后重评估六基因模型（敏感性，如实报告）
# ============================================================
suppressPackageStartupMessages({
  library(pROC); library(sva); library(rms)
})
BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "06_external_validation", "reviewer_w5")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
HUB6 <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")

# ---------- loaders（镜像 step49 / step2b2） ----------
load_243066 <- function() {
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
  cpm <- t(t(mat) / colSums(mat) * 1e6); lc <- log2(cpm + 1)
  keep <- g %in% c("ARDS", "Control")
  list(expr = lc[, keep, drop = FALSE], grp = g[keep])
}
load_100159 <- function() {
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
  mat <- read.delim(textConnection(sm[(tb + 2):(te - 1)]), header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE, check.names = FALSE)
  colnames(mat) <- hdr2; rownames(mat) <- mat[, 1]
  expr <- as.matrix(mat[, -1, drop = FALSE]); storage.mode(expr) <- "double"
  i <- which(grepl("sample class:", sm, fixed = TRUE))[1]
  cls <- strsplit(gsub("!Sample_characteristics_ch1\t|\"", "", sm[i]), "\t")[[1]]
  cls <- cls[cls != ""]
  samples <- colnames(expr); cls <- cls[1:length(samples)]; names(cls) <- samples
  case <- cls %in% c("sample class: melioidosis", "sample class: other infection")
  ctrl <- cls == "sample class: Control"
  grp <- rep(NA, length(cls)); grp[case] <- "Sepsis"; grp[ctrl] <- "Control"
  keep <- grp %in% c("Sepsis", "Control")
  expr2 <- expr[, keep, drop = FALSE]
  common <- intersect(rownames(expr2), ann$probe)
  expr2 <- expr2[common, , drop = FALSE]
  gmap <- setNames(ann$gene[match(common, ann$probe)], common)
  gmean <- rowMeans(expr2)
  df <- data.frame(probe = common, gene = gmap, mean = gmean, stringsAsFactors = FALSE)
  df <- df[order(df$mean, decreasing = TRUE), ]
  df <- df[!duplicated(df$gene), ]
  e <- expr2[df$probe, , drop = FALSE]
  rownames(e) <- df$gene                      # 行名映射为基因符号
  list(expr = e, grp = grp[keep])
}

d76293  <- list(expr = readRDS(file.path(BASE, "01_DEG", "GSE76293_gene_expr.rds")),
                grp  = readRDS(file.path(BASE, "01_DEG", "GSE76293_grp.rds"))$group)
d69528  <- list(expr = readRDS(file.path(BASE, "01_DEG", "GSE69528_gene_expr.rds")),
                grp  = readRDS(file.path(BASE, "01_DEG", "GSE69528_grp.rds"))$group)
d243066 <- load_243066()
d100159 <- load_100159()
cat("GSE76293:", paste(table(d76293$grp), collapse = "/"),
    "| GSE243066:", paste(table(d243066$grp), collapse = "/"),
    "| GSE69528:", paste(table(d69528$grp), collapse = "/"),
    "| GSE100159:", paste(table(d100159$grp), collapse = "/"), "\n")

# ---------- 1) 单基因 AUC + DeLong CI ----------
auc_ci <- function(expr, grp, genes, tag) {
  y <- as.numeric(grp != "Control")
  do.call(rbind, lapply(genes, function(g) {
    if (!g %in% rownames(expr)) return(NULL)
    v <- as.numeric(expr[g, ])
    ro <- roc(y, v, quiet = TRUE)
    ci <- ci.auc(ro, method = "delong")
    data.frame(cohort = tag, gene = g, AUC = round(as.numeric(auc(ro)), 3),
               ci_low = round(ci[1], 3), ci_high = round(ci[3], 3))
  }))
}
t1 <- rbind(auc_ci(d243066$expr, d243066$grp, HUB6, "GSE243066 (ARDS validation, 34/15)"),
            auc_ci(d100159$expr, d100159$grp, HUB6, "GSE100159 (Sepsis validation, 33/12)"))
write.csv(t1, file.path(OUT, "single_gene_AUC_CI_validation.csv"), row.names = FALSE)
print(t1)

# ---------- 2) 分布重叠与效应量对比（跨平台，队列内 z 标准化后比较） ----------
# OVL: 每基因在各自队列内 z 标准化（全部样本）后，训练-验证标准化分布的重叠
# SMD_shift: 标准化分布的位移；Cohen_d_train / Cohen_d_val: 队列内疾病 vs 对照效应量
ovl <- function(x, y) {
  rng <- range(c(x, y)); grid <- seq(rng[1], rng[2], length.out = 512)
  fx <- density(x, n = 512, from = rng[1], to = rng[2])
  fy <- density(y, n = 512, from = rng[1], to = rng[2])
  f1 <- approx(fx$x, fx$y, grid, rule = 2)$y
  f2 <- approx(fy$x, fy$y, grid, rule = 2)$y
  dx <- grid[2] - grid[1]
  sum(pmin(f1, f2)) * dx
}
zsc <- function(v) (v - mean(v)) / sd(v)
cohend <- function(x, y) (mean(x) - mean(y)) / sqrt((sd(x)^2 + sd(y)^2) / 2)

ovl_tab <- do.call(rbind, lapply(HUB6, function(g) {
  do.call(rbind, lapply(list(
    list(pair = "ARDS (GSE76293->GSE243066)", tr = d76293, va = d243066, lv = "ARDS"),
    list(pair = "Sepsis (GSE69528->GSE100159)", tr = d69528, va = d100159, lv = "Sepsis")), function(cfg) {
    if (!g %in% rownames(cfg$tr$expr) || !g %in% rownames(cfg$va$expr)) return(NULL)
    ztr <- zsc(as.numeric(cfg$tr$expr[g, ])); zva <- zsc(as.numeric(cfg$va$expr[g, ]))
    dD_tr <- cohend(as.numeric(cfg$tr$expr[g, cfg$tr$grp == cfg$lv]),
                    as.numeric(cfg$tr$expr[g, cfg$tr$grp == "Control"]))
    dD_va <- cohend(as.numeric(cfg$va$expr[g, cfg$va$grp == cfg$lv]),
                    as.numeric(cfg$va$expr[g, cfg$va$grp == "Control"]))
    data.frame(pair = cfg$pair, gene = g,
               OVL_standardized = round(ovl(ztr, zva), 3),
               SMD_shift = round(mean(zva) - mean(ztr), 2),
               Cohen_d_train = round(dD_tr, 2), Cohen_d_validation = round(dD_va, 2),
               effect_attenuation = round(dD_va / dD_tr, 2))
  }))
}))
write.csv(ovl_tab, file.path(OUT, "distribution_overlap_train_vs_validation.csv"), row.names = FALSE)
print(ovl_tab)

# ---------- 3) MetaDE 式合并（单基因 Hedges g 跨两验证队列） ----------
hedges_g <- function(x, y) {                       # x = disease, y = control
  n1 <- length(x); n0 <- length(y)
  s  <- sqrt(((n1 - 1) * sd(x)^2 + (n0 - 1) * sd(y)^2) / (n1 + n0 - 2))
  d  <- (mean(x) - mean(y)) / s
  J  <- 1 - 3 / (4 * (n1 + n0) - 9)
  v  <- (n1 + n0) / (n1 * n0) + d^2 / (2 * (n1 + n0))
  c(g = J * d, var = J^2 * v)
}
meta_tab <- do.call(rbind, lapply(HUB6, function(g) {
  e1 <- hedges_g(as.numeric(d243066$expr[g, d243066$grp == "ARDS"]),
                as.numeric(d243066$expr[g, d243066$grp == "Control"]))
  e2 <- hedges_g(as.numeric(d100159$expr[g, d100159$grp == "Sepsis"]),
                as.numeric(d100159$expr[g, d100159$grp == "Control"]))
  gi <- c(e1[1], e2[1]); vi <- c(e1[2], e2[2]); wi <- 1 / vi
  gfe <- sum(wi * gi) / sum(wi); sefe <- sqrt(1 / sum(wi))
  Q <- sum(wi * (gi - gfe)^2); df <- 1
  tau2 <- max(0, (Q - df) / (sum(wi) - sum(wi^2) / sum(wi)))   # DerSimonian-Laird
  wr <- 1 / (vi + tau2); gre <- sum(wr * gi) / sum(wr); sere <- sqrt(1 / sum(wr))
  I2 <- max(0, (Q - df) / Q)
  data.frame(gene = g,
             g_GSE243066 = round(gi[1], 2), g_GSE100159 = round(gi[2], 2),
             pooled_g_FE = round(gfe, 2), pooled_g_FE_ci = sprintf("%.2f-%.2f", gfe - 1.96*sefe, gfe + 1.96*sefe),
             p_FE = signif(2 * pnorm(-abs(gfe) / sefe), 3),
             pooled_g_RE = round(gre, 2), p_RE = signif(2 * pnorm(-abs(gre) / sere), 3),
             I2 = round(I2, 2))
}))
write.csv(meta_tab, file.path(OUT, "metaDE_single_gene_validation.csv"), row.names = FALSE)
print(meta_tab)

# ---------- 4) ComBat 批次校正后重评估六基因模型（敏感性分析） ----------
transfer_auc <- function(tr, va, genes) {
  trm <- t(tr$expr[genes, , drop = FALSE]); y_tr <- as.numeric(tr$grp != "Control")
  m <- colMeans(trm); s <- apply(trm, 2, sd)
  trz <- sweep(sweep(trm, 2, m), 2, s, "/")
  f <- lrm(y_tr ~ ., data = as.data.frame(trz), x = TRUE, y = TRUE, penalty = 1)
  vam <- t(va$expr[genes, , drop = FALSE])
  vaz <- sweep(sweep(vam, 2, m), 2, s, "/")
  pv <- plogis(predict(f, as.data.frame(vaz)))
  as.numeric(auc(roc(as.numeric(va$grp != "Control"), pv, quiet = TRUE)))
}
combat_eval <- function(tr, va, lv, tag) {
  genes_all <- intersect(rownames(tr$expr), rownames(va$expr))
  mcomb <- cbind(tr$expr[genes_all, , drop = FALSE], va$expr[genes_all, , drop = FALSE])
  batch <- c(rep("train", ncol(tr$expr)), rep("valid", ncol(va$expr)))
  cb <- ComBat(dat = mcomb, batch = batch)
  trc <- list(expr = cb[, colnames(tr$expr), drop = FALSE], grp = tr$grp)
  vac <- list(expr = cb[, colnames(va$expr), drop = FALSE], grp = va$grp)
  a0 <- transfer_auc(tr, va, HUB6); a1 <- transfer_auc(trc, vac, HUB6)
  data.frame(pair = tag, genes_used = length(genes_all),
             transfer_AUC_original = round(a0, 3), transfer_AUC_combat = round(a1, 3))
}
cres <- rbind(combat_eval(d76293, d243066, "ARDS", "ARDS (GSE76293->GSE243066)"),
              combat_eval(d69528, d100159, "Sepsis", "Sepsis (GSE69528->GSE100159)"))
write.csv(cres, file.path(OUT, "combat_sensitivity_model_transfer.csv"), row.names = FALSE)
print(cres)

cat("\n=== Step 52 完成 ===\n")
