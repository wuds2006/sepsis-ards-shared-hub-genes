#!/usr/bin/env Rscript
# Figure 2 panels (ggplot2 + ggsci NPG):
#   A-D: six hub-gene expression boxplots in 4 cohorts
#   E:   ROC curves (6 genes + combined signature) per cohort
#   F:   11-classifier CV-AUC bars + SHAP mean|SHAP| bars (both discovery cohorts)
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
suppressPackageStartupMessages({ library(dplyr); library(pROC) })

BASE <- PROJ

## ---------- data loaders (ported from step2d_expression_ROC.R) ----------
load_gene_expr <- function(ds) {
  expr <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_gene_expr.rds")))
  grp  <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_grp.rds")))
  list(expr = expr, grp = grp)
}
d1 <- load_gene_expr("GSE76293")   # ARDS discovery
d2 <- load_gene_expr("GSE69528")   # Sepsis discovery

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
g3 <- ifelse(grepl("_A[0-9]", sn), "ARDS", ifelse(grepl("_H[0-9]", sn), "Control", "Exclude"))
cpm <- t(t(mat) / colSums(mat) * 1e6); logcpm <- log2(cpm + 1)
keep <- g3 %in% c("ARDS", "Control"); logcpm <- logcpm[, keep, drop = FALSE]; g3 <- g3[keep]
d3 <- list(expr = logcpm, grp = data.frame(group = g3))

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
d4 <- list(expr = e4, grp = data.frame(group = g4))

datasets <- list(
  A = list(id = "GSE76293",  label = "GSE76293 ARDS discovery (n = 12 vs 12)",   disease = "ARDS",   d = d1),
  B = list(id = "GSE243066", label = "GSE243066 ARDS validation (n = 34 vs 15)", disease = "ARDS",   d = d3),
  C = list(id = "GSE69528",  label = "GSE69528 sepsis discovery (n = 83 vs 28)", disease = "Sepsis", d = d2),
  D = list(id = "GSE100159", label = "GSE100159 sepsis validation (n = 33 vs 12)", disease = "Sepsis", d = d4)
)

## ---------- A-D: boxplots ----------
mk_box <- function(z, panel_tag) {
  g <- intersect(HUB6, rownames(z$d$expr))
  yy <- z$d$grp$group
  dis <- z$disease
  long <- do.call(rbind, lapply(g, function(gene)
    data.frame(gene = gene, expr = z$d$expr[gene, ], group = yy)))
  long$gene <- factor(long$gene, levels = HUB6)
  long$group <- factor(long$group, levels = c("Control", dis))
  pval <- sapply(g, function(gene)
    wilcox.test(z$d$expr[gene, yy == "Control"], z$d$expr[gene, yy == dis],
                exact = FALSE, correct = FALSE)$p.value)
  star <- setNames(sig_stars(pval), g)
  ann <- data.frame(gene = factor(g, levels = HUB6),
                    y = sapply(g, function(gene) max(z$d$expr[gene, ]) * 1.04),
                    lab = star[g])
  ggplot(long, aes(gene, expr, fill = group)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.3, alpha = 0.9, width = 0.65) +
    geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.65),
               size = 0.35, alpha = 0.5, color = "grey25", stroke = 0) +
    geom_text(data = ann, aes(x = gene, y = y, label = lab), inherit.aes = FALSE, size = 2.6) +
    scale_fill_manual(values = setNames(c(COL_CONTROL, COL_DISEASE), c("Control", dis))) +
    scale_x_discrete(labels = function(x) parse(text = gene_italic(x))) +
    labs(title = paste0(panel_tag, "  ", z$label), x = NULL, y = "Expression") +
    theme_ajrccm() + theme(legend.position = "top", legend.title = element_blank())
}
for (tag in names(datasets)) {
  p <- mk_box(datasets[[tag]], tag)
  save_pdf(p, file.path(OUT_MAIN, paste0("Fig2", tag, "_boxplot_", datasets[[tag]]$id, ".pdf")), 82, 62)
}

## ---------- E: ROC curves per cohort ----------
roc_one <- function(z) {
  g <- intersect(HUB6, rownames(z$d$expr))
  yy <- as.numeric(z$d$grp$group != "Control")
  rr <- do.call(rbind, lapply(g, function(gene) {
    ro <- roc(yy, z$d$expr[gene, ], quiet = TRUE)
    data.frame(gene = gene, FPR = 1 - ro$specificities, TPR = ro$sensitivities,
               AUC = as.numeric(ro$auc))
  }))
  X <- t(z$d$expr[g, , drop = FALSE])
  m <- glm(yy ~ X, family = binomial)
  sc <- predict(m, type = "response")
  ro <- roc(yy, sc, quiet = TRUE)
  rbind(rr, data.frame(gene = "Combined", FPR = 1 - ro$specificities,
                       TPR = ro$sensitivities, AUC = as.numeric(ro$auc)))
}
mk_roc <- function(z, tag) {
  rl <- roc_one(z)
  rl$gene <- factor(rl$gene, levels = c(HUB6, "Combined"))
  auc_lab <- rl %>% distinct(gene, AUC) %>%
    mutate(lab = ifelse(gene == "Combined",
                        paste0("'Combined ('*", sprintf("%.3f", AUC), "*')'"),
                        paste0("italic('", gene, "')*' (", sprintf("%.3f", AUC), ")'")))
  cols <- c(setNames(pal_npg(6), HUB6), Combined = "grey20")
  lab_vec <- setNames(auc_lab$lab, auc_lab$gene)
  ggplot(rl, aes(FPR, TPR, color = gene)) +
    geom_line(linewidth = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    scale_color_manual(values = cols, labels = function(x) parse(text = lab_vec[x])) +
    coord_equal() +
    labs(title = paste0(tag, "  ", z$label), x = "1 - Specificity", y = "Sensitivity", color = NULL) +
    theme_ajrccm() +
    theme(legend.position = "right",
          legend.text = element_text(size = 5.5))
}
for (tag in names(datasets)) {
  p <- mk_roc(datasets[[tag]], paste0("E", match(tag, names(datasets))))
  save_pdf(p, file.path(OUT_MAIN, paste0("Fig2E", match(tag, names(datasets)), "_ROC_", datasets[[tag]]$id, ".pdf")), 88, 68)
}

## ---------- F: ML CV-AUC bars + SHAP bars ----------
auc_ards   <- read.csv(file.path(BASE, "05_ML_SHAP/ML_AUC_results.csv"))         %>% mutate(cohort = "ARDS (GSE76293)")
auc_sepsis <- read.csv(file.path(BASE, "05_ML_SHAP/sepsis/ML_AUC_results_sepsis.csv")) %>% mutate(cohort = "Sepsis (GSE69528)")
auc_all <- rbind(auc_ards, auc_sepsis)
pF1 <- ggplot(auc_all, aes(reorder(model, mean_AUC), mean_AUC, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62, color = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = pmax(mean_AUC - sd_AUC, 0), ymax = pmin(mean_AUC + sd_AUC, 1)),
                position = position_dodge(width = 0.7), width = 0.25, linewidth = 0.25) +
  scale_fill_npg() + coord_cartesian(ylim = c(0.5, 1.02)) +
  labs(title = "F1  Cross-validated AUC of 11 classifiers",
       subtitle = paste("Error bars: SD of cross-validated AUC.",
                        "Apparent AUC = 1.000 corresponds to an exact",
                        "one-sided 97.5% lower bound of 0.975",
                        "(ARDS training, 12 vs 12); label-permutation",
                        "P \u2264 0.014 for all classifiers",
                        "(Supplementary Tables S13, S14)",
                        sep = "\n"),
       x = NULL, y = "Cross-validated AUC", fill = NULL) +
  theme_ajrccm() + theme(legend.position = "top",
                         plot.subtitle = element_text(size = 5.2, lineheight = 0.9))
save_pdf(pF1, file.path(OUT_MAIN, "Fig2F1_ML_AUC_barplot.pdf"), 82, 62)

sh_ards   <- read.csv(file.path(BASE, "05_ML_SHAP/SHAP_importance.csv"))         %>% mutate(cohort = "ARDS (GSE76293)")
sh_sepsis <- read.csv(file.path(BASE, "05_ML_SHAP/sepsis/SHAP_importance_sepsis.csv")) %>% mutate(cohort = "Sepsis (GSE69528)")
sh_all <- rbind(sh_ards, sh_sepsis)
top_g <- sh_all %>% group_by(gene) %>% summarise(m = max(mean_abs_shap), .groups = "drop") %>%
  arrange(desc(m)) %>% slice_head(n = 10) %>% pull(gene)
sh_top <- sh_all %>% filter(gene %in% top_g) %>%
  mutate(gene = factor(gene, levels = rev(top_g)))
pF2 <- ggplot(sh_top, aes(gene, mean_abs_shap, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62, color = "black", linewidth = 0.2) +
  scale_fill_npg() + coord_flip() +
  scale_x_discrete(labels = function(x) parse(text = gene_italic(x))) +
  labs(title = "F2  SHAP feature importance (top 10)", x = NULL, y = "Mean |SHAP value|", fill = NULL) +
  theme_ajrccm() + theme(legend.position = "top")
save_pdf(pF2, file.path(OUT_MAIN, "Fig2F2_SHAP_barplot.pdf"), 82, 68)

cat("Figure 2 panels done.\n")
