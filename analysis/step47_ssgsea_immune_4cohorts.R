# ============================================================
# Step 47 —— ssGSEA 免疫浸润分析（四队列）
# 17 免疫细胞 marker 集（沿用 step2c，canonical marker sets）
# 每队列: ssGSEA 评分 → Wilcoxon 组间比较（BH 校正）→ 箱线图面板
#         6 hub 基因 × 17 细胞 Spearman 相关 → 显著性热图面板
# 输出: 07_hub_immune/ssgsea_immune/*.csv
#       AJRCCM_submission/03_figures/Figures_ggplot2_NPG/supp_panels/FigS15_*.pdf
# ============================================================
suppressPackageStartupMessages({
  library(GSVA)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
source("scripts/common.R")

BASE  <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
SUB   <- file.path(BASE, "AJRCCM_submission")
OUT   <- file.path(BASE, "07_hub_immune", "ssgsea_immune")
PANEL <- file.path(SUB, "03_figures", "Figures_ggplot2_NPG", "supp_panels")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

HUB6 <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")

# ---- 17 免疫细胞 canonical marker 集（与 step2c 完全一致） ----
immune_sets <- list(
  CD8_T = c("CD8A","CD8B","GZMA","GZMB"),
  CD4_T = c("CD4","IL7R","CD3D","CD3E"),
  Treg = c("FOXP3","IL2RA","CTLA4","TIGIT"),
  Th1 = c("TBX21","IFNG","STAT4","IL12RB2"),
  Th2 = c("GATA3","IL4","IL5","IL13"),
  Th17 = c("RORC","IL17A","IL22","IL23R"),
  NK = c("NKG7","KLRD1","NCAM1","GNLY"),
  B_cell = c("MS4A1","CD79A","CD19","CD79B"),
  Plasma = c("MZB1","XBP1","TNFRSF17","SDC1"),
  Monocyte = c("CD14","LYZ","FCN1","S100A12"),
  M1_Macro = c("IRF5","NOS2","TNF","IL1B"),
  M2_Macro = c("CD163","MRC1","ARG1","MSR1"),
  Neutrophil = c("FCGR3B","CSF3R","S100A8","CEACAM8"),
  DC = c("ITGAX","CD1C","FCER1A","CST3"),
  pDC = c("LILRA4","CLEC4C","IRF7","TCF4"),
  Mast = c("KIT","TPSAB1","CPA3","MS4A2"),
  Eosinophil = c("SIGLEC8","CCR3","RNASE2","IL5RA"),
  MDSC = c("CD33","FUT4","ITGAM","S100A9")
)
CELL17 <- names(immune_sets)

## ---------------- dataset loaders (ported from step41) ----------------
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
  cpm <- t(t(mat) / colSums(mat) * 1e6)
  logcpm <- log2(cpm + 1)
  keep <- g %in% c("ARDS", "Control")
  list(expr = logcpm[, keep, drop = FALSE], grp = data.frame(group = g[keep]))
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

d1 <- load_disc("GSE76293")
d2 <- load_disc("GSE69528")
d3 <- load_GSE243066()
d4 <- load_GSE100159()

COHORTS <- list(
  list(id = "GSE76293",  dis = "ARDS",   label = "GSE76293 (ARDS, PMN)",
       expr = d1$expr, grp = d1$grp$group),
  list(id = "GSE69528",  dis = "Sepsis", label = "GSE69528 (Sepsis, whole blood)",
       expr = d2$expr, grp = d2$grp$group),
  list(id = "GSE243066", dis = "ARDS",   label = "GSE243066 (ARDS, PBMC)",
       expr = d3$expr, grp = d3$grp$group),
  list(id = "GSE100159", dis = "Sepsis", label = "GSE100159 (Sepsis, whole blood)",
       expr = d4$expr, grp = d4$grp$group)
)

sig_stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 1e-4, "****", ifelse(p < 1e-3, "***",
         ifelse(p < 1e-2, "**", ifelse(p < 0.05, "*", "ns")))))
}

pdf_dev <- function(file, w_in, h_in) {
  grDevices::pdf(file, width = w_in, height = h_in, family = "Helvetica",
                 useDingbats = FALSE)
}

MM2IN <- 1 / 25.4
summary_rows <- list()

for (z in COHORTS) {
  expr <- as.matrix(z$expr)
  grp  <- z$grp; names(grp) <- colnames(expr)
  cat(sprintf("\n===== %s (%d samples: %s vs %s) =====\n", z$id, ncol(expr),
              z$dis, "Control"))

  ## ---- ssGSEA ----
  sets <- lapply(immune_sets, function(g) intersect(g, rownames(expr)))
  sets <- sets[lengths(sets) >= 2]
  cat("cell types scored:", length(sets), "of", length(CELL17), "\n")
  sc <- gsva(ssgseaParam(exprData = expr, geneSets = sets, minSize = 2))
  sc_df <- data.frame(sample = colnames(sc), group = grp[colnames(sc)],
                      t(as.matrix(sc)), check.names = FALSE)
  write.csv(sc_df, file.path(OUT, paste0(z$id, "_ssgsea_scores.csv")), row.names = FALSE)

  ## ---- Wilcoxon 组间比较 ----
  cells <- rownames(sc)
  wil <- do.call(rbind, lapply(cells, function(cell) {
    x <- sc[cell, grp == z$dis]; y <- sc[cell, grp == "Control"]
    p <- tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
    data.frame(cell = cell, median_disease = median(x), median_control = median(y),
               mean_disease = mean(x), mean_control = mean(y),
               delta = mean(x) - mean(y), P = p)
  }))
  wil$padj <- p.adjust(wil$P, method = "BH")
  wil <- wil[order(wil$P), ]
  wil$cell <- factor(wil$cell, levels = CELL17)
  wil <- wil[order(wil$cell), ]
  write.csv(wil, file.path(OUT, paste0(z$id, "_ssgsea_wilcox.csv")), row.names = FALSE)
  cat("significant (P<0.05):", sum(wil$P < 0.05, na.rm = TRUE), "of", nrow(wil), "\n")
  print(head(wil[, c("cell", "mean_disease", "mean_control", "P", "padj")], 6))

  ## ---- 箱线图面板（17 细胞 facet，星号=名义 P） ----
  long <- as.data.frame(t(sc)) |> tibble::rownames_to_column("sample")
  long$group <- factor(grp[long$sample], levels = c("Control", z$dis))
  long <- pivot_longer(long, cols = all_of(cells), names_to = "cell", values_to = "score")
  long$cell <- factor(long$cell, levels = CELL17)
  st <- wil[, c("cell", "P")]; st$star <- sig_stars(st$P)
  long$star <- st$star[match(long$cell, st$cell)]
  ymax <- tapply(long$score, long$cell, max, na.rm = TRUE)
  yrng <- tapply(long$score, long$cell, function(v) diff(range(v, na.rm = TRUE)))
  anno <- data.frame(cell = names(ymax),
                     y = ymax + 0.12 * ifelse(yrng > 0, yrng, abs(ymax) + 0.05),
                     star = st$star[match(names(ymax), st$cell)])
  anno <- anno[!is.na(anno$star) & !is.na(anno$y) & anno$star != "ns" & anno$star != "", ]
  anno$cell <- factor(anno$cell, levels = CELL17)
  anno$group <- factor(z$dis, levels = c("Control", z$dis))

  p_box <- ggplot(long, aes(x = group, y = score, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.62, linewidth = 0.3) +
    geom_jitter(width = 0.16, size = 0.4, alpha = 0.55) +
    geom_text(data = anno, aes(x = group, y = y, label = star),
              inherit.aes = FALSE, size = 2.6, hjust = 0.5) +
    facet_wrap(~cell, scales = "free_y", ncol = 6) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.22))) +
    scale_fill_manual(values = setNames(c("#8FA2D6", "#D99AAE"), c("Control", z$dis))) +
    labs(title = z$label, y = "ssGSEA score", x = NULL) +
    theme_pub(base_size = 7) +
    theme(legend.position = "none",
          plot.title = element_text(size = 8, face = "bold", hjust = 0),
          strip.text = element_text(size = 6.2, face = "bold"),
          axis.text.x = element_text(size = 5.6, angle = 30, hjust = 1),
          axis.text.y = element_text(size = 5.6))
  bw <- 162 * MM2IN; bh <- 52 * MM2IN
  pdf_dev(file.path(PANEL, paste0("FigS15_box_", z$id, ".pdf")), bw, bh)
  print(p_box); dev.off()
  ggsave(file.path(PANEL, paste0("FigS15_box_", z$id, ".png")), p_box,
         width = bw, height = bh, dpi = 300)

  ## ---- Spearman 相关: 6 hub × 17 细胞 ----
  genes <- intersect(HUB6, rownames(expr))
  cor_l <- do.call(rbind, lapply(genes, function(g) {
    do.call(rbind, lapply(cells, function(cell) {
      ct <- suppressWarnings(cor.test(expr[g, ], sc[cell, ], method = "spearman",
                                      exact = FALSE))
      data.frame(gene = g, cell = cell, rho = unname(ct$estimate), P = ct$p.value)
    }))
  }))
  cor_l$padj <- p.adjust(cor_l$P, method = "BH")
  write.csv(cor_l, file.path(OUT, paste0(z$id, "_ssgsea_hub_spearman.csv")),
            row.names = FALSE)
  cor_l$gene <- factor(cor_l$gene, levels = rev(HUB6))
  cor_l$cell <- factor(cor_l$cell, levels = CELL17)
  cor_l$star <- sig_stars(cor_l$P)

  p_heat <- ggplot(cor_l, aes(x = cell, y = gene, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = ifelse(star != "ns", star, "")), size = 2.4, vjust = 0.7) +
    scale_fill_gradient2(low = "#3C5488FF", mid = "white", high = "#E64B35FF",
                         midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho") +
    labs(title = z$label, x = NULL, y = NULL) +
    theme_pub(base_size = 7) +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0),
          axis.text.x = element_text(size = 5.6, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 6.4, face = "italic"),
          legend.position = "right", legend.key.height = grid::unit(7, "mm"),
          legend.title = element_text(size = 6), legend.text = element_text(size = 5.6),
          panel.grid = element_blank())
  hw <- 79 * MM2IN; hh <- 50 * MM2IN
  pdf_dev(file.path(PANEL, paste0("FigS15_cor_", z$id, ".pdf")), hw, hh)
  print(p_heat); dev.off()
  ggsave(file.path(PANEL, paste0("FigS15_cor_", z$id, ".png")), p_heat,
         width = hw, height = hh, dpi = 300)

  ## ---- 汇总行（供正文撰写） ----
  w2 <- wil[!is.na(wil$P) & wil$P < 0.05, ]
  up <- as.character(w2$cell[w2$delta > 0]); dn <- as.character(w2$cell[w2$delta < 0])
  summary_rows[[z$id]] <- data.frame(
    cohort = z$id,
    n_up = length(up), up = paste(up, collapse = ";"),
    n_down = length(dn), down = paste(dn, collapse = ";"))
}
summ <- do.call(rbind, summary_rows)
write.csv(summ, file.path(OUT, "ssgsea_direction_summary.csv"), row.names = FALSE)
cat("\n===== Step47 完成 =====\n")
print(summ[, c("cohort", "n_up", "n_down")])
