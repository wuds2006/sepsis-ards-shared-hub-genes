# ============================================================
# Step 59 —— 四 bulk 队列 PCA 补充图（审稿详细评语：批量效应评估）
# 每队列独立 PCA（top-2000 MAD 基因），样本按疾病/对照着色
# 目的：展示每队列内部疾病-对照分离度，无隐藏批次结构
# 产出：supp_panels/FigS16_cohort_PCA.pdf (+ PNG 300dpi + TIFF 1000dpi)
#       06_external_validation/reviewer_w5/cohort_PCA_coords.csv
# ============================================================
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
suppressPackageStartupMessages({ library(patchwork) })

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUTD <- file.path(BASE, "06_external_validation", "reviewer_w5")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

# ---------- loaders（镜像 step52） ----------
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
  rownames(e) <- df$gene
  list(expr = e, grp = grp[keep])
}

cohorts <- list(
  list(tag = "GSE76293", dis = "ARDS", platform = "Affymetrix HG-U133 Plus 2.0",
       expr = readRDS(file.path(BASE, "01_DEG", "GSE76293_gene_expr.rds")),
       grp = readRDS(file.path(BASE, "01_DEG", "GSE76293_grp.rds"))$group),
  list(tag = "GSE69528", dis = "Sepsis", platform = "Affymetrix HG-U133A",
       expr = readRDS(file.path(BASE, "01_DEG", "GSE69528_gene_expr.rds")),
       grp = readRDS(file.path(BASE, "01_DEG", "GSE69528_grp.rds"))$group),
  list(tag = "GSE243066", dis = "ARDS", platform = "Illumina NovaSeq 6000",
       expr = NULL, grp = NULL),
  list(tag = "GSE100159", dis = "Sepsis", platform = "Illumina HumanHT-12 v4",
       expr = NULL, grp = NULL)
)
v <- load_243066(); cohorts[[3]]$expr <- v$expr; cohorts[[3]]$grp <- v$grp
v <- load_100159(); cohorts[[4]]$expr <- v$expr; cohorts[[4]]$grp <- v$grp

pca_panel <- function(cfg, letter) {
  expr <- cfg$expr; grp <- factor(cfg$grp)
  mads <- apply(expr, 1, mad, na.rm = TRUE)
  sel <- order(mads, decreasing = TRUE)[seq_len(min(2000, nrow(expr)))]
  pc <- prcomp(t(expr[sel, ]), center = TRUE, scale. = FALSE)
  pv <- summary(pc)$importance[2, 1:2] * 100
  df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], group = grp)
  n_dis <- sum(grp != "Control"); n_ctl <- sum(grp == "Control")
  ggplot(df, aes(PC1, PC2, color = group)) +
    geom_point(size = 1.4, alpha = 0.85, stroke = 0.2) +
    scale_color_manual(values = c("Control" = COL_CONTROL,
                                  setNames(COL_DISEASE, setdiff(levels(grp), "Control")))) +
    labs(title = sprintf("%s  %s (%s)", letter, cfg$tag, cfg$dis),
         subtitle = sprintf("%s; %d %s / %d control", cfg$platform, n_dis, cfg$dis, n_ctl),
         x = sprintf("PC1 (%.1f%%)", pv[1]), y = sprintf("PC2 (%.1f%%)", pv[2])) +
    theme_ajrccm() +
    theme(legend.position = "bottom", legend.title = element_blank(),
          plot.subtitle = element_text(size = 5.5),
          aspect.ratio = 1)
}

plist <- lapply(seq_along(cohorts), function(i) pca_panel(cohorts[[i]], LETTERS[i]))
p <- (plist[[1]] + plist[[2]]) / (plist[[3]] + plist[[4]]) +
  plot_annotation(title = "Per-cohort principal-component analysis (batch-effect assessment)",
                  theme = theme(plot.title = element_text(size = 8, face = "bold")))

save_pdf(p, file.path(OUT_SUPP, "FigS16_cohort_PCA.pdf"), 178, 160)
print(p)

# 坐标 CSV（供补充材料/审稿人核查）
coords <- do.call(rbind, lapply(seq_along(cohorts), function(i) {
  cfg <- cohorts[[i]]; expr <- cfg$expr
  mads <- apply(expr, 1, mad, na.rm = TRUE)
  sel <- order(mads, decreasing = TRUE)[seq_len(min(2000, nrow(expr)))]
  pc <- prcomp(t(expr[sel, ]), center = TRUE, scale. = FALSE)
  data.frame(cohort = cfg$tag, disease = cfg$dis, group = cfg$grp,
             PC1 = pc$x[, 1], PC2 = pc$x[, 2])
}))
write.csv(coords, file.path(OUTD, "cohort_PCA_coords.csv"), row.names = FALSE)
cat("\n=== Step 59 完成 ===\n")
