#!/usr/bin/env Rscript
# step42b: Figure S4 v2 panels - hub-gene x CIBERSORT immune-cell correlation heatmaps
# Replaces the former ssGSEA-based S4 (and absorbs old Fig5C/D content, now
# recomputed against CIBERSORT fractions for methodological coherence).
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

BASE <- PROJ
CBD <- file.path(BASE, "07_hub_immune", "cibersort")

CELL22 <- c("B cells naive","B cells memory","Plasma cells","T cells CD8",
  "T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated",
  "T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta",
  "NK cells resting","NK cells activated","Monocytes","Macrophages M0",
  "Macrophages M1","Macrophages M2","Dendritic cells resting","Dendritic cells activated",
  "Mast cells resting","Mast cells activated","Eosinophils","Neutrophils")
CELL_SHORT <- c("B naive","B memory","Plasma","CD8 T","CD4 T naive","CD4 T mem rest",
  "CD4 T mem act","Tfh","Tregs","Tgd","NK rest","NK act","Monocytes","M0","M1","M2",
  "DC rest","DC act","Mast rest","Mast act","Eosinophils","Neutrophils")
names(CELL_SHORT) <- CELL22

load_disc <- function(ds) readRDS(file.path(BASE, "01_DEG", paste0(ds, "_gene_expr.rds")))

load_GSE243066_logcpm <- function() {
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
  g <- ifelse(grepl("_A[0-9]", colnames(mat)), "ARDS",
       ifelse(grepl("_H[0-9]", colnames(mat)), "Control", "Exclude"))
  logcpm <- log2(t(t(mat) / colSums(mat) * 1e6) + 1)
  keep <- g %in% c("ARDS", "Control")
  logcpm[, keep, drop = FALSE]
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
  common <- intersect(rownames(e4), ann$probe); e4 <- e4[common, , drop = FALSE]
  gmap <- setNames(ann$gene[match(common, ann$probe)], common)
  dfp <- data.frame(probe = common, gene = gmap, mean = rowMeans(e4), stringsAsFactors = FALSE)
  dfp <- dfp[order(dfp$mean, decreasing = TRUE), ]; dfp <- dfp[!duplicated(dfp$gene), ]
  e4 <- e4[dfp$probe, , drop = FALSE]; rownames(e4) <- dfp$gene
  e4
}

expr_list <- list(GSE76293 = load_disc("GSE76293"),
                  GSE69528 = load_disc("GSE69528"),
                  GSE243066 = load_GSE243066_logcpm(),
                  GSE100159 = load_GSE100159())
PANEL <- c(GSE76293 = "A", GSE69528 = "B", GSE243066 = "C", GSE100159 = "D")
TAG <- c(GSE76293 = "ARDS discovery (PMN fractions)", GSE69528 = "Sepsis discovery (whole blood)",
         GSE243066 = "ARDS validation", GSE100159 = "Sepsis validation")

for (ds in names(expr_list)) {
  expr <- as.matrix(expr_list[[ds]])
  res <- read.csv(file.path(CBD, paste0(ds, "_CIBERSORT_results.csv")), check.names = FALSE)
  samp <- intersect(colnames(expr), res$sample)
  expr <- expr[, samp, drop = FALSE]
  rr <- res[match(samp, res$sample), , drop = FALSE]
  genes <- intersect(HUB6, rownames(expr))
  cor_m <- matrix(NA_real_, nrow = length(genes), ncol = length(CELL22),
                  dimnames = list(genes, CELL22))
  for (g in genes) for (ct in CELL22)
    cor_m[g, ct] <- suppressWarnings(cor(as.numeric(expr[g, ]), rr[[ct]], method = "spearman"))
  write.csv(cor_m, file.path(CBD, paste0(ds, "_hub_ciber_cor.csv")))

  long <- as.data.frame(as.table(cor_m))
  names(long) <- c("gene", "cell", "rho")
  long$gene <- factor(long$gene, levels = rev(genes))
  long$cell <- factor(long$cell, levels = CELL22)
  p <- ggplot(long, aes(cell, gene, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#3C5488FF", mid = "white", high = "#E64B35FF",
                         midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho",
                         na.value = "grey88") +
    scale_x_discrete(breaks = CELL22, labels = CELL_SHORT[CELL22]) +
    scale_y_discrete(labels = function(x) parse(text = gene_italic(x))) +
    labs(title = paste0(ds, " (", TAG[[ds]], ")"), x = NULL, y = NULL) +
    theme_ajrccm(base_size = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
          legend.position = "right", legend.key.height = unit(9, "mm"),
          legend.title = element_text(size = 5.5), legend.text = element_text(size = 4.8))
  save_pdf(p, file.path(OUT_SUPP, paste0("FigS4", PANEL[[ds]], "_hub_ciber_cor_", ds, ".pdf")), 88, 50)
}
cat("step42b: Figure S4 v2 panels done.\n")
