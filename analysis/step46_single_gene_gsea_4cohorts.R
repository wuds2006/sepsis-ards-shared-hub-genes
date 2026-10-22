#!/usr/bin/env Rscript
# step46: single-gene ssGSEA-style annotation (fgsea, MSigDB Hallmark v2024.1) for
# the six hub genes, extended from the two discovery cohorts to ALL FOUR bulk
# cohorts (GSE76293, GSE69528, GSE243066, GSE100159) per the revised Methods.
# Rank metric: Pearson correlation of every other gene with the hub gene.
# Outputs (12_single_gene_gsea/):
#   single_gene_gsea_all_results_4cohorts.csv / single_gene_gsea_top2_4cohorts.csv
# Outputs (supp_panels + composites): rebuilt Figure S8 (6 genes x 4 cohorts grid)
suppressMessages({
  library(fgsea); library(ggplot2); library(data.table); library(patchwork)
})
set.seed(42)
root <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
outdir <- file.path(root, "12_single_gene_gsea")
panel_dir <- file.path(root, "AJRCCM_submission/03_figures/Figures_ggplot2_NPG/supp_panels")
dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)

hub <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")
cohort_lab <- c(GSE76293 = "GSE76293 (ARDS discovery)",
                GSE69528 = "GSE69528 (Sepsis discovery)",
                GSE243066 = "GSE243066 (ARDS validation)",
                GSE100159 = "GSE100159 (Sepsis validation)")

pathways <- gmtPathways(file.path(root, "00_rawdata/hallmark/h.all.v2024.1.Hs.symbols.gmt"))
pal <- c("#6A5C95", "#D99AAE", "#E7A889", "#A9CBB8", "#CFE3E6", "#8FA2D6", "#5B608C")

# ---- loaders (validation cohorts ported from step41; log scale for ranking) ----
load_GSE243066 <- function() {
  dd <- file.path(root, "00_rawdata/validation/GSE243066")
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
  logcpm[, keep, drop = FALSE]
}
load_GSE100159 <- function() {
  DD <- file.path(root, "00_rawdata/validation/GSE100159")
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
  keep <- g4 %in% c("Sepsis", "Control"); e4 <- e4[, keep, drop = FALSE]
  common <- intersect(rownames(e4), ann$probe); e4 <- e4[common, , drop = FALSE]
  gmap <- setNames(ann$gene[match(common, ann$probe)], common)
  gmean <- rowMeans(e4)
  dfp <- data.frame(probe = common, gene = gmap, mean = gmean, stringsAsFactors = FALSE)
  dfp <- dfp[order(dfp$mean, decreasing = TRUE), ]; dfp <- dfp[!duplicated(dfp$gene), ]
  e4 <- e4[dfp$probe, , drop = FALSE]; rownames(e4) <- dfp$gene
  e4   # keep linear scale: raw matrix contains negative intensities (CD177/UPP1),
       # log2(x+1) would produce NaN; all prior GSE100159 analyses use this matrix
}

exprs <- list(
  GSE76293  = as.matrix(readRDS(file.path(root, "01_DEG/GSE76293_gene_expr.rds"))),
  GSE69528  = as.matrix(readRDS(file.path(root, "01_DEG/GSE69528_gene_expr.rds"))),
  GSE243066 = load_GSE243066(),
  GSE100159 = load_GSE100159()
)

run_es <- function(stats, geneset) {
  inset <- names(stats) %in% geneset
  N <- length(stats); Nh <- sum(inset)
  if (Nh == 0) return(rep(0, N))
  hit <- ifelse(inset, abs(stats) / sum(abs(stats)[inset]), 0)
  miss <- ifelse(inset, 0, -1 / (N - Nh))
  cumsum(hit + miss)
}

all_res <- list(); plots <- list()
pi <- 0
# discovery cohorts: reuse ORIGINAL results (archived run) so manuscript numbers
# (Results text NES values) stay exact; validation cohorts: fresh fgsea run
orig_all <- fread(file.path(outdir, "single_gene_gsea_all_results.csv"))
orig_top <- fread(file.path(outdir, "single_gene_gsea_top2_pathways.csv"))
DISC <- c("GSE76293", "GSE69528")
for (ds in names(cohort_lab)) {
  m <- exprs[[ds]]
  stopifnot(all(hub %in% rownames(m)))
  for (g in hub) {
    stats <- suppressWarnings(cor(t(m), m[g, ], method = "pearson"))
    stats <- as.numeric(stats); names(stats) <- rownames(m)
    stats <- stats[!is.na(stats)]
    stats <- sort(stats[names(stats) != g], decreasing = TRUE)

    if (ds %in% DISC) {
      fg <- as.data.frame(orig_all[cohort == ds & hub_gene == g])
      top2 <- as.data.frame(orig_top[cohort == ds & hub_gene == g])
    } else {
      fg <- as.data.frame(fgsea(pathways = pathways, stats = stats, minSize = 10, maxSize = 500))
      fg$cohort <- ds; fg$hub_gene <- g; fg$rank_metric <- "Pearson r"
      top2 <- fg[order(-abs(fg$NES)), ][1:2, ]
    }
    all_res[[paste(ds, g)]] <- fg

    cat(sprintf("%s %s | top2: %s (NES=%.2f, p=%.2g) | %s (NES=%.2f, p=%.2g)\n",
                ds, g, top2$pathway[1], top2$NES[1], top2$pval[1],
                top2$pathway[2], top2$NES[2], top2$pval[2]))

    pi <- pi + 1
    es_df <- do.call(rbind, lapply(1:2, function(k) {
      es <- run_es(stats, pathways[[top2$pathway[k]]])
      data.frame(rank = seq_along(es), ES = es,
                 pathway = top2$pathway[k], NES = top2$NES[k], col = pal[k])
    }))
    met_df <- data.frame(rank = seq_along(stats), r = stats,
                         hit = names(stats) %in% unlist(pathways[top2$pathway]))
    met_df$col <- ifelse(met_df$hit, "#C0392B", "grey70")

    neg <- top2$NES[1] < 0
    txt_y <- if (neg) max(es_df$ES) * 0.96 else min(es_df$ES) * 0.96
    txt_x <- max(es_df$rank) * 0.98
    ttl <- paste0(LETTERS[pi], "  ", g, "\n", ds)
    p_es <- ggplot(es_df, aes(rank, ES, colour = col)) +
      geom_line(linewidth = 0.4) +
      scale_colour_identity() +
      geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.25) +
      ggtitle(ttl) +
      annotate("text", x = txt_x, y = txt_y, hjust = 1,
               vjust = if (neg) 1 else 0, size = 1.1, lineheight = 0.92,
               label = paste0(gsub("_", " ", gsub("HALLMARK_", "", top2$pathway)),
                              " (", sprintf("%.2f", top2$NES), ")",
                              collapse = "\n")) +
      labs(x = NULL, y = "Running ES") +
      theme_classic(base_size = 6) +
      theme(legend.position = "none", plot.margin = margin(1.5, 1.5, 0, 1.5),
            plot.title = element_text(size = 5.2, face = "bold", lineheight = 0.95,
                                      hjust = 0),
            axis.text = element_text(size = 4.6), axis.title = element_text(size = 5.2))

    p_met <- ggplot(met_df, aes(rank, r)) +
      geom_area(aes(fill = col), colour = NA, alpha = 0.9) +
      scale_fill_identity() +
      labs(x = "Gene rank", y = "Pearson r") +
      theme_classic(base_size = 6) +
      theme(legend.position = "none", plot.margin = margin(0, 1.5, 1.5, 1.5),
            axis.text = element_text(size = 4.6), axis.title = element_text(size = 5.2))

    plots[[pi]] <- p_es / p_met + plot_layout(heights = c(2, 1))
  }
}

res <- rbindlist(all_res, fill = TRUE)
res <- res[order(cohort, hub_gene, -abs(NES))]
fwrite(res, file.path(outdir, "single_gene_gsea_all_results_4cohorts.csv"))
top_tbl <- res[, .SD[1:2], by = .(cohort, hub_gene)]
fwrite(top_tbl, file.path(outdir, "single_gene_gsea_top2_4cohorts.csv"))

# Figure S8 (rebuilt): 6 genes (cols) x 4 cohorts (rows), 178 x 205 mm
fig <- wrap_plots(plots, ncol = 6) +
  plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2)))
grDevices::pdf(file.path(panel_dir, "FigureS8_4cohorts.pdf"),
               width = 178 / 25.4, height = 205 / 25.4,
               family = "Helvetica", useDingbats = FALSE)
print(fig); invisible(dev.off())
cat("step46 done. 24 panels ->", file.path(panel_dir, "FigureS8_4cohorts.pdf"), "\n")
