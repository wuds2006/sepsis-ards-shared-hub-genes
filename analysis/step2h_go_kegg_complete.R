# ============================================================
# Step 2h: GO 富集(BP+MF+CC) 桑葚气泡图 + KEGG 弦图（本地数据）
# 两组: 179 shared CGs + 6 hub genes
# ============================================================
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db)
  library(circlize); library(ggplot2); library(dplyr)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "AJRCCM_submission", "07_figure_exports", "Figure3_GO_sankey_chord")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

hub <- readLines(file.path(BASE, "07_hub_immune", "hub_genes.txt"))
cgs <- readLines(file.path(BASE, "03_intersection", "Shared_CGs.txt"))

# ---- GO(BP/MF/CC) 富集（本地 org.Hs.eg.db）----
go_enrich <- function(genes) {
  eg <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  list(BP = enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05, readable = TRUE),
       MF = enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05, readable = TRUE),
       CC = enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05, readable = TRUE))
}

# ---- 桑葚气泡图：GO BP+MF+CC 合并 ----
sankey_go <- function(ego, tag, topEach = 8) {
  dfs <- list(BP = as.data.frame(ego$BP), MF = as.data.frame(ego$MF), CC = as.data.frame(ego$CC))
  dfs <- Filter(function(d) nrow(d) > 0, dfs)   # 去掉空 ontology(如 6hub 的 MF=0)
  tops <- lapply(names(dfs), function(o){
    d <- dfs[[o]][order(dfs[[o]]$pvalue), ]
    d <- d[1:min(topEach, nrow(d)), ]; d$ont <- o; d
  })
  go <- do.call(rbind, tops)

  long <- do.call(rbind, lapply(1:nrow(go), function(i){
    gs <- sort(strsplit(go$geneID[i], "/")[[1]])
    gs <- gs[1:min(6, length(gs))]
    data.frame(gene = gs, term = paste0("[", go$ont[i], "] ", go$Description[i]),
               ont = go$ont[i], count = go$Count[i], neglogp = -log10(go$pvalue[i]),
               stringsAsFactors = FALSE)
  }))
  terms <- unique(long$term); genes <- unique(long$gene)
  long$term_f <- factor(long$term, levels = rev(terms))
  long$gene_f <- factor(long$gene, levels = rev(genes))
  long$y_term <- as.numeric(long$term_f); long$y_gene <- as.numeric(long$gene_f)
  ti <- long[!duplicated(long$term), c("term_f","y_term","count","neglogp","ont")]
  gi <- long[!duplicated(long$gene), c("gene_f","y_gene")]
  ont_col <- c(BP = "#6A5C95", MF = "#A9CBB8", CC = "#8FA2D6")

  p <- ggplot() +
    geom_curve(data = long, aes(x = 1, y = y_gene, xend = 1.8, yend = y_term, color = ont),
               curvature = -0.18, alpha = 0.32, linewidth = 0.4) +
    geom_point(data = gi, aes(x = 1, y = y_gene), size = 1.4, color = "grey40") +
    geom_text(data = gi, aes(x = 0.97, y = y_gene, label = gene_f),
              size = 2.5, hjust = 1, check_overlap = TRUE) +
    geom_point(data = ti, aes(x = 1.8, y = y_term, size = count, fill = neglogp),
               shape = 21, color = "grey30", alpha = 0.92) +
    geom_text(data = ti, aes(x = 1.83, y = y_term, label = term_f, color = ont),
              size = 2.8, hjust = 0, fontface = "bold") +
    scale_size_continuous(range = c(3, 15), name = "Gene count") +
    scale_fill_gradient(low = "#CFE3E6", high = "#6A5C95", name = "-log10(P)") +
    scale_color_manual(values = ont_col, guide = "none") +
    scale_x_continuous(limits = c(0.15, 4.6)) +
    theme_void(base_size = 9) +
    theme(legend.position = "right", plot.title = element_text(size = 11, face = "bold", hjust = 0.5))
  ggsave(file.path(OUT, paste0(tag, "_GO_sankey_bubble.pdf")), p, width = 15, height = 10)
  ggsave(file.path(OUT, paste0(tag, "_GO_sankey_bubble.png")), p, width = 15, height = 10, dpi = 300)
  cat(tag, "GO 桑葚气泡图已保存\n")
}

# ---- KEGG 弦图（从 term + geneID 数据框）----
chord_kegg_df <- function(terms, gsets, tag, topN = 15) {
  n <- length(terms)
  shared <- matrix(0, n, n)
  for (i in 1:n) for (j in 1:n) shared[i, j] <- length(intersect(gsets[[i]], gsets[[j]]))
  rownames(shared) <- colnames(shared) <- terms
  col_pal <- c("#6A5C95","#D99AAE","#E7A889","#A9CBB8","#CFE3E6","#8FA2D6","#5B608C","#C1A0C8",
               "#E8C9A0","#9AB0C8","#B5C9A0","#D6A0A8","#C8B0C8","#A0C8C8","#C8C8A0")
  path_colors <- setNames(rep(col_pal, length.out = n), terms)
  col_mat <- matrix(NA, n, n)
  for (i in 1:n) for (j in 1:n) col_mat[i, j] <- path_colors[i]
  for (dev in c("pdf","png")) {
    if (dev == "pdf") pdf(file.path(OUT, paste0(tag, "_KEGG_chord.pdf")), width = 11, height = 11)
    else png(file.path(OUT, paste0(tag, "_KEGG_chord.png")), width = 11, height = 11, units = "in", res = 300)
    par(mar = c(2,2,2,2))
    chordDiagram(shared, grid.col = path_colors, col = col_mat, transparency = 0.15,
                 annotationTrack = c("grid","axis"),
                 preAllocateTracks = list(track.height = 0.04, track.margin = c(0.01, 0.01)))
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      s <- get.cell.meta.data("sector.index")
      circos.text(mean(get.cell.meta.data("xlim")), get.cell.meta.data("ylim")[1] + 0.5, s,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6)
    }, bg.border = NA)
    dev.off()
  }
  cat(tag, "KEGG 弦图已保存\n")
}

# ---- 执行（支持命令行参数 179 / 6hub 分开运行，避免 OOM）----
args <- commandArgs(trailingOnly = TRUE)
which <- if (length(args) > 0) args[1] else "both"

if (which %in% c("both", "179")) {
  e179 <- go_enrich(cgs)
  sankey_go(e179, "179CGs")
  k179 <- read.csv(file.path(BASE, "03_intersection", "Shared_CGs_KEGG.csv"), stringsAsFactors = FALSE)
  k179 <- k179[order(k179$pvalue), ][1:min(15, nrow(k179)), ]
  chord_kegg_df(k179$Description, lapply(k179$geneID, function(s) strsplit(s, "/")[[1]]), "179CGs")
}
if (which %in% c("both", "6hub")) {
  e6 <- go_enrich(hub)
  sankey_go(e6, "6hub")
  k6 <- read.csv(file.path(BASE, "11_network_pharma", "KEGG_2021_Human_enrichment.csv"), stringsAsFactors = FALSE)
  k6 <- k6[order(k6$pvalue), ][1:min(15, nrow(k6)), ]
  chord_kegg_df(k6$term, lapply(k6$overlapping_genes, function(s) strsplit(s, ",")[[1]]), "6hub")
}
cat("\n完成 (", which, ")\n")
