# ============================================================
# Step 2g: GO 富集 - 桑葚气泡图(sankey bubble) + 弦图(chord)
# 两组: 179 shared CGs + 6 hub genes
# 参考: Kir2.1 (Nat Commun 2022) fig1i 微生信桑基气泡图
# ============================================================
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db)
  library(circlize); library(ggplot2); library(dplyr)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "AJRCCM_submission", "07_figure_exports", "Figure3_GO_sankey_chord")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- 1. 数据 ----
go179 <- read.csv(file.path(BASE, "03_intersection", "Shared_CGs_GO_BP.csv"), stringsAsFactors = FALSE)
hub <- readLines(file.path(BASE, "07_hub_immune", "hub_genes.txt"))
eg_hub <- bitr(hub, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ego6 <- enrichGO(gene = eg_hub$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                 pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
go6 <- as.data.frame(ego6)
cat("6 hub GO BP:", nrow(go6), "条; 179 CGs GO BP:", nrow(go179), "条\n")

make_long <- function(go, topN = 12, topGenes = 8) {
  go <- go[order(go$pvalue), ]
  go <- go[1:min(topN, nrow(go)), ]
  do.call(rbind, lapply(1:nrow(go), function(i){
    gs <- sort(strsplit(go$geneID[i], "/")[[1]])
    gs <- gs[1:min(topGenes, length(gs))]
    data.frame(gene = gs, term = go$Description[i],
               count = go$Count[i], neglogp = -log10(go$pvalue[i]),
               stringsAsFactors = FALSE)
  }))
}

# ---- 2. 桑葚气泡图（左 gene 流带 -> 右 GO term 气泡）----
sankey_bubble <- function(long, tag) {
  terms <- unique(long$term); genes <- unique(long$gene)
  long$term_f <- factor(long$term, levels = rev(terms))  # 最显著在上
  long$gene_f <- factor(long$gene, levels = rev(genes))
  long$y_term <- as.numeric(long$term_f)
  long$y_gene <- as.numeric(long$gene_f)
  ti <- long[!duplicated(long$term), c("term_f","y_term","count","neglogp")]
  gi <- long[!duplicated(long$gene), c("gene_f","y_gene")]

  col_pal <- c("#6A5C95","#D99AAE","#E7A889","#A9CBB8","#CFE3E6","#8FA2D6","#5B608C","#C1A0C8",
               "#E8C9A0","#9AB0C8","#B5C9A0","#D6A0A8")
  term_col <- setNames(rep(col_pal, length.out = length(terms)), rev(terms))

  p <- ggplot() +
    # gene -> term 流带（sigmoid）
    geom_curve(data = long,
               aes(x = 1, y = y_gene, xend = 1.8, yend = y_term, color = term),
               curvature = -0.18, alpha = 0.35, linewidth = 0.4) +
    # 左列 gene 标签（小圆点 + 文字）
    geom_point(data = gi, aes(x = 1, y = y_gene), size = 1.5, color = "grey40") +
    geom_text(data = gi, aes(x = 0.97, y = y_gene, label = gene_f),
              size = 2.6, hjust = 1, check_overlap = TRUE) +
    # 右列 GO term 气泡（size=count, fill=-log10P）
    geom_point(data = ti, aes(x = 1.8, y = y_term, size = count, fill = neglogp),
               shape = 21, color = "grey30", alpha = 0.9) +
    geom_text(data = ti, aes(x = 1.83, y = y_term, label = term_f),
              size = 2.8, hjust = 0) +
    scale_size_continuous(range = c(3, 16), name = "Gene count") +
    scale_fill_gradient(low = "#CFE3E6", high = "#6A5C95", name = "-log10(P)") +
    scale_color_manual(values = term_col, guide = "none") +
scale_x_continuous(limits = c(0.2, 4.5)) +
  theme_void(base_size = 9) +
  theme(legend.position = "right", plot.title = element_text(size = 11, face = "bold", hjust = 0.5))
ggsave(file.path(OUT, paste0(tag, "_sankey_bubble.pdf")), p, width = 15, height = 9)
ggsave(file.path(OUT, paste0(tag, "_sankey_bubble.png")), p, width = 15, height = 9, dpi = 300)
  cat(tag, "桑葚气泡图已保存\n")
}

# ---- 3. 弦图（GO term 间 shared genes）----
chord_plot <- function(go, tag, topN = 12) {
  go <- go[order(go$pvalue), ][1:min(topN, nrow(go)), ]
  terms <- go$Description; n <- length(terms)
  gsets <- lapply(go$geneID, function(s) strsplit(s, "/")[[1]])
  shared <- matrix(0, n, n)
  for (i in 1:n) for (j in 1:n) shared[i, j] <- length(intersect(gsets[[i]], gsets[[j]]))
  rownames(shared) <- colnames(shared) <- terms
  col_pal <- c("#6A5C95","#D99AAE","#E7A889","#A9CBB8","#CFE3E6","#8FA2D6","#5B608C","#C1A0C8",
               "#E8C9A0","#9AB0C8","#B5C9A0","#D6A0A8")
  path_colors <- setNames(rep(col_pal, length.out = n), terms)
  col_mat <- matrix(NA, n, n)
  for (i in 1:n) for (j in 1:n) col_mat[i, j] <- path_colors[i]
  for (dev in c("pdf","png")) {
    if (dev == "pdf") pdf(file.path(OUT, paste0(tag, "_chord.pdf")), width = 11, height = 11)
    else png(file.path(OUT, paste0(tag, "_chord.png")), width = 11, height = 11, units = "in", res = 300)
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
  cat(tag, "弦图已保存\n")
}

# ---- 4. 执行 4 张图 ----
sankey_bubble(make_long(go179, 12), "179CGs")
chord_plot(go179, "179CGs", 12)
sankey_bubble(make_long(go6, 12), "6hub")
chord_plot(go6, "6hub", 12)
cat("\n完成\n")
