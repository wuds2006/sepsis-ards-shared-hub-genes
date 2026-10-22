# ============================================================
# Step 2f: GO/KEGG 可视化 - 桑葚图 (sankey/alluvial) + 弦图 (chord)
# 数据: 11_network_pharma/KEGG_2021_Human_enrichment.csv
# ============================================================
suppressPackageStartupMessages({
  library(ggalluvial); library(circlize); library(ggplot2); library(dplyr)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "AJRCCM_submission", "07_figure_exports", "Figure3_GO_sankey_chord")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# 读 KEGG 数据(Enrichr 格式)
kegg <- read.csv(file.path(BASE, "11_network_pharma/KEGG_2021_Human_enrichment.csv"),
                 stringsAsFactors = FALSE)
cat("KEGG pathways:", nrow(kegg), "\n")

# 取 top 12 个 pathway
kegg <- kegg[1:min(12, nrow(kegg)), ]
cat("使用 top", nrow(kegg), "pathways\n")

# ---- 拆 overlapping_genes 为 long format ----
parse_genes <- function(s) trimws(strsplit(s, ",")[[1]])
long <- do.call(rbind, lapply(1:nrow(kegg), function(i){
  data.frame(gene    = parse_genes(kegg$overlapping_genes[i]),
             pathway = kegg$term[i],
             pvalue  = kegg$pvalue[i],
             neglogp = -log10(kegg$pvalue[i]),
             zscore  = kegg$zscore[i],
             stringsAsFactors = FALSE)
}))
long$pathway <- factor(long$pathway, levels = rev(kegg$term))
long$gene    <- factor(long$gene, levels = rev(unique(long$gene)))

# ---- 桑葚图 (ggalluvial) ----
# 三层: Gene -> Pathway -> (右列节点用 -log10P)
p_sankey <- ggplot(long, aes(axis1 = gene, axis2 = pathway, y = neglogp)) +
  geom_alluvium(aes(fill = pathway), width = 0.25, alpha = 0.85, color = "white", linewidth = 0.1) +
  geom_stratum(width = 0.25, fill = "white", color = "grey40", linewidth = 0.2) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.6, check_overlap = TRUE) +
  scale_x_discrete(limits = c("Gene", "KEGG Pathway"), expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = rep(c("#6A5C95","#D99AAE","#E7A889","#A9CBB8","#CFE3E6","#8FA2D6","#5B608C"), 3)) +
  theme_classic(base_size = 9) +
  theme(legend.position = "none",
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size = 9, face = "bold"),
        panel.grid = element_blank())
ggsave(file.path(OUT, "Figure3_Kegg_Sankey.pdf"), p_sankey, width = 14, height = 9)
ggsave(file.path(OUT, "Figure3_Kegg_Sankey.png"), p_sankey, width = 14, height = 9, dpi = 300)
cat("桑葚图已保存\n")

# ---- 弦图 (circlize) ----
# 计算 pathway-pathway 的 shared gene 数
terms <- kegg$term
n <- nrow(kegg)
shared <- matrix(0, n, n)
for (i in 1:n) {
  gi <- parse_genes(kegg$overlapping_genes[i])
  for (j in i:n) {
    gj <- parse_genes(kegg$overlapping_genes[j])
    s <- length(intersect(gi, gj))
    shared[i, j] <- s; shared[j, i] <- s
  }
}
rownames(shared) <- colnames(shared) <- terms
cat("shared 矩阵 (前 5x5):\n"); print(shared[1:5, 1:5])

# 配色 (与桑葚图一致)
ncol_pal <- 7
col_pal <- c("#6A5C95","#D99AAE","#E7A889","#A9CBB8","#CFE3E6","#8FA2D6","#5B608C")
path_colors <- setNames(rep(col_pal, length.out = n), terms)

col_mat <- matrix(NA, n, n)
for (i in 1:n) for (j in 1:n) col_mat[i, j] <- path_colors[i]

# 弦图
pdf(file.path(OUT, "Figure3_Kegg_Chord.pdf"), width = 9, height = 9)
par(mar = c(2, 2, 2, 2))
chordDiagram(shared,
             grid.col = path_colors,
             col = col_mat,
             transparency = 0.15,
             annotationTrack = c("grid", "axis"),
             preAllocateTracks = list(
               track.height = 0.04,
               track.margin = c(0.01, 0.01)
             ))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1] + 0.5, sector.name,
              facing = "clockwise", niceFacing = TRUE,
              adj = c(0, 0.5), cex = 0.7)
}, bg.border = NA)
dev.off()
png(file.path(OUT, "Figure3_Kegg_Chord.png"), width = 9, height = 9, units = "in", res = 300)
par(mar = c(2, 2, 2, 2))
chordDiagram(shared,
             grid.col = path_colors,
             col = col_mat,
             transparency = 0.15,
             annotationTrack = c("grid", "axis"),
             preAllocateTracks = list(track.height = 0.04, track.margin = c(0.01, 0.01)))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim"); ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1] + 0.5, sector.name,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.7)
}, bg.border = NA)
dev.off()
cat("弦图已保存\n")
cat("\n完成，输出到:", OUT, "\n")
