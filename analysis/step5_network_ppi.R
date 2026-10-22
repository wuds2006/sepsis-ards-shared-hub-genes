# Step: PPI 网络 + 通路富集（STRING + clusterProfiler）
suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

setwd("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis")
OUT <- "11_network_pharma"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

HUB <- c("CDK5RAP2","CD177","GRB10","UPP1","EXOSC4","VNN1")

# ---- 1. 读 STRING 网络 ----
net <- read.delim(file.path(OUT, "string_functional.tsv"), stringsAsFactors = FALSE)
cat("网络边数:", nrow(net), "\n")
edges <- net[, c("preferredName_A","preferredName_B")]
edges$score <- as.numeric(net$score)

# 去重 + 自环
edges <- edges[edges$preferredName_A != edges$preferredName_B, ]
g <- graph_from_data_frame(edges[,1:2], directed = FALSE)
E(g)$score <- edges$score
nodes <- V(g)$name
cat("网络节点数:", length(nodes), "\n")

# ---- 2. 度数中心性（识别 hub 节点）----
deg <- degree(g)
deg_sorted <- sort(deg, decreasing = TRUE)
cat("\n===== 度数 Top15 节点 =====\n")
for (i in 1:min(15, length(deg_sorted))) {
  cat(sprintf("  %-12s degree=%d\n", names(deg_sorted)[i], deg_sorted[i]))
}

# 6 个 Hub 基因的度数
cat("\n===== 6 个 Hub 基因度数 =====\n")
for (h in HUB) {
  if (h %in% names(deg)) cat(sprintf("  %-12s degree=%d\n", h, deg[h])) else cat(sprintf("  %-12s 不在网络\n", h))
}

# ---- 3. 富集基因集: 6 Hub + 直接邻居 ----
hub_idx <- which(V(g)$name %in% HUB)
neigh_idx <- unique(unlist(lapply(hub_idx, function(i) as.integer(neighbors(g, i)))))
neighbors <- V(g)$name[neigh_idx]
enrich_genes <- unique(c(HUB, neighbors))
cat("\n富集基因集大小(6 Hub + 直接邻居):", length(enrich_genes), "\n")

# 符号 -> ENTREZ
eg <- bitr(enrich_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat("映射到 ENTREZ 的基因数:", nrow(eg), "\n")

# ---- 4. GO 富集 ----
go_bp <- enrichGO(gene = eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                  readable = TRUE)
if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
  go_df <- as.data.frame(go_bp)
  write.csv(go_df, file.path(OUT, "GO_BP_enrichment.csv"), row.names = FALSE)
  cat("\n===== GO BP Top15 =====\n")
  for (i in 1:min(15, nrow(go_df))) {
    cat(sprintf("  %-55s p.adjust=%.2e 基因数=%d\n", go_df$Description[i], go_df$p.adjust[i], go_df$Count[i]))
  }
  # 点图
  p_go <- dotplot(go_bp, showCategory = 15) + ggtitle("GO Biological Process Enrichment") + theme_minimal()
  ggsave(file.path(OUT, "GO_BP_dotplot.png"), p_go, width = 9, height = 7, dpi = 300)
  ggsave(file.path(OUT, "GO_BP_dotplot.pdf"), p_go, width = 9, height = 7)
} else {
  cat("GO BP 无显著富集\n")
}

# ---- 5. KEGG 富集（rest.kegg.jp 可能不通，用 tryCatch 保护，失败则用 Enrichr 补做）----
kegg <- tryCatch(
  enrichKEGG(gene = eg$ENTREZID, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
  error = function(e) { cat("KEGG API 失败:", conditionMessage(e), "（改用 Enrichr 补做）\n"); NULL })
if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
  kegg_df <- as.data.frame(kegg)
  write.csv(kegg_df, file.path(OUT, "KEGG_enrichment.csv"), row.names = FALSE)
  cat("\n===== KEGG Top15 =====\n")
  for (i in 1:min(15, nrow(kegg_df))) {
    cat(sprintf("  %-55s p.adjust=%.2e 基因数=%d\n", kegg_df$Description[i], kegg_df$p.adjust[i], kegg_df$Count[i]))
  }
  p_kegg <- dotplot(kegg, showCategory = 15) + ggtitle("KEGG Pathway Enrichment") + theme_minimal()
  ggsave(file.path(OUT, "KEGG_dotplot.png"), p_kegg, width = 9, height = 7, dpi = 300)
  ggsave(file.path(OUT, "KEGG_dotplot.pdf"), p_kegg, width = 9, height = 7)
} else {
  cat("KEGG 无显著富集（或 API 失败）\n")
}

# ---- 6. 网络可视化（igraph + ggplot2）----
set.seed(42)
layout_fr <- layout_with_fr(g, weights = E(g)$score)
lay <- data.frame(x = layout_fr[,1], y = layout_fr[,2], name = V(g)$name, stringsAsFactors = FALSE)
lay$degree <- degree(g)
lay$is_hub <- lay$name %in% HUB
lay$hub_label <- ifelse(lay$is_hub, lay$name, "")

# 边
edge_df <- data.frame(
  x = layout_fr[as.numeric(as.factor(edges$preferredName_A)), 1],
  y = layout_fr[as.numeric(as.factor(edges$preferredName_A)), 2],
  xend = layout_fr[as.numeric(as.factor(edges$preferredName_B)), 1],
  yend = layout_fr[as.numeric(as.factor(edges$preferredName_B)), 2],
  score = edges$score
)

p_net <- ggplot() +
  geom_segment(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend, alpha = score),
               color = "grey70", linewidth = 0.3) +
  geom_point(data = lay, aes(x = x, y = y, size = degree, color = is_hub)) +
  scale_color_manual(values = c("TRUE" = "#D9534F", "FALSE" = "#5B7FA6")) +
  scale_size_continuous(range = c(2, 10)) +
  ggrepel::geom_text_repel(data = lay[lay$is_hub, ], aes(x = x, y = y, label = hub_label),
                            size = 4, fontface = "bold", color = "#D9534F", max.overlaps = 20) +
  theme_void() + theme(legend.position = "none") +
  ggtitle("PPI Network: 6 Hub Genes + STRING Neighbors (functional)")
ggsave(file.path(OUT, "PPI_network.png"), p_net, width = 9, height = 8, dpi = 300)
ggsave(file.path(OUT, "PPI_network.pdf"), p_net, width = 9, height = 8)

# 保存网络边 + 节点属性
write.csv(edges, file.path(OUT, "PPI_edges.csv"), row.names = FALSE)
node_df <- data.frame(gene = lay$name, degree = lay$degree, is_hub = lay$is_hub, stringsAsFactors = FALSE)
write.csv(node_df, file.path(OUT, "PPI_nodes.csv"), row.names = FALSE)

cat("\n===== PPI 网络 + 富集完成 =====\n")
