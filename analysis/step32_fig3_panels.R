#!/usr/bin/env Rscript
# Figure 3 panels (ggplot2 + ggsci NPG):
#   A: STRING PPI network (igraph)
#   B: GO BP dotplot (179 shared candidate genes)
#   C: GO sankey-bubble, 179 CGs (BP/MF/CC)
#   D: KEGG chord, 179 CGs
#   E: GO sankey-bubble, 6 hub genes (BP/CC)
#   F: KEGG chord, 6 hub genes
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
suppressPackageStartupMessages({
  library(dplyr); library(igraph); library(circlize)
  library(clusterProfiler); library(org.Hs.eg.db)
})

## ---------- A: PPI network ----------
nodes <- read.csv(file.path(PROJ, "11_network_pharma/PPI_nodes.csv"))
edges <- read.csv(file.path(PROJ, "11_network_pharma/PPI_edges.csv"))
g <- graph_from_data_frame(edges[, c("preferredName_A", "preferredName_B")], directed = FALSE,
                           vertices = nodes$gene)
deg <- nodes$degree[match(V(g)$name, nodes$gene)]
deg[is.na(deg)] <- 1
hub_flag <- nodes$is_hub[match(V(g)$name, nodes$gene)]
hub_flag[is.na(hub_flag)] <- FALSE
mod <- ifelse(V(g)$name %in% c("UPP1", "UMPS", "CMPK1", "NME1", "NME2", "DTYMK"), "pyrimidine",
       ifelse(V(g)$name %in% c("EXOSC4", "EXOSC1", "EXOSC2", "EXOSC3", "EXOSC5", "EXOSC6",
                               "EXOSC7", "EXOSC8", "EXOSC9", "EXOSC10", "DIS3", "DIS3L"), "exosome",
       ifelse(V(g)$name %in% c("VNN1", "VNN2", "VNN3", "PANTETHEINASE"), "pantetheinase", "other")))
mod[hub_flag] <- "hub"
set.seed(42)
lay <- layout_with_fr(g)
grDevices::pdf(file.path(OUT_MAIN, "Fig3A_PPI_network.pdf"), width = 120/25.4, height = 100/25.4,
               family = "Helvetica", useDingbats = FALSE)
par(mar = c(0.5, 0.5, 0.5, 0.5))
cols <- c(hub = "#E64B35FF", pyrimidine = "#3C5488FF", exosome = "#00A087FF",
          pantetheinase = "#F39B7FFF", other = "#8491B4B3")
plot(g, layout = lay, vertex.size = pmax(2.5, sqrt(deg) * 2.2),
     vertex.color = cols[mod], vertex.frame.color = "white",
     vertex.label = ifelse(hub_flag, V(g)$name, NA),
     vertex.label.cex = 0.62, vertex.label.color = "black",
     vertex.label.family = "Helvetica", vertex.label.font = 3,
     edge.color = adjustcolor("grey70", alpha.f = 0.5), edge.width = 0.4)
legend("bottomleft", bty = "n", cex = 0.6, pt.cex = 1.4,
       legend = c("Hub gene", "Pyrimidine metabolism", "RNA exosome", "Pantetheinase", "Other"),
       pch = 21, pt.bg = cols[c("hub", "pyrimidine", "exosome", "pantetheinase", "other")])
dev.off()

## ---------- B: GO BP dotplot (179 CGs) ----------
go <- read.csv(file.path(PROJ, "03_intersection/Shared_CGs_GO_BP.csv"))
go_top <- go %>% arrange(p.adjust) %>% slice_head(n = 20) %>%
  mutate(Description = factor(Description, levels = rev(Description)),
         GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
pB <- ggplot(go_top, aes(GeneRatio_num, Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "#E64B35FF", high = "#3C5488FF", trans = "log10",
                       name = expression(adjusted~italic(P))) +
  scale_size(range = c(1.2, 5), name = "Count") +
  labs(title = "GO biological-process enrichment",
       x = "Gene ratio", y = NULL) +
  theme_ajrccm() + theme(axis.text.y = element_text(size = 6))
save_pdf(pB, file.path(OUT_MAIN, "Fig3B_GO_BP_dotplot.pdf"), 120, 95)

## ---------- enrichment recompute (BP/MF/CC) ----------
cg179 <- readLines(file.path(PROJ, "03_intersection/Shared_CGs.txt"))
cg179 <- trimws(cg179[nzchar(trimws(cg179))])
ent179 <- bitr(cg179, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ent6   <- bitr(HUB6, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

ego179 <- lapply(c("BP", "MF", "CC"), function(ont)
  as.data.frame(enrichGO(ent179$ENTREZID, OrgDb = org.Hs.eg.db, ont = ont,
                         pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)))
names(ego179) <- c("BP", "MF", "CC")
ego6 <- lapply(c("BP", "CC"), function(ont)
  as.data.frame(enrichGO(ent6$ENTREZID, OrgDb = org.Hs.eg.db, ont = ont,
                         pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)))
names(ego6) <- c("BP", "CC")

## sankey-bubble: genes (left) -> top terms (right); links curved; bubble = Count, color = -log10(padj)
mk_sankey_bubble <- function(enrich_list, title_txt, top_n = 8) {
  d <- do.call(rbind, lapply(names(enrich_list), function(ont) {
    x <- enrich_list[[ont]]
    if (nrow(x) == 0) return(NULL)
    x <- x %>% arrange(p.adjust) %>% slice_head(n = top_n)
    x$ONTOLOGY <- ont
    x
  }))
  stopifnot(!is.null(d))
  d$term_lab <- paste0(d$Description, " (", d$ONTOLOGY, ")")
  d <- d %>% arrange(p.adjust)
  terms <- unique(d$term_lab)
  long <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
    genes <- trimws(strsplit(d$geneID[i], "/")[[1]])
    data.frame(gene = genes, term_lab = d$term_lab[i], ONTOLOGY = d$ONTOLOGY[i],
               Count = d$Count[i], neglogp = -log10(d$p.adjust[i]), stringsAsFactors = FALSE)
  }))
  genes <- unique(long$gene)
  gy <- setNames(seq(1, length(genes), length.out = max(length(genes), 2)), genes)
  ty <- setNames(seq(1, length(terms), length.out = max(length(terms), 2)), terms)
  long$x1 <- 0; long$y1 <- gy[long$gene]
  long$x2 <- 1; long$y2 <- ty[long$term_lab]
  ont_cols <- setNames(pal_npg(3), c("BP", "MF", "CC"))
  term_pt <- long %>% distinct(term_lab, ONTOLOGY, Count, neglogp) %>%
    mutate(x = 1, y = ty[term_lab])
  ggplot(long) +
    geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2, color = ONTOLOGY),
               curvature = -0.35, alpha = 0.25, linewidth = 0.25) +
    geom_point(data = data.frame(gene = genes, y = gy), aes(x = 0, y = y),
               size = 0.8, color = "grey40") +
    geom_point(data = term_pt, aes(x = x, y = y, size = Count, fill = neglogp),
               shape = 21, color = "black", stroke = 0.25) +
    geom_text(data = data.frame(gene = genes, y = gy), aes(x = -0.02, y = y, label = gene),
              hjust = 1, size = 1.6, fontface = "italic") +
    geom_text(data = term_pt, aes(x = 1.02, y = y, label = term_lab), hjust = 0, size = 1.9) +
    scale_color_manual(values = ont_cols, guide = "none") +
    scale_fill_gradient(low = "#3C5488FF", high = "#E64B35FF",
                        name = expression(-log[10]~adjusted~italic(P))) +
    scale_size(range = c(1.5, 6), name = "Count") +
    scale_x_continuous(limits = c(-0.45, 2.6), expand = c(0, 0)) +
    labs(title = title_txt) +
    theme_void(base_size = 7, base_family = "Helvetica") +
    theme(plot.title = element_text(size = 7, face = "bold", family = "Helvetica"),
          legend.text = element_text(size = 6), legend.title = element_text(size = 6),
          legend.key.size = unit(3, "mm"))
}
pC <- mk_sankey_bubble(ego179, "GO enrichment sankey–bubble (179 shared candidates)")
save_pdf(pC, file.path(OUT_MAIN, "Fig3C_GO_sankey_bubble_179.pdf"), 150, 110)
pE <- mk_sankey_bubble(ego6, "GO enrichment sankey–bubble (six hub genes)", top_n = 12)
save_pdf(pE, file.path(OUT_MAIN, "Fig3E_GO_sankey_bubble_6hub.pdf"), 130, 80)

## ---------- D/F: KEGG chord ----------
kegg179 <- read.csv(file.path(PROJ, "03_intersection/Shared_CGs_KEGG.csv"))
kegg179_top <- kegg179 %>% arrange(p.adjust) %>% slice_head(n = 12)
map179 <- setNames(ent179$SYMBOL, ent179$ENTREZID)

# 6-hub KEGG (recompute; allow no-hit fallback to pyrimidine terms only)
ek6 <- tryCatch(as.data.frame(enrichKEGG(ent6$ENTREZID, organism = "hsa",
                                         pvalueCutoff = 0.05, pAdjustMethod = "BH")),
                error = function(e) NULL)
if (!is.null(ek6) && nrow(ek6) > 0) {
  kegg6_top <- ek6 %>% arrange(p.adjust) %>% slice_head(n = 10)
} else {
  ek6u <- tryCatch(as.data.frame(enrichKEGG(ent6$ENTREZID, organism = "hsa",
                                            pvalueCutoff = 1, pAdjustMethod = "none", qvalueCutoff = 1)),
                   error = function(e) NULL)
  kegg6_top <- if (!is.null(ek6u)) ek6u %>% arrange(pvalue) %>% slice_head(n = 10) else NULL
}

mk_chord <- function(kegg_df, sym_map, title_txt, pdf_file, w = 95, h = 95) {
  terms <- kegg_df$Description
  gid_lists <- lapply(strsplit(kegg_df$geneID, "/"), function(v)
    unname(na.omit(sym_map[trimws(v)])))
  mat <- matrix(0, length(terms), length(terms), dimnames = list(terms, terms))
  for (i in seq_along(terms)) for (j in seq_along(terms)) {
    if (i != j) mat[i, j] <- length(intersect(gid_lists[[i]], gid_lists[[j]]))
  }
  keep <- rowSums(mat) + colSums(mat) > 0
  if (sum(keep) < 3) { cat("chord skipped (too few connected terms):", pdf_file, "\n"); return(invisible(NULL)) }
  mat <- mat[keep, keep]
  grDevices::pdf(pdf_file, width = w/25.4, height = h/25.4, family = "Helvetica", useDingbats = FALSE)
  circos.par(start.degree = 90, gap.degree = 4)
  chordDiagram(mat, grid.col = colorRampPalette(pal_npg(10))(nrow(mat)),
               transparency = 0.35, annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.12),
               link.sort = TRUE, link.decreasing = TRUE)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    nm <- get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.3,
                strwrap(nm, 18)[1], facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = 0.42, family = "Helvetica")
  }, bg.border = NA)
  title(title_txt, cex.main = 0.55, family = "Helvetica")
  circos.clear()
  dev.off()
}
mk_chord(kegg179_top, map179, "KEGG pathway–pathway chord (179 shared candidates)",
         file.path(OUT_MAIN, "Fig3D_KEGG_chord_179.pdf"), 110, 110)
if (!is.null(kegg6_top))
  mk_chord(kegg6_top, map179, "KEGG pathway–pathway chord (six hub genes)",
           file.path(OUT_MAIN, "Fig3F_KEGG_chord_6hub.pdf"), 100, 100)

cat("Figure 3 panels done.\n")
