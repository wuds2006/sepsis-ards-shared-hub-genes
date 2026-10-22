# ============================================================
# Step 22 — scTenifoldKnk 敲除 top10 基因 logFC + GO/KEGG 富集
# 输入: 08_singlecell/knockout (PBMC, UPP1/EXOSC4)
#       08_singlecell/GSE167363/knockout (CD14+ 单核, UPP1/CDK5RAP2/EXOSC4/VNN1)
# 输出: KO_top10_logFC_summary.csv / KO_top10_GO_KEGG.csv
#       单面板 + FigureS6_KO_top10_enrichment.{pdf,png,tiff}
# 配色: Nature 7 色; line art TIFF 1000dpi / 彩色 500dpi
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2); library(clusterProfiler); library(org.Hs.eg.db)
  library(patchwork); library(ragg)
})

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "08_singlecell/GSE167363/knockout")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

PURPLE <- "#6A5C95"; PINK <- "#D99AAE"; ORANGE <- "#E7A889"
GREEN  <- "#A9CBB8"; MINT <- "#CFE3E6"; BLUE <- "#8FA2D6"; NAVY <- "#5B608C"

theme_ko <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(axis.title = element_text(size = 8.5),
          axis.text  = element_text(size = 7.5),
          plot.title = element_text(size = 9, face = "bold", hjust = 0),
          plot.subtitle = element_text(size = 7, colour = "#555555"),
          legend.title = element_text(size = 7.5),
          legend.text  = element_text(size = 7),
          axis.line = element_line(colour = "#cccccc"))
}

save_fig <- function(p, name, w, h, png_dpi = 500, tiff_dpi = 500, dir = OUT) {
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h,
         device = "pdf")
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = png_dpi, bg = "white")
  agg_tiff(file.path(dir, paste0(name, ".tiff")), width = w, height = h,
           units = "in", res = tiff_dpi, compression = "lzw", background = "white")
  print(p); invisible(dev.off())
  invisible(p)
}

ko_sets <- list(
  "PBMC|UPP1"    = file.path(BASE, "08_singlecell/knockout/KO_UPP1_diffRegulation.csv"),
  "PBMC|EXOSC4"  = file.path(BASE, "08_singlecell/knockout/KO_EXOSC4_diffRegulation.csv"),
  "CD14+|UPP1"   = file.path(OUT, "KO_UPP1_diffRegulation.csv"),
  "CD14+|CDK5RAP2" = file.path(OUT, "KO_CDK5RAP2_diffRegulation.csv"),
  "CD14+|EXOSC4" = file.path(OUT, "KO_EXOSC4_diffRegulation.csv"),
  "CD14+|VNN1"   = file.path(OUT, "KO_VNN1_diffRegulation.csv")
)

top10_list <- list()
for (k in names(ko_sets)) {
  d <- read.csv(ko_sets[[k]], stringsAsFactors = FALSE)
  d <- d[order(d$distance, decreasing = TRUE), ]
  gene_ko <- sub(".*\\|", "", k)
  t10 <- head(d[d$gene != gene_ko, ], 10)      # 排除敲除基因自身
  t10$cohort <- sub("\\|.*", "", k)
  t10$KO     <- gene_ko
  top10_list[[k]] <- t10
}
top10_all <- do.call(rbind, top10_list)
write.csv(top10_all, file.path(OUT, "KO_top10_logFC_summary.csv"),
          row.names = FALSE)
cat("top10 汇总已写出\n")

# ------------------------------------------------------------- bar panels
bar_panel <- function(k) {
  d <- top10_list[[k]]
  d$gene <- factor(d$gene, levels = rev(d$gene))
  d$dir  <- ifelse(d$FC > 0, "up", "down")
  d$sig  <- ifelse(d$p.adj < 0.05, "*", "")
  gene_ko <- d$KO[1]; cohort <- d$cohort[1]
  ggplot(d, aes(x = FC, y = gene, fill = dir)) +
    geom_col(width = 0.72, colour = NA) +
    geom_text(aes(label = sig, x = FC), hjust = -0.25, size = 3.4,
              colour = "#2b2b2b") +
    scale_fill_manual(values = c(up = PINK, down = BLUE)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(title = paste0("KO: ", gene_ko),
         subtitle = paste0(cohort, " cohort, scTenifoldKnk"),
         x = "log2 FC (edge weight)", y = NULL) +
    theme_ko() +
    theme(axis.text.y = element_text(face = "italic"),
          legend.position = "none")
}
panel_names <- c("CD14+|UPP1", "CD14+|CDK5RAP2", "CD14+|EXOSC4", "CD14+|VNN1",
                 "PBMC|UPP1", "PBMC|EXOSC4")
bars <- lapply(panel_names, bar_panel)
names(bars) <- panel_names
# 单面板导出 (line art: TIFF 1000 dpi)
for (k in panel_names) {
  save_fig(bars[[k]], paste0("KO_", sub("\\|", "_", k), "_top10_logFC"),
           w = 4.6, h = 3.4, png_dpi = 500, tiff_dpi = 1000)
}

# ------------------------------------------------------------- GO / KEGG
run_enrich <- function(genes, k, kind = "GO") {
  out <- tryCatch({
    if (kind == "GO") {
      e <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                    ont = "BP", pAdjustMethod = "BH",
                    pvalueCutoff = 1, qvalueCutoff = 1)
    } else {
      ez <- suppressMessages(bitr(genes, fromType = "SYMBOL",
                                  toType = "ENTREZID", OrgDb = org.Hs.eg.db))
      e <- enrichKEGG(gene = ez$ENTREZID, organism = "hsa",
                      keyType = "ncbi-geneid", pAdjustMethod = "BH",
                      pvalueCutoff = 1, qvalueCutoff = 1)
    }
    if (is.null(e) || nrow(as.data.frame(e)) == 0) return(NULL)
    df <- as.data.frame(e)
    df <- df[order(df$p.adjust), ]
    df <- head(df, 5)
    df$KO <- sub(".*\\|", "", k); df$cohort <- sub("\\|.*", "", k)
    df
  }, error = function(e) { cat(kind, k, "失败:", conditionMessage(e), "\n"); NULL })
  out
}

# enrichKEGG keyType: 用 SYMBOL 更稳妥 (clusterProfiler 新版本支持 ncbi-geneid/uniprot)
go_list  <- lapply(names(top10_list), function(k)
  run_enrich(top10_list[[k]]$gene, k, "GO"))
kegg_list <- lapply(names(top10_list), function(k)
  run_enrich(top10_list[[k]]$gene, k, "KEGG"))
go_all   <- do.call(rbind, Filter(Negate(is.null), go_list))
kegg_all <- do.call(rbind, Filter(Negate(is.null), kegg_list))
cat("GO terms:", nrow(go_all), "| KEGG terms:", nrow(kegg_all), "\n")
write.csv(go_all,   file.path(OUT, "KO_top10_GO_BP.csv"),   row.names = FALSE)
write.csv(kegg_all, file.path(OUT, "KO_top10_KEGG.csv"),    row.names = FALSE)

dot_panel <- function(df, db_label, ko_order) {
  if (is.null(df) || nrow(df) == 0) {
    return(ggplot() + annotate("text", 0.5, 0.5, label = paste(db_label, ": 无结果"),
                               size = 3) + theme_void())
  }
  df$KOp <- factor(paste0(df$cohort, " KO: ", df$KO),
                   levels = ko_order)
  df <- df[order(df$p.adjust), ]
  df$Description <- factor(df$Description,
                           levels = rev(unique(df$Description)))
  df$Count <- as.numeric(sub("/.*", "", df$GeneRatio))
  ggplot(df, aes(x = KOp, y = Description)) +
    geom_point(aes(size = Count, colour = -log10(p.adjust))) +
    scale_colour_gradient(low = MINT, high = PURPLE,
                          name = expression(-log[10](P[adj]))) +
    scale_size_continuous(range = c(1.6, 5), name = "Count") +
    labs(title = db_label, x = NULL, y = NULL) +
    theme_ko() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7))
}
ko_order <- paste0(rep(c("CD14+", "PBMC"), c(4, 2)), " KO: ",
                   c("UPP1", "CDK5RAP2", "EXOSC4", "VNN1", "UPP1", "EXOSC4"))
p_go   <- dot_panel(go_all,   "GO biological process (top 5 / KO)", ko_order)
p_kegg <- dot_panel(kegg_all, "KEGG pathway (top 5 / KO)", ko_order)
save_fig(p_go,   "KO_top10_GO_dotplot",   w = 7.2, h = 4.6)
save_fig(p_kegg, "KO_top10_KEGG_dotplot", w = 7.2, h = 4.0)

# ------------------------------------------------------------- composite
lay <- "
ABCDEF
GGGGHH
"
mk <- function(k, ttl) {
  bars[[k]] + labs(title = ttl, subtitle = NULL)
}
pA <- mk("CD14+|UPP1",     "KO: UPP1 (CD14+ monocytes)")
pB <- mk("CD14+|CDK5RAP2", "KO: CDK5RAP2 (CD14+ monocytes)")
pC <- mk("CD14+|EXOSC4",   "KO: EXOSC4 (CD14+ monocytes)")
pD <- mk("CD14+|VNN1",     "KO: VNN1 (CD14+ monocytes)")
pE <- mk("PBMC|UPP1",      "KO: UPP1 (PBMCs)")
pF <- mk("PBMC|EXOSC4",    "KO: EXOSC4 (PBMCs)")

comp <- (pA + pB + pC + pD + pE + pF + p_go + p_kegg) +
  plot_layout(design = lay, heights = c(1.15, 1.6)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, face = "bold", colour = "#2b2b2b"))

save_fig(comp, "FigureS6_KO_top10_enrichment", w = 17, h = 10.5,
         png_dpi = 500, tiff_dpi = 500)
cat("\n===== Step 22 完成 =====\n")
