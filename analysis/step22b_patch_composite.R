# ============================================================
# Step 22b — 仅重拼 FigureS6 composite（不重跑富集，避免 KEGG 在线查询）
# 修正: 顶部 6 个 bar 面板标题互相遮挡 -> 标题缩短为 KO: GENE，
#       队列信息移入 plot_annotation caption
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ragg)
})

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "08_singlecell/GSE167363/knockout")

PURPLE <- "#6A5C95"; PINK <- "#D99AAE"; ORANGE <- "#E7A889"
GREEN  <- "#A9CBB8"; MINT <- "#CFE3E6"; BLUE <- "#8FA2D6"; NAVY <- "#5B608C"

theme_ko <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(axis.title = element_text(size = 8.5),
          axis.text  = element_text(size = 7.5),
          plot.title = element_text(size = 9, face = "bold", hjust = 0),
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
  "PBMC|UPP1"      = file.path(BASE, "08_singlecell/knockout/KO_UPP1_diffRegulation.csv"),
  "PBMC|EXOSC4"    = file.path(BASE, "08_singlecell/knockout/KO_EXOSC4_diffRegulation.csv"),
  "CD14+|UPP1"     = file.path(OUT, "KO_UPP1_diffRegulation.csv"),
  "CD14+|CDK5RAP2" = file.path(OUT, "KO_CDK5RAP2_diffRegulation.csv"),
  "CD14+|EXOSC4"   = file.path(OUT, "KO_EXOSC4_diffRegulation.csv"),
  "CD14+|VNN1"     = file.path(OUT, "KO_VNN1_diffRegulation.csv")
)
top10_list <- list()
for (k in names(ko_sets)) {
  d <- read.csv(ko_sets[[k]], stringsAsFactors = FALSE)
  d <- d[order(d$distance, decreasing = TRUE), ]
  gene_ko <- sub(".*\\|", "", k)
  t10 <- head(d[d$gene != gene_ko, ], 10)
  t10$cohort <- sub("\\|.*", "", k); t10$KO <- gene_ko
  top10_list[[k]] <- t10
}

bar_panel <- function(k) {
  d <- top10_list[[k]]
  d$gene <- factor(d$gene, levels = rev(d$gene))
  d$dir  <- ifelse(d$FC > 0, "up", "down")
  d$sig  <- ifelse(d$p.adj < 0.05, "*", "")
  ggplot(d, aes(x = FC, y = gene, fill = dir)) +
    geom_col(width = 0.72, colour = NA) +
    geom_text(aes(label = sig, x = FC), hjust = -0.25, size = 3.4,
              colour = "#2b2b2b") +
    scale_fill_manual(values = c(up = PINK, down = BLUE)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(title = paste0("KO: ", d$KO[1]),
         x = "log2 FC (edge weight)", y = NULL) +
    theme_ko() +
    theme(axis.text.y = element_text(face = "italic"),
          legend.position = "none")
}

go_all   <- read.csv(file.path(OUT, "KO_top10_GO_BP.csv"), stringsAsFactors = FALSE)
kegg_all <- read.csv(file.path(OUT, "KO_top10_KEGG.csv"), stringsAsFactors = FALSE)

dot_panel <- function(df, db_label, ko_order) {
  df$KOp <- factor(paste0(df$cohort, " KO: ", df$KO), levels = ko_order)
  df <- df[order(df$p.adjust), ]
  df$Description <- factor(df$Description, levels = rev(unique(df$Description)))
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

lay <- "
ABCDEF
GGGGHH
"
comp <- (bar_panel("CD14+|UPP1") + bar_panel("CD14+|CDK5RAP2") +
         bar_panel("CD14+|EXOSC4") + bar_panel("CD14+|VNN1") +
         bar_panel("PBMC|UPP1") + bar_panel("PBMC|EXOSC4") +
         p_go + p_kegg) +
  plot_layout(design = lay, heights = c(1.15, 1.6)) +
  plot_annotation(
    tag_levels = "A",
    caption = paste0("Panels A-D: CD14+ monocytes (GSE167363). ",
                     "Panels E-F: PBMCs (GSE151263). ",
                     "Bars: log2 fold change of network edge weight (scTenifoldKnk); ",
                     "* FDR < 0.05. Dots: over-representation analysis of each KO's top-10 ",
                     "affected genes (clusterProfiler, BH-adjusted).")) &
  theme(plot.tag = element_text(size = 13, face = "bold", colour = "#2b2b2b"),
        plot.caption = element_text(size = 7, colour = "#555555", hjust = 0))

save_fig(comp, "FigureS6_KO_top10_enrichment", w = 17, h = 10.5,
         png_dpi = 500, tiff_dpi = 500)
cat("===== Step 22b 完成 =====\n")
