# ---------------------------------------------------------------------------
# 25_supp_S6.R -- Figure S6. Single-gene GSEA of the six hub genes, four cohorts
#   tile grid of the two leading hallmark pathways per gene and cohort, with
#   normalised enrichment scores and adjusted P values
# Input  : data/frozen/12_single_gene_gsea/single_gene_gsea_top2_4cohorts.csv
# Output : Figures/supplementary/FigureS6.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

g <- read.csv(file.path("data", "frozen", "12_single_gene_gsea",
                        "single_gene_gsea_top2_4cohorts.csv"), stringsAsFactors = FALSE)

coh <- c(GSE76293 = "ARDS\ndiscovery", GSE69528 = "Sepsis\ndiscovery",
         GSE243066 = "ARDS\nvalidation", GSE100159 = "Sepsis\nvalidation")
g$coh <- factor(coh[g$cohort], levels = coh)
g$gene <- factor(g$hub_gene, levels = HUB6)
g$pw <- sub("^HALLMARK_", "", g$pathway)
g$pw <- gsub("_", " ", g$pw)
g$pw <- tolower(g$pw)
g$pw <- paste0(toupper(substr(g$pw, 1, 1)), substr(g$pw, 2, nchar(g$pw)))
g$sig <- ifelse(g$padj < 0.05, "*", "")

# order pathways within each gene by their mean enrichment across cohorts
ord <- aggregate(NES ~ gene + pw, data = g, FUN = mean)
ord <- ord[order(ord$gene, -ord$NES), ]
g$pw <- factor(g$pw, levels = rev(unique(ord$pw)))

ggplot(g, aes(coh, pw)) +
  geom_tile(aes(fill = NES), colour = "white", linewidth = 0.45) +
  geom_text(aes(label = sig), family = FONT_BODY, size = 9 / .pt,
            colour = pal[["ink"]], vjust = 0.78) +
  facet_grid(gene ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradientn(colours = ramp_div, limits = c(-3, 3),
                       oob = scales::squish, name = "NES") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_ajr() +
  theme(axis.text.x = element_text(lineheight = 0.9, size = PT_BODY),
          axis.text.y = element_text(size = PT_BODY),
          axis.line = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, size = PT_BODY,
                                           face = "italic", hjust = 1),
          strip.background = element_blank(),
          panel.spacing.y = unit(1.4, "mm"),
          legend.key.height = unit(6, "mm"), legend.key.width = unit(2.2, "mm"),
          legend.title = element_text(size = PT_BODY),
          plot.margin = margin(2, 3, 2, 3))

fig <- ggplot2::last_plot()
save_fig(fig, "FigureS6", dir = out_dir("supplementary"), height_mm = 208)
msg("FigureS6 written")
