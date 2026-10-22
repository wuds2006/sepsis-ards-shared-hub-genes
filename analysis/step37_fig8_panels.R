#!/usr/bin/env Rscript
# Figure 8 panels (drug prediction + molecular docking, ggplot2 + ggsci NPG):
#   A: docking affinity heatmap (6 hub genes x 3 natural compounds)
#   C: CMap/DSigDB drug prediction barplot (top terms by combined score)
# (Panel B PyMOL binding-pose renders are kept from the original pipeline.)
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

DD <- file.path(PROJ, "09_docking_drug")

## A: docking heatmap
dk <- read.csv(file.path(DD, "docking_results.csv"))
dk$gene <- factor(dk$gene, levels = rev(HUB6))
dk$ligand <- factor(dk$ligand, levels = c("EGCG", "quercetin", "tetrahydroalstonine"))
pA <- ggplot(dk, aes(ligand, gene, fill = affinity)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", affinity)), size = 2.2) +
  scale_fill_gradient(low = "#E64B35FF", high = "#3C5488FF",
                      name = "Affinity\n(kcal/mol)") +
  scale_y_discrete(labels = function(x) parse(text = gene_italic(x))) +
  labs(title = "A  Molecular docking (AutoDock Vina)", x = NULL, y = NULL) +
  theme_ajrccm() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.title = element_text(size = 6), legend.key.height = unit(4, "mm"))
save_pdf(pA, file.path(OUT_MAIN, "Fig8A_docking_heatmap.pdf"), 95, 68)

## C: drug prediction barplots (Old CMap down + up, DSigDB)
mk_bar <- function(f, ttl, top_n = 10) {
  x <- read.csv(f)
  x <- x[order(x$combined_score, decreasing = TRUE), ][seq_len(min(top_n, nrow(x))), ]
  x$term <- gsub("-\\d+$| CTD \\d+$| HL60 (UP|DOWN)$", "", x$term)
  x <- x[order(x$combined_score), ]
  x$term <- factor(x$term, levels = x$term)
  ggplot(x, aes(term, combined_score)) +
    geom_col(fill = "#3C5488FF", width = 0.7) + coord_flip() +
    labs(title = ttl, x = NULL, y = "Combined score") +
    theme_ajrccm() + theme(axis.text.y = element_text(size = 5.5))
}
pC1 <- mk_bar(file.path(DD, "Old_CMAP_down_hub_prediction.csv"), "C1  CMap (down)")
pC2 <- mk_bar(file.path(DD, "Old_CMAP_up_hub_prediction.csv"), "C2  CMap (up)")
pC3 <- mk_bar(file.path(DD, "DSigDB_hub_prediction.csv"), "C3  DSigDB")
save_pdf(pC1 + pC2 + pC3 + plot_layout(nrow = 1),
         file.path(OUT_MAIN, "Fig8C_drug_prediction.pdf"), 178, 68)

cat("Figure 8 panels done.\n")
