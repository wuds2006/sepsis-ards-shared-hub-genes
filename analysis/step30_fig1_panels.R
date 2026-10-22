#!/usr/bin/env Rscript
# Figure 1 panels (ggplot2 + ggsci NPG): B volcano, C module-trait, D 4-set Venn
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
suppressPackageStartupMessages({ library(dplyr); library(ggrepel); library(VennDiagram); library(grid) })

## ---------- Fig 1B: volcano plots ----------
mk_volcano <- function(deg_file, title_txt) {
  d <- read.csv(deg_file)
  d$grp <- ifelse(d$adj.P.Val < 0.05 & d$logFC > 1, "Up",
           ifelse(d$adj.P.Val < 0.05 & d$logFC < -1, "Down", "NS"))
  d$grp <- factor(d$grp, levels = c("Up", "Down", "NS"))
  d$mlog10 <- -log10(d$adj.P.Val)
  lab <- d[d$gene %in% HUB6, ]
  n_up <- sum(d$grp == "Up"); n_dn <- sum(d$grp == "Down")
  ggplot(d, aes(logFC, mlog10, color = grp)) +
    geom_point(size = 0.6, alpha = 0.8, stroke = 0) +
    scale_color_manual(values = c(Up = COL_UP, Down = COL_DOWN, NS = COL_NS),
                       labels = c(paste0("Up (n = ", n_up, ")"),
                                  paste0("Down (n = ", n_dn, ")"), "NS")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3, color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3, color = "grey50") +
    geom_text_repel(data = lab, aes(label = gene), fontface = "italic",
                    size = 2.2, color = "black", max.overlaps = 20,
                    box.padding = 0.3, point.padding = 0.2, segment.size = 0.2) +
    labs(title = title_txt, x = expression(log[2]~"fold change"),
         y = expression(-log[10]~"adjusted"~italic(P)), color = NULL) +
    theme_ajrccm() + theme(legend.position = "top", legend.key.height = unit(2.5, "mm"))
}

pB1 <- mk_volcano(file.path(PROJ, "01_DEG/GSE76293_DEG_all.csv"),
                  "GSE76293 (ARDS vs control)")
pB2 <- mk_volcano(file.path(PROJ, "01_DEG/GSE69528_DEG_all.csv"),
                  "GSE69528 (sepsis vs control)")
save_pdf(pB1, file.path(OUT_MAIN, "Fig1B1_volcano_GSE76293.pdf"), 78, 72)
save_pdf(pB2, file.path(OUT_MAIN, "Fig1B2_volcano_GSE69528.pdf"), 78, 72)

## ---------- Fig 1C: WGCNA module-trait (disease association) ----------
mk_mt <- function(rds_file, title_txt) {
  mt <- readRDS(rds_file)$mt
  mt <- mt %>% arrange(cor) %>%
    mutate(module = factor(module, levels = module),
           sig = sig_stars(p), direction = ifelse(cor > 0, "Positive", "Negative"))
  ggplot(mt, aes(cor, module, fill = direction)) +
    geom_col(width = 0.65, color = "black", linewidth = 0.25) +
    geom_text(aes(label = paste0(sprintf("%.2f", cor), " ", sig),
                  hjust = ifelse(cor > 0, -0.08, 1.08)), size = 2.1) +
    scale_fill_manual(values = c(Positive = COL_DISEASE, Negative = COL_DOWN)) +
    scale_x_continuous(limits = c(-1.25, 1.35), breaks = seq(-1, 1, 0.5)) +
    labs(title = title_txt, x = "Module-trait correlation (disease vs control)", y = NULL,
         fill = "Direction") +
    theme_ajrccm() + theme(legend.position = "top")
}

pC1 <- mk_mt(file.path(PROJ, "02_WGCNA/GSE76293_wgcna.rds"), "GSE76293 (ARDS)")
pC2 <- mk_mt(file.path(PROJ, "02_WGCNA/GSE69528_wgcna.rds"), "GSE69528 (sepsis)")
save_pdf(pC1, file.path(OUT_MAIN, "Fig1C1_module_trait_GSE76293.pdf"), 70, 85)
save_pdf(pC2, file.path(OUT_MAIN, "Fig1C2_module_trait_GSE69528.pdf"), 70, 85)

## ---------- Fig 1D: 4-set Venn ----------
read_lines <- function(f) trimws(readLines(f))
deg_ards   <- read.csv(file.path(PROJ, "01_DEG/GSE76293_DEG.csv"))$gene
deg_sepsis <- read.csv(file.path(PROJ, "01_DEG/GSE69528_DEG.csv"))$gene
wg_ards    <- read_lines(file.path(PROJ, "02_WGCNA/GSE76293_disease_module_genes.txt"))
wg_sepsis  <- read_lines(file.path(PROJ, "02_WGCNA/GSE69528_disease_module_genes.txt"))
sets <- list(`ARDS DEGs` = deg_ards, `ARDS WGCNA` = wg_ards,
             `Sepsis DEGs` = deg_sepsis, `Sepsis WGCNA` = wg_sepsis)
n_shared <- length(Reduce(intersect, sets))
cat("4-set intersection n =", n_shared, "\n")

vp <- venn.diagram(sets, filename = NULL, disable.logging = TRUE,
                   fill = pal_npg(4), alpha = 0.55,
                   col = pal_npg(4), lwd = 0.8,
                   cex = 0.8, cat.cex = 0.8, cat.fontface = "bold",
                   # cat.pos/dist tuned 2026-09-18: Sepsis pair pushed apart horizontally,
                   # ARDS WGCNA lifted above its 989 count to avoid label overlap
                   cat.pos = c(-30, 32, -42, 42), cat.dist = c(0.25, 0.27, 0.18, 0.18),
                   fontfamily = "Helvetica", cat.fontfamily = "Helvetica",
                   margin = 0.27)
grDevices::pdf(file.path(OUT_MAIN, "Fig1D_venn_4sets.pdf"), width = 78/25.4, height = 72/25.4,
               family = "Helvetica", useDingbats = FALSE)
grid.draw(vp)
dev.off()

cat("Figure 1 panels done.\n")
