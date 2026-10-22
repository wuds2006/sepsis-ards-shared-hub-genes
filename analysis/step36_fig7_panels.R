#!/usr/bin/env Rscript
# Figure 7 panels (single-cell, ggplot2 + ggsci NPG):
#   A: GSE151263 UMAP by cell type
#   C: GSE151263 monocyte pseudotime (UMAP + violin) + UPP1-program AUCell violin
#   D: GSE151263 CellOracle TF-KO delta barplot
#   E: GSE167363 UMAP by cell type
#   G: GSE167363 monocyte pseudotime (UMAP + violin) + UPP1-program AUCell violin
#   H: GSE167363 CellOracle TF-KO delta barplot
# (Panels B, F CellChat network renders are kept from the original pipeline.)
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

SC <- file.path(PROJ, "08_singlecell")

umap_plot <- function(umap_f, meta_f, umap_cols, cell_col, title) {
  u <- read.csv(umap_f, row.names = 1, check.names = FALSE)
  m <- read.csv(meta_f, row.names = 1, check.names = FALSE)
  colnames(u) <- umap_cols
  ct <- m[[cell_col]][match(rownames(u), rownames(m))]
  d <- data.frame(u, celltype = ct)
  d <- d[!is.na(d$celltype), ]
  set.seed(1); d <- d[sample(nrow(d)), ]
  ggplot(d, aes(UMAP1, UMAP2, color = celltype)) +
    geom_point(size = 0.4, alpha = 0.85, stroke = 0) +
    scale_color_manual(values = pal_npg(length(unique(d$celltype)))) +
    guides(color = guide_legend(override.aes = list(size = 1.6, alpha = 1), ncol = 1)) +
    labs(title = title, color = NULL) +
    theme_ajrccm() + theme(legend.text = element_text(size = 5.5),
                           legend.key.size = unit(2.5, "mm"))
}

# Trajectory panels use monocle3 pseudotime (the source of the manuscript statistics)
# plotted on the full-dataset UMAP embedding.
traj_panels <- function(tag, m3_f, umap_f, umap_cols, aucell_f, groups, letter) {
  m3 <- read.csv(m3_f)
  u <- read.csv(umap_f, row.names = 1, check.names = FALSE)
  colnames(u) <- umap_cols
  d <- data.frame(u[match(m3$cell, rownames(u)), ],
                  pseudotime = m3$pseudotime, group = m3$group)
  d <- d[!is.na(d$pseudotime) & !is.na(d$group), ]
  d$group <- factor(d$group, levels = groups)
  colnames(d)[1:2] <- c("UMAP1", "UMAP2")

  set.seed(1); d1 <- d[sample(nrow(d)), ]
  p1 <- ggplot(d1, aes(UMAP1, UMAP2, color = pseudotime)) +
    geom_point(size = 0.4, alpha = 0.85, stroke = 0) +
    scale_color_gradient(low = "#4DBBD5FF", high = "#E64B35FF", name = "Pseudotime") +
    labs(title = paste0(letter, "1  ", tag, " monocle3 trajectory")) +
    theme_ajrccm() + theme(legend.key.height = unit(4, "mm"), legend.title = element_text(size = 6))

  pval <- suppressWarnings(wilcox.test(pseudotime ~ group, data = d, exact = FALSE)$p.value)
  med <- tapply(d$pseudotime, d$group, median)
  p2 <- ggplot(d, aes(group, pseudotime, fill = group)) +
    geom_violin(linewidth = 0.3, scale = "width") +
    geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.28, fill = "white") +
    annotate("text", x = 1.5, y = max(d$pseudotime) * 1.02,
             label = paste0("P = ", format.pval(pval, digits = 2)), size = 2.2) +
    scale_fill_manual(values = pal_npg(2)) +
    labs(title = paste0(letter, "2  Pseudotime by group"), x = NULL, y = "Pseudotime") +
    theme_ajrccm() + theme(legend.position = "none")
  cat(sprintf("%s pseudotime medians: %s | wilcox P = %.4g\n",
              tag, paste(sprintf("%.2f", med), collapse = " vs "), pval))

  ac <- read.csv(aucell_f, check.names = FALSE)
  mono_ct <- grep("Mono", unique(ac$celltype), value = TRUE)
  acm <- ac[ac$celltype %in% mono_ct, ]
  acm$group <- factor(acm$group, levels = groups)
  p3v <- suppressWarnings(wilcox.test(UPP1_program ~ group, data = acm, exact = FALSE)$p.value)
  p3 <- ggplot(acm, aes(group, UPP1_program, fill = group)) +
    geom_violin(linewidth = 0.3, scale = "width") +
    geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.28, fill = "white") +
    annotate("text", x = 1.5, y = max(acm$UPP1_program) * 1.02,
             label = paste0("P = ", format.pval(p3v, digits = 2)), size = 2.2) +
    scale_fill_manual(values = pal_npg(2)) +
    labs(title = substitute(paste(l3, "  ", italic("UPP1"), " program (AUCell)"),
                            list(l3 = paste0(letter, "3"))),
         x = NULL, y = "AUCell score") +
    theme_ajrccm() + theme(legend.position = "none")
  cat(sprintf("%s UPP1 program medians: %s | wilcox P = %.4g\n",
              tag, paste(sprintf("%.3f", tapply(acm$UPP1_program, acm$group, median)),
                         collapse = " vs "), p3v))
  list(p1 = p1, p2 = p2, p3 = p3)
}

ko_bar <- function(ko_f, tag, letter, top_n = 8) {
  ko <- read.csv(ko_f)
  ko <- ko[order(ko$delta), ]
  tf_keep <- ko$KO_TF[ko$target == "UPP1"][seq_len(min(top_n, sum(ko$target == "UPP1")))]
  ko <- ko[ko$KO_TF %in% tf_keep, ]
  ko$KO_TF <- factor(ko$KO_TF, levels = rev(unique(tf_keep)))
  ko$target <- factor(ko$target, levels = HUB6)
  ggplot(ko, aes(KO_TF, delta, fill = target)) +
    geom_col(position = position_dodge(0.75), width = 0.7, linewidth = 0.2, color = "grey30") +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    coord_flip() +
    scale_fill_manual(values = pal_npg(6), labels = function(x) parse(text = gene_italic(x))) +
    labs(title = paste0(letter, "  ", tag, " CellOracle TF knockout"), x = NULL,
         y = "Simulated expression change (delta)", fill = "Target") +
    theme_ajrccm() + theme(legend.text = element_text(size = 5.5))
}

## A: GSE151263 UMAP
pA <- umap_plot(file.path(SC, "cellchat_input/umap.csv"), file.path(SC, "cellchat_input/meta.csv"),
                c("UMAP1", "UMAP2"), "celltype", "A  GSE151263 PBMC")
save_pdf(pA, file.path(OUT_MAIN, "Fig7A_UMAP_GSE151263.pdf"), 100, 78)

## C: GSE151263 trajectory + AUCell
c151 <- traj_panels("GSE151263",
                    file.path(SC, "monocle3/monocle3_pseudotime_GSE151263.csv"),
                    file.path(SC, "cellchat_input/umap.csv"), c("UMAP1", "UMAP2"),
                    file.path(SC, "aucell/GSE151263_AUCell_scores.csv"),
                    c("Sepsis", "ARDS"), "C")
save_pdf(c151$p1 + c151$p2 + c151$p3 + plot_layout(nrow = 1),
         file.path(OUT_MAIN, "Fig7C_trajectory_AUCell_GSE151263.pdf"), 178, 62)

## D: GSE151263 CellOracle
pD <- ko_bar(file.path(SC, "celloracle/CellOracle_TF_KO_hub_effect.csv"), "GSE151263", "D")
save_pdf(pD, file.path(OUT_MAIN, "Fig7D_CellOracle_GSE151263.pdf"), 100, 68)

## E: GSE167363 UMAP
pE <- umap_plot(file.path(SC, "GSE167363/input/umap.csv"), file.path(SC, "GSE167363/input/meta.csv"),
                c("UMAP1", "UMAP2"), "celltype", "E  GSE167363 CD14+ monocytes")
save_pdf(pE, file.path(OUT_MAIN, "Fig7E_UMAP_GSE167363.pdf"), 100, 78)

## G: GSE167363 trajectory + AUCell
g167 <- traj_panels("GSE167363",
                    file.path(SC, "monocle3/monocle3_pseudotime_GSE167363.csv"),
                    file.path(SC, "GSE167363/input/umap.csv"), c("UMAP1", "UMAP2"),
                    file.path(SC, "GSE167363/aucell/GSE167363_AUCell_scores.csv"),
                    c("HC", "Sepsis"), "G")
save_pdf(g167$p1 + g167$p2 + g167$p3 + plot_layout(nrow = 1),
         file.path(OUT_MAIN, "Fig7G_trajectory_AUCell_GSE167363.pdf"), 178, 62)

## H: GSE167363 CellOracle
pH <- ko_bar(file.path(SC, "GSE167363/celloracle/CellOracle_TF_KO_hub_effect.csv"), "GSE167363", "H")
save_pdf(pH, file.path(OUT_MAIN, "Fig7H_CellOracle_GSE167363.pdf"), 100, 68)

cat("Figure 7 panels done.\n")
