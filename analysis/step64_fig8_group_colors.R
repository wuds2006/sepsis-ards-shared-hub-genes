#!/usr/bin/env Rscript
# Step 64 (subject-reviewer revision): Figure 8 group-color contrast.
# Reviewer request (minor comments): "图8的配色方案：ARDS（红色）和脓毒症（蓝色）之间的对比可以更强烈。"
# Change: trajectory/AUCell violin panels now use explicit group colors with stronger contrast
#   Sepsis = NPG navy (#3C5488), ARDS = NPG red (#E64B35), HC = NPG grey-blue (#8491B4)
# Panels affected: Fig7C_trajectory_AUCell_GSE151263.pdf (panels C1-C3) and
#                  Fig7G_trajectory_AUCell_GSE167363.pdf (panels G1-G3).
# All statistics, layouts, and other Figure 8 panels are unchanged.
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
suppressPackageStartupMessages({ library(dplyr) })

SC <- file.path(PROJ, "08_singlecell")

COL_GROUP <- c(Sepsis = "#3C5488FF", ARDS = "#E64B35FF", HC = "#8491B4FF")

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
    scale_fill_manual(values = COL_GROUP[groups]) +
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
    scale_fill_manual(values = COL_GROUP[groups]) +
    labs(title = substitute(paste(l3, "  ", italic("UPP1"), " program (AUCell)"),
                            list(l3 = paste0(letter, "3"))),
         x = NULL, y = "AUCell score") +
    theme_ajrccm() + theme(legend.position = "none")
  cat(sprintf("%s UPP1 program medians: %s | wilcox P = %.4g\n",
              tag, paste(sprintf("%.3f", tapply(acm$UPP1_program, acm$group, median)),
                         collapse = " vs "), p3v))
  list(p1 = p1, p2 = p2, p3 = p3)
}

c151 <- traj_panels("GSE151263",
                    file.path(SC, "monocle3/monocle3_pseudotime_GSE151263.csv"),
                    file.path(SC, "cellchat_input/umap.csv"), c("UMAP1", "UMAP2"),
                    file.path(SC, "aucell/GSE151263_AUCell_scores.csv"),
                    c("Sepsis", "ARDS"), "C")
save_pdf(c151$p1 + c151$p2 + c151$p3 + plot_layout(nrow = 1),
         file.path(OUT_MAIN, "Fig7C_trajectory_AUCell_GSE151263.pdf"), 178, 62)

g167 <- traj_panels("GSE167363",
                    file.path(SC, "monocle3/monocle3_pseudotime_GSE167363.csv"),
                    file.path(SC, "GSE167363/input/umap.csv"), c("UMAP1", "UMAP2"),
                    file.path(SC, "GSE167363/aucell/GSE167363_AUCell_scores.csv"),
                    c("HC", "Sepsis"), "G")
save_pdf(g167$p1 + g167$p2 + g167$p3 + plot_layout(nrow = 1),
         file.path(OUT_MAIN, "Fig7G_trajectory_AUCell_GSE167363.pdf"), 178, 62)

cat("Step 64 done.\n")
