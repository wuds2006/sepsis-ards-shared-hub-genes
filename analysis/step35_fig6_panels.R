#!/usr/bin/env Rscript
# Figure 6 panels (GSVA top-20 immune-pathway heatmaps, NPG palette):
#   A: GSE76293 (ARDS)   B: GSE69528 (sepsis)
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
suppressPackageStartupMessages(library(pheatmap))

BASE <- PROJ
GD <- file.path(BASE, "07_hub_immune/gsva_top20")

clean_pw <- function(x) {
  x <- gsub("^GOBP_|^GOMF_|^GOCC_|^HALLMARK_|^KEGG_", "", x)
  x <- gsub("_", " ", x)
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

gsva_heat <- function(ds, dis, panel, w_mm = 150, h_mm = 110) {
  sc <- as.matrix(read.csv(file.path(GD, paste0(ds, "_GSVA_scores.csv")), row.names = 1))
  top <- read.csv(file.path(GD, paste0(ds, "_GSVA_top20.csv")))
  grp <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_grp.rds")))
  grp$group <- factor(grp$group, levels = c("Control", dis))
  grp <- grp[order(grp$group), , drop = FALSE]
  pw <- top$pathway
  pw <- pw[pw %in% rownames(sc)]
  samp <- intersect(grp$sample, colnames(sc))
  grp <- grp[match(samp, grp$sample), , drop = FALSE]
  m <- sc[pw, samp, drop = FALSE]
  mz <- t(scale(t(m)))
  mz[mz > 2.5] <- 2.5; mz[mz < -2.5] <- -2.5
  rownames(mz) <- clean_pw(rownames(mz))
  ann <- data.frame(Group = grp$group); rownames(ann) <- samp
  ann_col <- list(Group = setNames(c(COL_CONTROL, COL_DISEASE), c("Control", dis)))
  grDevices::pdf(file.path(OUT_MAIN, paste0("Fig6", panel, "_GSVA_top20_", ds, ".pdf")),
                 width = w_mm/25.4, height = h_mm/25.4, family = "Helvetica", useDingbats = FALSE)
  pheatmap(mz, cluster_cols = FALSE, cluster_rows = TRUE,
           annotation_col = ann, annotation_colors = ann_col,
           color = colorRampPalette(c("#3C5488FF", "white", "#E64B35FF"))(101),
           show_colnames = FALSE, fontsize = 6, fontsize_row = 5.5,
           border_color = NA, treeheight_row = 12,
           main = paste0(panel, "  ", ds, " GSVA top 20 (row z-score)"))
  dev.off()
  invisible(NULL)
}

gsva_heat("GSE76293", "ARDS", "A")
gsva_heat("GSE69528", "Sepsis", "B")
cat("Figure 6 panels done.\n")
