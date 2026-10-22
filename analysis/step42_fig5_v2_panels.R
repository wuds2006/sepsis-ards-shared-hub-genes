#!/usr/bin/env Rscript
# step42: Figure 5 v2 panels (CIBERSORT immune infiltration, ggplot2 + ggsci NPG)
#   A:   six-hub-gene ssGSEA scores, disease vs control, two discovery cohorts (facets)
#   B-E: CIBERSORT per-sample stacked bar plots (four datasets, one PDF each)
#   F:   Spearman correlation heatmaps among 22 immune cell types (corrplot, 2x2)
#   G:   differential cell abundance violins, Wilcoxon stars (facet ggplot)
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(corrplot) })

BASE <- PROJ
CBD  <- file.path(BASE, "07_hub_immune", "cibersort")
IMM  <- file.path(BASE, "07_hub_immune")

DS_META <- list(
  list(id = "GSE76293",  dis = "ARDS",   tag = "ARDS discovery (PMN fractions)"),
  list(id = "GSE69528",  dis = "Sepsis", tag = "Sepsis discovery (whole blood)"),
  list(id = "GSE243066", dis = "ARDS",   tag = "ARDS validation"),
  list(id = "GSE100159", dis = "Sepsis", tag = "Sepsis validation")
)
PANEL_LETTER <- c(GSE76293 = "B", GSE69528 = "C", GSE243066 = "D", GSE100159 = "E")

CELL22 <- c("B cells naive","B cells memory","Plasma cells","T cells CD8",
  "T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated",
  "T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta",
  "NK cells resting","NK cells activated","Monocytes","Macrophages M0",
  "Macrophages M1","Macrophages M2","Dendritic cells resting","Dendritic cells activated",
  "Mast cells resting","Mast cells activated","Eosinophils","Neutrophils")
CELL_SHORT <- c("B naive","B memory","Plasma","CD8 T","CD4 T naive","CD4 T mem rest",
  "CD4 T mem act","Tfh","Tregs","Tgd","NK rest","NK act","Monocytes","M0","M1","M2",
  "DC rest","DC act","Mast rest","Mast act","Eosinophils","Neutrophils")
names(CELL_SHORT) <- CELL22

# 22 distinguishable colors: NPG(10) + AAAS(10) + 2
pal22 <- c(ggsci::pal_npg("nrc")(10), ggsci::pal_aaas("default")(10), "#7E6148FF", "#B09C85FF")
names(pal22) <- CELL22

# ---------- A: six-hub-gene ssGSEA (discovery cohorts only: GSE76293, GSE69528) ----------
ss <- do.call(rbind, lapply(DS_META[1:2], function(z) {
  df <- read.csv(file.path(IMM, paste0(z$id, "_hub6_ssGSEA.csv")))
  st <- read.csv(file.path(IMM, "hub6_ssGSEA_stats.csv"))
  p <- st$P[st$dataset == z$id]
  data.frame(dataset = paste0(z$id, " ", z$dis, " ",
                              ifelse(z$id %in% c("GSE76293", "GSE69528"), "discovery", "validation")),
             group = factor(df$group, levels = c("Control", "ARDS", "Sepsis")),
             score = df$score,
             plab = sprintf("italic(P) == %.2g", p))
}))
ds_levels <- unique(ss$dataset)
ss$dataset <- factor(ss$dataset, levels = ds_levels)
plab_df <- ss %>% distinct(dataset, plab) %>%
  mutate(y = max(ss$score) * 1.06, x = 1.5)

p5a <- ggplot(ss, aes(group, score, fill = group)) +
  geom_violin(color = "grey30", linewidth = 0.25, scale = "width", width = 0.8) +
  geom_boxplot(width = 0.14, linewidth = 0.25, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.05, size = 0.35, alpha = 0.5, color = "grey25") +
  geom_text(data = plab_df, aes(x, y, label = plab), inherit.aes = FALSE,
            parse = TRUE, size = 2.2) +
  facet_wrap(~dataset, nrow = 1, scales = "free_x") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  scale_fill_manual(values = c(Control = COL_CONTROL, ARDS = COL_DISEASE, Sepsis = COL_DISEASE)) +
  labs(title = "Six-hub-gene signature ssGSEA score", x = NULL,
       y = "ssGSEA enrichment score", fill = NULL) +
  theme_ajrccm() +
  theme(legend.position = "none", strip.text = element_text(size = 6))
save_pdf(p5a, file.path(OUT_MAIN, "Fig5A_hub6_ssGSEA.pdf"), 178, 50)

# ---------- B-E: CIBERSORT stacked bar plots ----------
for (z in DS_META) {
  res <- read.csv(file.path(CBD, paste0(z$id, "_CIBERSORT_results.csv")), check.names = FALSE)
  res <- res[order(res$group, res$sample), ]
  long <- pivot_longer(res, cols = all_of(CELL22), names_to = "cell", values_to = "frac")
  long$cell <- factor(long$cell, levels = rev(CELL22))
  p <- ggplot(long, aes(x = sample, y = frac, fill = cell)) +
    geom_col(width = 0.92, linewidth = 0) +
    scale_fill_manual(values = pal22, breaks = CELL22, labels = CELL_SHORT[CELL22],
                      name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    labs(title = paste0(z$id, " (", z$tag, ")"),
         x = NULL, y = "Cell fraction") +
    theme_ajrccm(base_size = 6) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = 4.6),
          legend.key.size = unit(2.2, "mm"),
          legend.margin = margin(-2, 0, 0, 0)) +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE))
  save_pdf(p, file.path(OUT_MAIN, paste0("Fig5", PANEL_LETTER[[z$id]], "_CIBERSORT_stacked_", z$id, ".pdf")), 88, 62)
}

# ---------- F: Spearman correlation heatmaps (corrplot 2x2) ----------
grDevices::pdf(file.path(OUT_MAIN, "Fig5F_spearman_corrplot.pdf"),
               width = 178 / 25.4, height = 172 / 25.4,
               family = "Helvetica", useDingbats = FALSE)
op <- par(mfrow = c(2, 2), mar = c(2, 1.5, 2.5, 1))
col_div <- colorRampPalette(c("#3C5488FF", "white", "#E64B35FF"))(200)
for (z in DS_META) {
  cor22 <- as.matrix(read.csv(file.path(CBD, paste0(z$id, "_spearman_cor.csv")),
                              row.names = 1, check.names = FALSE))
  dimnames(cor22) <- list(CELL22, CELL22)
  cor22[is.na(cor22)] <- 0  # zero-variance cell types: undefined rho shown as 0 (white)
  corrplot(cor22, method = "color", col = col_div, is.corr = TRUE,
           tl.col = "black", tl.cex = 0.42, tl.srt = 45,
           cl.cex = 0.42, cl.length = 3, mar = c(0, 0, 2, 0))
  mtext(paste0(z$id, " (", z$tag, ")"), side = 3, line = 1, cex = 0.62, font = 2)
}
par(op); dev.off()

# ---------- G: differential abundance violins ----------
allg <- do.call(rbind, lapply(DS_META, function(z) {
  res <- read.csv(file.path(CBD, paste0(z$id, "_CIBERSORT_results.csv")), check.names = FALSE)
  dif <- read.csv(file.path(CBD, paste0(z$id, "_celltype_wilcox.csv")))
  long <- pivot_longer(res, cols = all_of(CELL22), names_to = "cell", values_to = "frac")
  long <- long %>% mutate(cell = factor(cell, levels = CELL22),
                          group = factor(group, levels = c("Control", z$dis)))
  stars <- dif %>% transmute(cell, star = sig_stars(P)) %>% filter(star != "ns")
  long$star <- stars$star[match(long$cell, stars$cell)]
  long$dataset <- paste0(z$id, " ", z$dis, " ",
                         ifelse(z$id %in% c("GSE76293", "GSE69528"), "discovery", "validation"))
  long
}))
allg$dataset <- factor(allg$dataset, levels = vapply(DS_META, function(z)
  paste0(z$id, " ", z$dis, " ",
         ifelse(z$id %in% c("GSE76293", "GSE69528"), "discovery", "validation")), character(1)))
star_lab <- allg %>% filter(!is.na(star)) %>%
  group_by(dataset, cell, star) %>% summarise(y = max(frac, na.rm = TRUE) + 0.06, .groups = "drop")
p5g <- ggplot(allg, aes(x = cell, y = frac, fill = group)) +
  geom_violin(scale = "width", linewidth = 0.22, color = "grey30", width = 0.85,
              position = position_dodge(0.85)) +
  geom_text(data = star_lab, aes(x = cell, y = y, label = star), inherit.aes = FALSE,
            size = 2.1, vjust = 0.3) +
  facet_wrap(~dataset, ncol = 2) +
  scale_x_discrete(breaks = CELL22, labels = CELL_SHORT[CELL22]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  scale_fill_manual(values = c(Control = COL_CONTROL, ARDS = COL_DISEASE, Sepsis = COL_DISEASE),
                    name = NULL) +
  coord_flip() +
  labs(title = "Differential immune-cell abundance (CIBERSORT, two-sided Wilcoxon)",
       x = NULL, y = "Estimated cell fraction") +
  theme_ajrccm(base_size = 6) +
  theme(axis.text.y = element_text(size = 5.4),
        strip.text = element_text(size = 6),
        legend.position = "top",
        panel.spacing = unit(3, "mm"))
save_pdf(p5g, file.path(OUT_MAIN, "Fig5G_celltype_violins.pdf"), 178, 150)

cat("step42: Figure 5 v2 panels done.\n")
