# ============================================================
# Step 60 —— 供体级数值补充图 Figure S17（审稿 W4 / 详细评语 C2）
# A: 每供体 CellOracle CEBPB 敲除 → UPP1 delta（5 供体）
# B: 每供体 UPP1-拟时序 Spearman rho（GSE167363, 5 供体）
# C: 供体级 pseudobulk DESeq2 log2FC 热图（基因 x 细胞型, GSE167363）
# 产出: supp_panels/FigS17_donor_level.pdf → Figures_final/FigureS17_composite.*
# ============================================================
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
suppressPackageStartupMessages({ library(patchwork); library(dplyr) })

SC <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/08_singlecell"

## A: per-donor CellOracle delta
co <- read.csv(file.path(SC, "celloracle_per_donor", "CellOracle_CEBPB_KO_per_donor.csv"))
co$donor <- sub("^GSM[0-9]+_", "", co$donor)
co$group <- ifelse(grepl("^HC", co$donor), "Healthy", "Sepsis")
co$donor <- factor(co$donor, levels = c("HC1", "HC2", "S2_T0", "S2_T6", "S3_T0"))
pA <- ggplot(co, aes(donor, delta, fill = group)) +
  geom_col(width = 0.62) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(Healthy = COL_CONTROL, Sepsis = COL_DISEASE)) +
  labs(title = "A  Per-donor CellOracle: CEBPB knockout effect on UPP1",
       subtitle = "independent GRN per donor; all five donors negative",
       x = NULL, y = expression(Delta*" simulated UPP1 (per donor)")) +
  theme_ajrccm() + theme(legend.position = "bottom", legend.title = element_blank(),
                         plot.subtitle = element_text(size = 5.5))

## B: per-donor UPP1-pseudotime rho
dl <- read.csv(file.path(SC, "donor_level", "GSE167363_donor_level_UPP1_CEBPB.csv"))
dl$group <- ifelse(dl$group == "HC", "Healthy", "Sepsis")
dl$lab_y <- dl$UPP1_pseudotime_rho + ifelse(dl$UPP1_pseudotime_rho > 0, 0.05, -0.05)
pB <- ggplot(dl, aes(donor, UPP1_pseudotime_rho, fill = group)) +
  geom_col(width = 0.62) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_text(aes(label = sprintf("rho = %.2f", UPP1_pseudotime_rho), y = lab_y),
            size = 1.9) +
  scale_fill_manual(values = c(Healthy = COL_CONTROL, Sepsis = COL_DISEASE)) +
  scale_y_continuous(limits = c(-0.55, 0.55)) +
  labs(title = "B  Per-donor UPP1-pseudotime association",
       subtitle = "Spearman rho per donor; cell-level P < 0.001 for each",
       x = NULL, y = "Spearman rho") +
  theme_ajrccm() + theme(legend.position = "none", plot.subtitle = element_text(size = 5.5))

## C: pseudobulk DESeq2 log2FC heatmap（GSE167363 Sepsis vs HC）
pb <- read.csv(file.path(SC, "donor_level", "pseudobulk_DESeq2_GSE167363_Sepsis_vs_HC.csv"))
pb$gene <- factor(pb$gene, levels = HUB6)
pb$sig <- cut(pb$padj, c(-Inf, 0.001, 0.01, 0.05, Inf),
              labels = c("***", "**", "*", ""))
pC <- ggplot(pb[!is.na(pb$gene), ], aes(celltype, gene, fill = log2FoldChange)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = ifelse(is.na(log2FoldChange), "NA",
                              ifelse(padj < 0.05 & !is.na(padj), "*", ""))),
            size = 3.4, vjust = 0.75) +
  scale_fill_gradient2(low = "#3C5488FF", mid = "white", high = "#E64B35FF",
                       midpoint = 0, name = "log2 FC", na.value = "grey90") +
  scale_y_discrete(labels = function(x) parse(text = gene_italic(x))) +
  labs(title = "C  Donor-level pseudobulk DESeq2 (GSE167363, sepsis vs control)",
       subtitle = "raw counts aggregated per donor x cell type; 2 vs 2 donors; * P < 0.05",
       x = NULL, y = NULL) +
  theme_ajrccm() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6),
        plot.subtitle = element_text(size = 5.5))

p <- pA + pB + pC + plot_layout(nrow = 3, heights = c(1, 1, 1.25)) +
  plot_annotation(title = "Donor-level sensitivity analyses of single-cell findings",
                  theme = theme(plot.title = element_text(size = 8, face = "bold")))
save_pdf(p, file.path(OUT_SUPP, "FigS17_donor_level.pdf"), 140, 210)
print(p)
cat("\n=== Step 60 完成 ===\n")
