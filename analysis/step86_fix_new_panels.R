#!/usr/bin/env Rscript
# Step 86 —— 修复 3 个新面板：R4E 标签重叠/NA、R6B 图例重复、R7D 多行 top-SNP
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(ggsci)
})
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
PROJ <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(PROJ, "AJRCCM_submission/03_figures_RESULTS_NPG/panels")

# ---- R4E: GSEA NES heatmap（标签旋转 + NA 明示 + star 条件化） ----
gw <- read.csv(file.path(PROJ, "12_single_gene_gsea/single_gene_gsea_all_results_4cohorts.csv"))
path_sel <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
              "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
              "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_COMPLEMENT",
              "HALLMARK_MYC_TARGETS_V1", "HALLMARK_DNA_REPAIR")
coh_lab <- c("GSE76293" = "ARDS discovery", "GSE243066" = "ARDS validation",
             "GSE69528" = "Sepsis discovery", "GSE100159" = "Sepsis validation")
gw <- gw %>% filter(pathway %in% path_sel, hub_gene %in% HUB6) %>%
  mutate(pathway = factor(pathway, levels = path_sel),
         cohort = factor(coh_lab[cohort], levels = unname(coh_lab)),
         gene = factor(hub_gene, levels = rev(HUB6)),
         star = ifelse(!is.na(NES) & !is.na(padj) & padj < 0.05, "*", ""))
pl_lab <- parse(text = c(
  "paste('TNF/',NF*kappa*B)", "'IL6/STAT3'", "'Inflamm.'", "IFN*'-'*alpha",
  "IFN*'-'*gamma", "'Compl.'", "'MYC targets'", "'DNA repair'"))
p4e <- ggplot(gw, aes(pathway, gene, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = star), size = 3.2, vjust = 0.65) +
  scale_fill_gradient2(low = "#3C5488FF", mid = "white", high = "#E64B35FF",
                       midpoint = 0, limits = c(-3.2, 3.2), breaks = c(-3, 0, 3),
                       na.value = "grey92") +
  scale_x_discrete(labels = pl_lab) +
  facet_wrap(~cohort, nrow = 1) +
  labs(x = NULL, y = NULL, fill = "NES") +
  theme_ajrccm(6.5) +
  theme(axis.text.x = element_text(size = 4.6, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 6.2, face = "bold"),
        legend.key.height = unit(6, "mm"), legend.key.width = unit(2.2, "mm"),
        panel.border = element_blank(), panel.grid = element_blank())
save_pdf(p4e, file.path(OUT, "R4E_gsea_NES_heatmap.pdf"), 178, 74)

# ---- R6B: meta forest（单层映射 + 统一图例） ----
md <- read.csv(file.path(PROJ, "06_external_validation/reviewer_w5/metaDE_single_gene_validation.csv"))
ci <- strcapture("([0-9.]+)-([0-9.]+)", md$pooled_g_FE_ci,
                 data.frame(lo = numeric(), hi = numeric()))
md$lo <- ci$lo; md$hi <- ci$hi
md$gene <- factor(md$gene, levels = rev(HUB6))
per <- bind_rows(
  data.frame(gene = md$gene, g = md$g_GSE243066, what = "GSE243066 (ARDS validation)"),
  data.frame(gene = md$gene, g = md$g_GSE100159, what = "GSE100159 (Sepsis validation)"))
per$what <- factor(per$what, levels = c("GSE243066 (ARDS validation)", "GSE100159 (Sepsis validation)"))
pool <- data.frame(gene = md$gene, g = md$pooled_g_FE, lo = md$lo, hi = md$hi)
p6b <- ggplot(per, aes(g, gene, color = what)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_point(size = 1.6, shape = 16, position = position_dodge(width = 0.55)) +
  geom_errorbarh(data = pool, aes(x = g, y = gene, xmin = lo, xmax = hi),
                 inherit.aes = FALSE, height = 0, linewidth = 0.55, color = "grey25") +
  geom_point(data = pool, aes(x = g, y = gene), inherit.aes = FALSE,
             shape = 18, size = 2.6, color = "grey20") +
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF"),
                     guide = guide_legend(title = NULL, nrow = 1,
                                          override.aes = list(shape = 16, size = 1.8))) +
  scale_x_continuous(limits = c(-0.3, 3.9), breaks = 0:3) +
  annotate("text", x = 3.75, y = 6, label = "diamond = pooled (fixed effect)",
           size = 2.1, family = AJRCCM_FONT, hjust = 1, color = "grey20") +
  labs(x = "Standardized effect (Hedges g)", y = NULL) +
  theme_ajrccm(6.8) + theme(legend.position = "top",
                            axis.text.y = element_text(face = "italic"))
save_pdf(p6b, file.path(OUT, "R6B_meta_forest.pdf"), 120, 62)

# ---- R7D: SMR forest（按 manuscript/SMR_coloc_summary.txt 指定的 topSNP 提取） ----
snp_map <- c(VNN1 = "rs2300077", UPP1 = "rs7459020", GRB10 = "rs73113863",
             CDK5RAP2 = "rs7037806", CD177 = "rs34627620")
smr <- do.call(rbind, lapply(names(snp_map), function(g) {
  do.call(rbind, lapply(c("sepsis", "ards"), function(oc) {
    f <- file.path(PROJ, sprintf("04_MR_SMR/smr_%s_%s_heidi.smr", tolower(g), oc))
    x <- read.delim(f)
    r <- x[x$topSNP == snp_map[g], , drop = FALSE]
    stopifnot(nrow(r) == 1)
    data.frame(gene = g, outcome = ifelse(oc == "sepsis", "Sepsis", "ARDS"),
               b_SMR = r$b_SMR, se_SMR = r$se_SMR, p_SMR = r$p_SMR, p_HEIDI = r$p_HEIDI)
  }))
}))
smr$gene <- factor(smr$gene, levels = rev(HUB6[HUB6 != "EXOSC4"]))
smr$outcome <- factor(smr$outcome, levels = c("Sepsis", "ARDS"))
cat("SMR top-SNP rows:\n"); print(smr)
stopifnot(all(smr$p_SMR > 0.1), all(smr$p_HEIDI > 0.05))
p7d <- ggplot(smr, aes(b_SMR, gene, color = outcome, shape = outcome)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = b_SMR - 1.96 * se_SMR, xmax = b_SMR + 1.96 * se_SMR),
                 height = 0, linewidth = 0.5,
                 position = position_dodge(width = 0.55), show.legend = FALSE) +
  geom_point(size = 1.9, position = position_dodge(width = 0.55)) +
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF")) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "SMR causal estimate (log odds per SD expression)", y = NULL,
       color = NULL, shape = NULL,
       subtitle = "All P_SMR > 0.1; P_HEIDI > 0.05 for all 10 tests; no colocalization (PP.H4 < 0.05)") +
  theme_ajrccm(6.8) + theme(legend.position = "top",
                            axis.text.y = element_text(face = "italic"),
                            plot.subtitle = element_text(size = 5.4, color = "grey25"))
save_pdf(p7d, file.path(OUT, "R7D_smr_forest.pdf"), 120, 62)
cat("Step 86 done\n")
