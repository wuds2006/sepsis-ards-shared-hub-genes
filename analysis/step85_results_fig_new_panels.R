#!/usr/bin/env Rscript
# Step 85 —— RESULTS 图组重设计：8 个新面板（R ggplot2 + ggsci NPG）
# 对应新 Figure 3C / 4E / 6A / 6B / 7A-D / 10A
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(ggsci); library(patchwork)
})
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
PROJ <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(PROJ, "AJRCCM_submission/03_figures_RESULTS_NPG/panels")

# ================= Fig 3C: 跨疾病候选验证（发现队列内交叉测试） =================
cd <- read.csv(file.path(PROJ, "05_ML_SHAP/discovery_crossdisease_AUC.csv"))
cd$list <- trimws(cd$list)
cd$status <- with(cd, case_when(
  gene %in% HUB6 ~ "Confirmed hub gene",
  AUC_GSE69528_SepsisDiscovery < 0.5 | AUC_GSE76293_ARDSdiscovery < 0.5 ~ "Direction-reversed (excluded)",
  TRUE ~ "Same-direction (not retained)"))
cd_ar <- cd[grepl("ARDS", cd$list), ]
cd_se <- cd[grepl("Sepsis", cd$list) & cd$gene != "TCTN1", ]
mk_cross <- function(df, auc_col, xlab_txt, title_txt) {
  df$auc <- df[[auc_col]]
  df <- df %>% arrange(desc(auc))
  df$gene <- factor(df$gene, levels = rev(df$gene))
  ggplot(df, aes(x = auc, y = gene, color = status)) +
    geom_segment(aes(x = 0.5, xend = auc, yend = gene), linewidth = 0.4,
                 show.legend = FALSE) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40", linewidth = 0.3) +
    geom_point(size = 1.9, aes(shape = status)) +
    scale_color_manual(values = c("Confirmed hub gene" = "#E64B35FF",
                                  "Same-direction (not retained)" = "grey55",
                                  "Direction-reversed (excluded)" = "#3C5488FF")) +
    scale_shape_manual(values = c("Confirmed hub gene" = 16,
                                  "Same-direction (not retained)" = 1,
                                  "Direction-reversed (excluded)" = 17)) +
    scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
    labs(x = xlab_txt, y = NULL, title = title_txt, color = NULL, shape = NULL) +
    theme_ajrccm(6.8) + theme(plot.title = element_text(size = 7, face = "bold", hjust = 0))
}
p_ar <- mk_cross(cd_ar, "AUC_GSE69528_SepsisDiscovery",
                 "Single-gene AUC in sepsis discovery cohort",
                 "ARDS-nominated candidates tested in sepsis discovery")
p_se <- mk_cross(cd_se, "AUC_GSE76293_ARDSdiscovery",
                 "Single-gene AUC in ARDS discovery cohort",
                 "Sepsis-nominated candidates tested in ARDS discovery")
p3c <- p_ar + p_se + patchwork::plot_layout(ncol = 2, widths = c(1.15, 1)) +
  patchwork::plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1)))
save_pdf(p3c, file.path(OUT, "R3C_crossdisease_AUC.pdf"), 178, 74)

# ================= Fig 4E: 单基因 GSEA NES 热图（4 队列 × 6 基因 × 8 通路） =================
gw <- read.csv(file.path(PROJ, "12_single_gene_gsea/single_gene_gsea_all_results_4cohorts.csv"))
path_sel <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
              "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
              "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_COMPLEMENT",
              "HALLMARK_MYC_TARGETS_V1", "HALLMARK_DNA_REPAIR")
lab <- c("TNF/\nNF-\u03baB", "IL6/JAK/\nSTAT3", "Inflam.", "IFN-\u03b1", "IFN-\u03b3",
         "Comple.", "MYC\ntargets", "DNA\nrepair")
coh_lab <- c("GSE76293" = "ARDS discovery", "GSE243066" = "ARDS validation",
             "GSE69528" = "Sepsis discovery", "GSE100159" = "Sepsis validation")
gw <- gw %>% filter(pathway %in% path_sel, hub_gene %in% HUB6) %>%
  mutate(pathway = factor(pathway, levels = path_sel),
         cohort = factor(coh_lab[cohort], levels = unname(coh_lab)),
         gene = factor(hub_gene, levels = rev(HUB6)),
         star = ifelse(padj < 0.05, "*", ""))
gw$pathway <- factor(gw$pathway, levels = path_sel, labels = lab)
p4e <- ggplot(gw, aes(pathway, gene, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = star), size = 3.2, vjust = 0.65) +
  scale_fill_gradient2(low = "#3C5488FF", mid = "white", high = "#E64B35FF",
                       midpoint = 0, limits = c(-3.2, 3.2), breaks = c(-3, 0, 3)) +
  facet_wrap(~cohort, nrow = 1) +
  labs(x = NULL, y = NULL, fill = "NES") +
  theme_ajrccm(6.5) +
  theme(axis.text.x = element_text(size = 5.6),
        strip.text = element_text(size = 6.2, face = "bold"),
        legend.key.height = unit(6, "mm"), legend.key.width = unit(2.2, "mm"),
        panel.border = element_blank(), panel.grid = element_blank())
save_pdf(p4e, file.path(OUT, "R4E_gsea_NES_heatmap.pdf"), 178, 62)

# ================= Fig 6A: 外部验证单基因 AUC（DeLong 95% CI） =================
va <- read.csv(file.path(PROJ, "06_external_validation/reviewer_w5/single_gene_AUC_CI_validation.csv"))
va$cohort <- gsub(" \\(.*\\)", "", va$cohort)
va$coh_lab <- paste0(va$cohort, c(GSE243066 = " (ARDS validation)", GSE100159 = " (Sepsis validation)")[va$cohort])
ord <- va %>% group_by(gene) %>% summarise(m = mean(AUC)) %>% arrange(m) %>% pull(gene)
va$gene <- factor(va$gene, levels = ord)
va$coh_lab <- factor(va$coh_lab, levels = c("GSE243066 (ARDS validation)", "GSE100159 (Sepsis validation)"))
p6a <- ggplot(va, aes(AUC, gene, color = coh_lab)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.5,
                 position = position_dodge(width = 0.55), show.legend = FALSE) +
  geom_point(size = 1.8, position = position_dodge(width = 0.55)) +
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF")) +
  scale_x_continuous(limits = c(0.45, 1.02), breaks = seq(0.5, 1.0, 0.1)) +
  labs(x = "AUC (disease vs. control, DeLong 95% CI)", y = NULL, color = NULL) +
  theme_ajrccm(6.8) + theme(legend.position = "top",
                            axis.text.y = element_text(face = "italic"))
save_pdf(p6a, file.path(OUT, "R6A_validation_AUC.pdf"), 120, 56)

# ================= Fig 6B: 跨队列 Meta 分析森林图（Hedges g） =================
md <- read.csv(file.path(PROJ, "06_external_validation/reviewer_w5/metaDE_single_gene_validation.csv"))
ci <- strcapture("([0-9.]+)-([0-9.]+)", md$pooled_g_FE_ci,
                 data.frame(lo = numeric(), hi = numeric()))
md$lo <- ci$lo; md$hi <- ci$hi
md$gene <- factor(md$gene, levels = rev(HUB6))
rows <- bind_rows(
  data.frame(gene = rep(md$gene, 2),
             g = c(md$g_GSE243066, md$g_GSE100159),
             lo = NA, hi = NA,
             what = rep(c("GSE243066 (ARDS validation)", "GSE100159 (Sepsis validation)"), each = nrow(md))),
  data.frame(gene = md$gene, g = md$pooled_g_FE, lo = md$lo, hi = md$hi,
             what = "Pooled (fixed effect)"))
rows$what <- factor(rows$what, levels = c("GSE243066 (ARDS validation)",
                                           "GSE100159 (Sepsis validation)",
                                           "Pooled (fixed effect)"))
p6b <- ggplot(rows, aes(g, gene, color = what, shape = what)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_errorbarh(data = ~ subset(.x, what == "Pooled (fixed effect)"),
                 aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.55,
                 color = "grey25", show.legend = FALSE) +
  geom_point(size = 1.7, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "grey20")) +
  scale_shape_manual(values = c(16, 16, 18)) +
  labs(x = "Standardized effect (Hedges g)", y = NULL, color = NULL, shape = NULL) +
  theme_ajrccm(6.8) + theme(legend.position = "top",
                            axis.text.y = element_text(face = "italic"))
save_pdf(p6b, file.path(OUT, "R6B_meta_forest.pdf"), 120, 58)

# ================= Fig 7A/B: MR 森林图（forward / reverse） =================
mk_mr_forest <- function(csv, title_txt) {
  mr <- read.csv(csv)
  mr$method <- c("Inverse variance weighted" = "IVW", "MR Egger" = "MR Egger",
                 "Weighted median" = "Weighted median", "Weighted mode" = "Weighted mode")[mr$method]
  mr <- mr %>% mutate(or = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),
                      P = sprintf("P = %.3f", pval)) %>%
    mutate(method = factor(method, levels = rev(c("IVW", "MR Egger", "Weighted median", "Weighted mode"))))
  ggplot(mr, aes(or, method)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.5, color = "grey30") +
    geom_point(aes(shape = method == "IVW", size = method == "IVW"), color = "#E64B35FF",
               show.legend = FALSE) +
    geom_text(aes(x = hi, label = paste0(sprintf("%.2f", or), " (", sprintf("%.2f", lo), "-", sprintf("%.2f", hi), ")  "),
                  hjust = 0), size = 1.9, family = AJRCCM_FONT, vjust = -0.9) +
    geom_text(aes(x = 4.6, label = P), size = 1.9, family = AJRCCM_FONT) +
    scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 16), guide = "none") +
    scale_size_manual(values = c(`TRUE` = 2.4, `FALSE` = 1.7), guide = "none") +
    scale_x_log10(limits = c(0.15, 5.2), breaks = c(0.25, 0.5, 1, 2, 4)) +
    labs(x = "Odds ratio (log scale)", y = NULL, title = title_txt) +
    theme_ajrccm(6.8) + theme(plot.title = element_text(size = 7, face = "bold", hjust = 0))
}
p7a <- mk_mr_forest(file.path(PROJ, "04_MR_SMR/MR_sepsis4980_ARDS_results.csv"),
                    "Sepsis \u2192 ARDS (forward, 18 instruments)")
p7b <- mk_mr_forest(file.path(PROJ, "04_MR_SMR/MR_reverse_ards_sepsis_results.csv"),
                    "ARDS \u2192 Sepsis (reverse, 7 instruments)")
save_pdf(p7a, file.path(OUT, "R7A_mr_forest_forward.pdf"), 120, 44)
save_pdf(p7b, file.path(OUT, "R7B_mr_forest_reverse.pdf"), 120, 44)

# ================= Fig 7C: MR 检验效能曲线 =================
pw <- read.csv(file.path(PROJ, "04_MR_SMR/MR_power_by_OR.csv"))
pw$dir_lab <- c("forward_sepsis_to_ARDS" = "Forward: sepsis \u2192 ARDS (n = 18 IVs)",
                "reverse_ARDS_to_sepsis" = "Reverse: ARDS \u2192 sepsis (n = 7 IVs)")[pw$direction]
p7c <- ggplot(pw, aes(OR, power_fixed, color = dir_lab)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
  annotate("text", x = 1.42, y = 0.83, label = "80% power", size = 2, family = AJRCCM_FONT,
           color = "grey30") +
  scale_color_manual(values = c("#E64B35FF", "#3C5488FF")) +
  scale_x_continuous(breaks = c(1.1, 1.2, 1.5, 2, 2.5, 3)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(x = "Causal odds ratio", y = "Power", color = NULL) +
  theme_ajrccm(6.8) + theme(legend.position = "top")
save_pdf(p7c, file.path(OUT, "R7C_mr_power.pdf"), 100, 62)

# ================= Fig 7D: SMR-HEIDI（5 基因 × 2 结局） =================
smr_files <- list.files(file.path(PROJ, "04_MR_SMR"), pattern = "smr_.*_heidi\\.smr$", full.names = TRUE)
smr <- do.call(rbind, lapply(smr_files, function(f) {
  x <- read.delim(f)
  gene <- sub("smr_(.*)_.*_heidi\\.smr$", "\\1", basename(f))
  outcome <- ifelse(grepl("_sepsis_", basename(f)), "Sepsis", "ARDS")
  data.frame(gene = toupper(gene), outcome, b_SMR = x$b_SMR, se_SMR = x$se_SMR,
             p_SMR = x$p_SMR, p_HEIDI = x$p_HEIDI)
}))
smr$gene <- factor(smr$gene, levels = rev(HUB6[HUB6 != "EXOSC4"]))
smr$outcome <- factor(smr$outcome, levels = c("Sepsis", "ARDS"))
p7d <- ggplot(smr, aes(b_SMR, gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = b_SMR - 1.96 * se_SMR, xmax = b_SMR + 1.96 * se_SMR,
                     color = outcome), height = 0, linewidth = 0.5,
                 position = position_dodge(width = 0.55), show.legend = FALSE) +
  geom_point(aes(color = outcome, shape = outcome), size = 1.8,
             position = position_dodge(width = 0.55)) +
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF")) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "SMR causal estimate (log odds per SD expression)", y = NULL,
       color = NULL, shape = NULL,
       subtitle = "All P_SMR > 0.1; P_HEIDI > 0.05 for all 10 tests; no colocalization (PP.H4 < 0.05)") +
  theme_ajrccm(6.8) + theme(legend.position = "top",
                            axis.text.y = element_text(face = "italic"),
                            plot.subtitle = element_text(size = 5.6, color = "grey25"))
save_pdf(p7d, file.path(OUT, "R7D_smr_forest.pdf"), 120, 58)

# ================= Fig 10A: HPA Blood Atlas 蛋白水平验证 =================
hpa <- read.csv(file.path(PROJ, "10_HPA_validation/HPA_CD177_VNN2_infection_vs_healthy.csv"))
hpa$Disease <- gsub("Bacterial gastroenteritis", "Bact. gastroenteritis", hpa$Disease)
hpa$Disease <- gsub("Mycoplasma pneumonia", "M. pneumonia", hpa$Disease)
ord_h <- hpa %>% filter(Gene == "CD177") %>% arrange(logFC) %>% pull(Disease)
hpa$Disease <- factor(hpa$Disease, levels = ord_h)
p10a <- ggplot(hpa, aes(logFC, Disease, color = Gene, shape = Gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_point(aes(size = neglogp), position = position_dodge(width = 0.6)) +
  scale_color_manual(values = c(CD177 = "#E64B35FF", VNN2 = "#3C5488FF")) +
  scale_shape_manual(values = c(CD177 = 16, VNN2 = 17)) +
  scale_size_continuous(range = c(1.2, 3.2), name = expression(-log[10]~italic(P))) +
  labs(x = "Protein log2 fold change vs. healthy\n(Human Protein Atlas Blood Atlas, n = 10 infectious diseases)",
       y = NULL, color = NULL, shape = NULL) +
  theme_ajrccm(6.8) + theme(legend.position = "top",
                            legend.box = "horizontal", legend.direction = "horizontal",
                            axis.text.y = element_text(size = 5.8))
save_pdf(p10a, file.path(OUT, "R10A_hpa_infection.pdf"), 128, 78)

cat("All 8 new panels written to", OUT, "\n")
