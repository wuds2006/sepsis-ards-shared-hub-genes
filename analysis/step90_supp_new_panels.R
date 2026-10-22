# ============================================================
# Step 90 —— Supplementary_NPG 新面板（R ggplot2 + ggsci NPG）
# 输出: AJRCCM_submission/05_supplementary_NPG/panels/
# 规范: base pdf() 设备 (cairo 本机不可用), Helvetica, 正文 7pt
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2); library(ggsci); library(dplyr); library(tidyr)
  library(stringr); library(readr); library(tibble)
})

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "AJRCCM_submission/05_supplementary_NPG/panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

NPG  <- ggsci::pal_npg("nrc")(10)
RED  <- NPG[1]; CYAN <- NPG[2]; GREEN <- NPG[3]; BLUE <- NPG[4]
CORAL <- NPG[5]; GREY6 <- NPG[6]
FONT <- "Helvetica"

theme_supp <- theme_classic(base_size = 7, base_family = FONT) +
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        axis.text = element_text(color = "grey10"),
        axis.title = element_text(color = "grey10"),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        strip.background = element_blank(),
        strip.text = element_text(size = 7, face = "bold"),
        plot.margin = margin(3, 4, 3, 3))

save_pdf <- function(p, name, w, h) {
  pdf(file.path(OUT, paste0(name, ".pdf")), width = w, height = h,
      family = FONT, useDingbats = FALSE)
  print(p); dev.off()
  cat("saved", name, "\n")
}

ital <- function(g) parse(text = paste0("italic('", g, "')"))

# ------------------------------------------------------------
# S2: hub-gene co-expression evidence (kME + degree in 179-CG network)
# ------------------------------------------------------------
kme <- read_csv(file.path(BASE, "02_WGCNA/hub_gene_kME.csv"),
                show_col_types = FALSE)
gene_ord <- c("VNN1", "UPP1", "CD177", "CDK5RAP2", "GRB10", "EXOSC4")
kme$gene <- factor(kme$gene, levels = gene_ord)

p <- ggplot(kme, aes(x = gene, y = kME, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = NA) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey45",
             linewidth = 0.3) +
  annotate("text", x = 0.55, y = 0.82, label = "kME = 0.8", hjust = 0,
           size = 2, color = "grey45", family = FONT) +
  scale_fill_manual(values = c("GSE76293" = RED, "GSE69528" = CYAN),
                    labels = c("GSE76293 (ARDS)", "GSE69528 (sepsis)")) +
  scale_x_discrete(labels = gene_ord) +
  labs(x = NULL, y = "Module membership (kME)", fill = NULL) +
  theme_supp + theme(legend.position = c(0.98, 0.02),
                     legend.justification = c(1, 0),
                     legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
                     axis.text.x = element_text(face = "italic"))
save_pdf(p, "RS2A_hub_kME", 3.5, 2.7)

deg <- kme %>%
  select(cohort, gene, `r0.5` = `degree_r0.5_in_179CG`,
         `r0.7` = `degree_r0.7_in_179CG`) %>%
  pivot_longer(cols = c(`r0.5`, `r0.7`), names_to = "threshold",
               values_to = "degree")
deg$threshold <- factor(deg$threshold, levels = c("r0.5", "r0.7"),
                        labels = c("cor > 0.5", "cor > 0.7"))
p <- ggplot(deg, aes(x = gene, y = degree, fill = threshold)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = NA) +
  facet_wrap(~cohort, ncol = 2,
             labeller = labeller(cohort = c("GSE76293" = "GSE76293 (ARDS)",
                                            "GSE69528" = "GSE69528 (sepsis)"))) +
  scale_fill_manual(values = c("cor > 0.5" = BLUE, "cor > 0.7" = CORAL)) +
  labs(x = NULL, y = "Co-expression degree (of 178)", fill = NULL) +
  theme_supp + theme(legend.position = c(0.98, 0.02),
                     legend.justification = c(1, 0),
                     legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
                     axis.text.x = element_text(face = "italic", angle = 30,
                                                hjust = 1))
save_pdf(p, "RS2B_hub_degree_179CG", 3.5, 2.7)

# ------------------------------------------------------------
# S3: functional enrichment supplements
# ------------------------------------------------------------
wrap20 <- function(x) str_wrap(x, width = 42)

kegg <- read_csv(file.path(BASE, "03_intersection/Shared_CGs_KEGG.csv"),
                 show_col_types = FALSE) %>%
  arrange(p.adjust) %>% slice_head(n = 15) %>%
  mutate(Description = factor(wrap20(Description),
                              levels = rev(wrap20(Description))))
p <- ggplot(kegg, aes(x = Count, y = Description, fill = -log10(p.adjust))) +
  geom_col(width = 0.65, color = NA) +
  scale_fill_gradient(low = CYAN, high = RED, name = expression(-log[10](P[adj]))) +
  labs(x = "Gene count", y = NULL) +
  theme_supp + theme(legend.position = c(0.98, 0.02),
                     legend.justification = c(1, 0),
                     legend.background = element_rect(fill = alpha("white", 0.8), color = NA))
save_pdf(p, "RS3A_KEGG_179_bar", 3.6, 3.4)

enrichr_dot <- function(path, name, top = 12) {
  d <- read_csv(path, show_col_types = FALSE) %>%
    mutate(combined_score = suppressWarnings(as.numeric(combined_score)),
           pvalue = as.numeric(pvalue),
           Count = lengths(str_split(overlapping_genes, ","))) %>%
    filter(!is.na(combined_score), combined_score > 0) %>%
    arrange(pvalue) %>% slice_head(n = top) %>%
    mutate(term = str_remove(term, " R-HSA-\\d+$"),
           term = str_remove(term, " WP\\d+$"),
           term = str_remove(term, " \\(GO:\\d+\\)$"),
           term = factor(wrap20(term), levels = rev(wrap20(term))))
  ggplot(d, aes(x = log10(combined_score), y = term)) +
    geom_point(aes(size = Count, color = -log10(pvalue))) +
    scale_color_gradient(low = BLUE, high = RED,
                         name = expression(-log[10](P))) +
    scale_size_continuous(range = c(1.2, 3.5), breaks = pretty(d$Count, 3)) +
    labs(x = expression(log[10]("combined score")), y = NULL) +
    theme_supp + theme(legend.position = "right",
                       legend.box = "vertical",
                       legend.spacing.y = unit(1, "mm"))
}
save_pdf(enrichr_dot(file.path(BASE, "11_network_pharma/Reactome_2022_enrichment.csv")),
         "RS3B_Reactome_hub_dot", 3.7, 3.0)
# WikiPathways_2023 富集无显著条目（源 csv 仅表头）-> 以 KEGG 2021 (hub) 替代
save_pdf(enrichr_dot(file.path(BASE, "11_network_pharma/KEGG_2021_Human_enrichment.csv")),
         "RS3C_KEGG2021_hub_dot", 3.7, 3.0)
save_pdf(enrichr_dot(file.path(BASE, "11_network_pharma/GO_Biological_Process_2023_enrichment.csv")),
         "RS3D_GOBP2023_hub_dot", 3.7, 3.0)

# ------------------------------------------------------------
# S5: ML robustness (apparent vs CV, LASSO bootstrap, Jaccard, permutation)
# ------------------------------------------------------------
appcv <- bind_rows(
  read_csv(file.path(BASE, "05_ML_SHAP/stability/apparent_vs_cv_AUC_ARDS.csv"),
           show_col_types = FALSE) %>% mutate(cohort = "GSE76293 (ARDS)"),
  read_csv(file.path(BASE, "05_ML_SHAP/stability/apparent_vs_cv_AUC_Sepsis.csv"),
           show_col_types = FALSE) %>% mutate(cohort = "GSE69528 (sepsis)"))
for (cc in unique(appcv$cohort)) {
  d <- filter(appcv, cohort == cc) %>%
    select(model, `5-fold CV` = cv_AUC, Apparent = apparent_AUC) %>%
    pivot_longer(-model, names_to = "estimate", values_to = "AUC")
  p <- ggplot(d, aes(x = reorder(model, -AUC), y = AUC, fill = estimate)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, color = NA) +
    coord_cartesian(ylim = c(0.5, 1.02)) +
    scale_fill_manual(values = c("5-fold CV" = GREEN, "Apparent" = GREY6)) +
    labs(x = NULL, y = "AUC", fill = NULL, title = cc) +
    theme_supp + theme(legend.position = c(0.98, 0.98),
                       legend.justification = c(1, 1),
                       legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
                       plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
                       axis.text.x = element_text(angle = 35, hjust = 1))
  save_pdf(p, paste0("RS5", if (grepl("ARDS", cc)) "A" else "B",
                     "_apparent_vs_cv_",
                     if (grepl("ARDS", cc)) "GSE76293" else "GSE69528"),
           3.5, 2.8)
}

for (cc in c("ARDS", "Sepsis")) {
  d <- read_csv(file.path(BASE, "05_ML_SHAP/stability",
                          paste0("lasso_bootstrap_stability_", cc, ".csv")),
                show_col_types = FALSE) %>%
    arrange(desc(boot_selection_freq)) %>% slice_head(n = 15) %>%
    mutate(is_hub = gene %in% gene_ord,
           gene = factor(gene, levels = rev(gene)))
  p <- ggplot(d, aes(x = boot_selection_freq, y = gene, fill = is_hub)) +
    geom_col(width = 0.65, color = NA) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.08))) +
    scale_fill_manual(values = c("TRUE" = RED, "FALSE" = "grey70"),
                      labels = c("TRUE" = "hub gene", "FALSE" = "other"),
                      guide = guide_legend(reverse = TRUE)) +
    labs(x = "Bootstrap selection frequency", y = NULL, fill = NULL,
         title = if (cc == "ARDS") "GSE76293 (ARDS)" else "GSE69528 (sepsis)") +
    theme_supp + theme(legend.position = c(0.98, 0.02),
                       legend.justification = c(1, 0),
                       legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
                       plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
                       axis.text.y = element_text(face = "italic"))
  save_pdf(p, paste0("RS5", if (cc == "ARDS") "C" else "D",
                     "_lasso_bootstrap_", cc), 3.5, 3.2)
}

for (cc in c("ARDS", "Sepsis")) {
  d <- read_csv(file.path(BASE, "05_ML_SHAP/stability",
                          paste0("threshold_perturbation_jaccard_", cc, ".csv")),
                show_col_types = FALSE)
  p <- ggplot(d, aes(x = factor(rf_k), y = factor(shap_m),
                     fill = jaccard_vs_base)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", jaccard_vs_base)),
              size = 2.1, family = FONT, color = "grey10") +
    scale_fill_gradient(low = tint <- "#F4F6FB", high = BLUE,
                        limits = c(0, 1), name = "Jaccard\nvs base") +
    labs(x = "RF top-k", y = "SHAP top-m",
         title = if (cc == "ARDS") "GSE76293 (ARDS)" else "GSE69528 (sepsis)") +
    theme_supp + theme(plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
                       panel.border = element_blank(),
                       axis.line = element_blank(),
                       axis.ticks = element_blank())
  save_pdf(p, paste0("RS5", if (cc == "ARDS") "E" else "F",
                     "_jaccard_", cc), 3.3, 2.8)
}

perm <- read_csv(file.path(BASE, "05_ML_SHAP/stability/permutation_500x_results.csv"),
                 show_col_types = FALSE)
nd <- readRDS(file.path(BASE, "05_ML_SHAP/stability/permutation_null_dist.rds"))
null_vec <- if (is.list(nd) && !is.data.frame(nd)) {
  nm <- names(nd); pick <- if (!is.null(nm) && "XGBoost" %in% nm) "XGBoost" else nm[1]
  nd[[pick]]
} else if (is.data.frame(nd)) {
  nd[[1]]
} else nd
best <- perm %>% arrange(desc(real_cv_AUC)) %>% slice(1)
pd <- tibble(null_auc = as.numeric(null_vec))
p <- ggplot(pd, aes(x = null_auc)) +
  geom_histogram(bins = 40, fill = GREY6, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = best$real_cv_AUC, color = RED, linewidth = 0.7,
             linetype = "dashed") +
  annotate("text", x = best$real_cv_AUC, y = Inf,
           label = paste0("observed CV AUC = ", sprintf("%.2f", best$real_cv_AUC),
                          "\npermutation P = ", format(best$perm_P, digits = 2)),
           vjust = 1.4, hjust = 1.05, size = 2.3, family = FONT, color = RED) +
  labs(x = "Permuted CV AUC (500 label permutations)",
       y = "Frequency") +
  theme_supp
save_pdf(p, "RS5G_permutation_null", 3.5, 2.6)

# ------------------------------------------------------------
# S7: external-validation diagnostics
# ------------------------------------------------------------
muc <- read_csv(file.path(BASE, "06_external_validation/muc1_sensitivity.csv"),
                show_col_types = FALSE)
p <- ggplot(muc, aes(x = str_wrap(model, 12), y = valid_auc)) +
  geom_col(width = 0.55, fill = GREEN, color = NA) +
  geom_errorbar(aes(ymin = valid_lo, ymax = valid_hi), width = 0.18,
                linewidth = 0.35, color = "grey20") +
  coord_cartesian(ylim = c(0.8, 1.01)) +
  labs(x = NULL, y = "Validation AUC (95% CI)",
       title = "GSE243066 (ARDS validation)") +
  theme_supp + theme(plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
                     axis.text.x = element_text(size = 5.5, lineheight = 0.9))
save_pdf(p, "RS7B_MUC1_sensitivity", 3.5, 2.8)

a1 <- read_csv(file.path(BASE, "06_external_validation/AUC1_exact_lower_bounds.csv"),
               show_col_types = FALSE) %>%
  mutate(cohort = str_wrap(cohort, 26),
         cohort = factor(cohort, levels = rev(cohort)))
p <- ggplot(a1, aes(x = exact_lower_bound_AUC, y = cohort)) +
  geom_col(width = 0.5, fill = BLUE, color = NA) +
  geom_text(aes(label = sprintf("%.4f  (n = %d vs %d)",
                                exact_lower_bound_AUC, n_case, n_ctrl)),
            hjust = 1.05, size = 2.2, family = FONT, color = "white") +
  coord_cartesian(xlim = c(0.9, 1.0)) +
  labs(x = "Exact one-sided lower 95% bound of AUC = 1 estimate", y = NULL) +
  theme_supp
save_pdf(p, "RS7C_AUC1_exact_lower_bounds", 3.5, 2.2)

sp <- read_csv(file.path(BASE, "06_external_validation/GSE243066_AUC_sepsis_candidates.csv"),
               show_col_types = FALSE) %>%
  arrange(desc(AUC)) %>% slice_head(n = 12) %>%
  mutate(gene = factor(gene, levels = rev(gene)))
p <- ggplot(sp, aes(x = AUC, y = gene)) +
  geom_col(width = 0.65, fill = CYAN, color = NA) +
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "grey45",
             linewidth = 0.3) +
  coord_cartesian(xlim = c(0.5, 1.0)) +
  labs(x = "AUC in GSE243066 (ARDS cohort)", y = NULL,
       title = "Sepsis-ML candidate genes scored in ARDS cohort") +
  theme_supp + theme(plot.title = element_text(size = 6.5, face = "bold", hjust = 0.5),
                     axis.text.y = element_text(face = "italic"))
save_pdf(p, "RS7D_sepsis_candidates_in_ARDS", 3.5, 3.0)

# ------------------------------------------------------------
# S14: AlphaFold pLDDT confidence
# ------------------------------------------------------------
pl <- read_csv(file.path(BASE, "09_docking_drug/alphafold_plddt_assessment.csv"),
               show_col_types = FALSE) %>%
  mutate(gene = factor(gene, levels = gene_ord)) %>%
  select(gene, `Whole structure` = mean_pLDDT, `Docking box` = box_mean_pLDDT) %>%
  pivot_longer(-gene, names_to = "region", values_to = "pLDDT")
p <- ggplot(pl, aes(x = gene, y = pLDDT, fill = region)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = NA) +
  geom_hline(yintercept = 70, linetype = "dashed", color = "grey45",
             linewidth = 0.3) +
  annotate("text", x = 0.55, y = 72, label = "pLDDT = 70 (confident)",
           hjust = 0, size = 2, color = "grey45", family = FONT) +
  scale_fill_manual(values = c("Whole structure" = GREY6, "Docking box" = GREEN)) +
  labs(x = NULL, y = "Mean AlphaFold pLDDT", fill = NULL) +
  theme_supp + theme(legend.position = c(0.98, 0.98),
                     legend.justification = c(1, 1),
                     legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
                     axis.text.x = element_text(face = "italic"))
save_pdf(p, "RS14D_alphafold_pLDDT", 3.5, 2.7)

cat("step90 done ->", OUT, "\n")
