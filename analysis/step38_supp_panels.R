#!/usr/bin/env Rscript
# Supplementary figure panels S1-S14 (ggplot2 + ggsci NPG, Illustrator-editable PDF)
# Specialized renders kept from the original pipeline: S5 (CellChat), S8 (running-ES
# GSEA curves), S11 (monocle3 trajectory graphs), S1 dendrograms.
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(pROC) })

B <- PROJ
SC <- file.path(B, "08_singlecell")

safe <- function(tag, expr) {
  tryCatch(expr, error = function(e) cat("SKIP", tag, ":", conditionMessage(e), "\n"))
}

## ---- S1: WGCNA soft-threshold curves (recomputed) ----
safe("S1", {
  library(WGCNA, quietly = TRUE)
  sft_plot <- function(ds, letter) {
    expr <- readRDS(file.path(B, "01_DEG", paste0(ds, "_gene_expr.rds")))
    dat <- t(as.matrix(expr))
    powers <- c(1:10, seq(12, 20, 2))
    sft <- pickSoftThreshold(dat, powerVector = powers, verbose = 0, networkType = "signed")
    fi <- sft$fitIndices
    d1 <- data.frame(power = fi$Power, fit = -sign(fi$slope) * fi$SFT.R.sq)
    d2 <- data.frame(power = fi$Power, conn = fi$mean.k.)
    p1 <- ggplot(d1, aes(power, fit)) + geom_hline(yintercept = 0.8, color = "#E64B35FF",
            linetype = "dashed", linewidth = 0.3) +
      geom_point(color = "#3C5488FF", size = 1.2) + geom_line(color = "#3C5488FF", linewidth = 0.4) +
      scale_x_continuous(breaks = powers) +
      labs(title = paste0(letter, "1  ", ds, " scale-free fit"), x = "Soft threshold (power)",
           y = expression("Scale-free topology fit " * R^2)) + theme_ajrccm()
    p2 <- ggplot(d2, aes(power, conn)) +
      geom_point(color = "#00A087FF", size = 1.2) + geom_line(color = "#00A087FF", linewidth = 0.4) +
      scale_x_continuous(breaks = powers) +
      labs(title = paste0(letter, "2  ", ds, " mean connectivity"), x = "Soft threshold (power)",
           y = "Mean connectivity") + theme_ajrccm()
    p1 + p2
  }
  save_pdf(sft_plot("GSE76293", "A"), file.path(OUT_SUPP, "FigS1A_soft_threshold_GSE76293.pdf"), 140, 62)
  save_pdf(sft_plot("GSE69528", "B"), file.path(OUT_SUPP, "FigS1B_soft_threshold_GSE69528.pdf"), 140, 62)
})

## ---- S2: ML out-of-fold ROC + confusion matrices ----
safe("S2", {
  oof_roc <- function(ds, letter) {
    o <- read.csv(file.path(B, "05_ML_SHAP", paste0(ds, "_OOF_predictions.csv")))
    y <- as.numeric(o$group != "Control")
    cls <- setdiff(colnames(o), c("sample", "group"))
    roc_list <- lapply(cls, function(m) roc(y, o[[m]], quiet = TRUE))
    names(roc_list) <- cls
    aucs <- sapply(roc_list, function(r) as.numeric(r$auc))
    ord <- names(sort(aucs, decreasing = TRUE))
    d <- do.call(rbind, lapply(ord, function(m) {
      r <- roc_list[[m]]
      data.frame(FPR = 1 - r$specificities, TPR = r$sensitivities,
                 model = paste0(m, " (", sprintf("%.2f", aucs[m]), ")"))
    }))
    d$model <- factor(d$model, levels = paste0(ord, " (", sprintf("%.2f", aucs[ord]), ")"))
    cols <- colorRampPalette(pal_npg(10))(length(ord))
    ggplot(d, aes(FPR, TPR, color = model)) + geom_line(linewidth = 0.5) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
      scale_color_manual(values = cols) + coord_equal() +
      labs(title = paste0(letter, "  ", ds, " out-of-fold ROC"), x = "1 - Specificity",
           y = "Sensitivity", color = NULL) +
      theme_ajrccm() + theme(legend.position = "right", legend.text = element_text(size = 5),
                             legend.key.size = unit(2.5, "mm"))
  }
  save_pdf(oof_roc("GSE76293", "A"), file.path(OUT_SUPP, "FigS2A_OOF_ROC_GSE76293.pdf"), 110, 80)
  save_pdf(oof_roc("GSE69528", "B"), file.path(OUT_SUPP, "FigS2B_OOF_ROC_GSE69528.pdf"), 110, 80)

  conf_heat <- function(ds, letter, model = "RF") {
    o <- read.csv(file.path(B, "05_ML_SHAP", paste0(ds, "_OOF_predictions.csv")))
    y <- factor(o$group); pred <- factor(ifelse(o[[model]] >= 0.5,
                                           levels(y)[2], levels(y)[1]), levels = levels(y))
    cm <- as.data.frame(table(Observed = y, Predicted = pred))
    ggplot(cm, aes(Predicted, Observed, fill = Freq)) +
      geom_tile(color = "white", linewidth = 0.4) +
      geom_text(aes(label = Freq), size = 2.8) +
      scale_fill_gradient(low = "white", high = "#3C5488FF") +
      labs(title = paste0(letter, "  ", ds, " ", model, " confusion matrix"), fill = "N") +
      theme_ajrccm() + theme(legend.position = "none")
  }
  save_pdf(conf_heat("GSE76293", "C"), file.path(OUT_SUPP, "FigS2C_confusion_GSE76293.pdf"), 62, 58)
  save_pdf(conf_heat("GSE69528", "D"), file.path(OUT_SUPP, "FigS2D_confusion_GSE69528.pdf"), 62, 58)
})

## ---- S3: MR method forest ----
safe("S3", {
  mr <- read.csv(file.path(B, "04_MR_SMR", "MR_sepsis4980_ARDS_results.csv"))
  mr$or <- exp(mr$b); mr$lo <- exp(mr$b - 1.96 * mr$se); mr$hi <- exp(mr$b + 1.96 * mr$se)
  mr <- mr[order(mr$b), ]
  mr$method <- factor(mr$method, levels = mr$method)
  mr$plab <- paste0("P = ", format.pval(mr$pval, digits = 2))
  p <- ggplot(mr, aes(or, method)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, linewidth = 0.4, color = "#3C5488FF") +
    geom_point(size = 1.8, color = "#E64B35FF") +
    geom_text(aes(x = hi, label = plab), hjust = -0.1, size = 2.2) +
    scale_x_log10(limits = c(min(mr$lo) * 0.8, max(mr$hi) * 2.2)) +
    labs(title = "MR: sepsis (ieu-b-4980) on ARDS, odds ratios by method",
         x = "OR (log scale, 95% CI)", y = NULL) +
    theme_ajrccm()
  save_pdf(p, file.path(OUT_SUPP, "FigS3_MR_method_forest.pdf"), 120, 68)
})

## ---- S4: ssGSEA score heatmap GSE69528 ----
safe("S4", {
  library(GSVA)
  immune_sets <- list(
    CD8_T = c("CD8A","CD8B","GZMA","GZMB"), CD4_T = c("CD4","IL7R","CD3D","CD3E"),
    Treg = c("FOXP3","IL2RA","CTLA4","TIGIT"), Th1 = c("TBX21","IFNG","STAT4","IL12RB2"),
    Th2 = c("GATA3","IL4","IL5","IL13"), Th17 = c("RORC","IL17A","IL22","IL23R"),
    NK = c("NKG7","KLRD1","NCAM1","GNLY"), B_cell = c("MS4A1","CD79A","CD19","CD79B"),
    Plasma = c("MZB1","XBP1","TNFRSF17","SDC1"), Monocyte = c("CD14","LYZ","FCN1","S100A12"),
    M1_Macro = c("IRF5","NOS2","TNF","IL1B"), M2_Macro = c("CD163","MRC1","ARG1","MSR1"),
    Neutrophil = c("FCGR3B","CSF3R","S100A8","CEACAM8"), DC = c("ITGAX","CD1C","FCER1A","CST3"),
    pDC = c("LILRA4","CLEC4C","IRF7","TCF4"), Mast = c("KIT","TPSAB1","CPA3","MS4A2"),
    Eosinophil = c("SIGLEC8","CCR3","RNASE2","IL5RA"), MDSC = c("CD33","FUT4","ITGAM","S100A9"))
  expr <- readRDS(file.path(B, "01_DEG", "GSE69528_gene_expr.rds"))
  grp <- readRDS(file.path(B, "01_DEG", "GSE69528_grp.rds"))
  sets <- lapply(immune_sets, function(g) intersect(g, rownames(expr)))
  sets <- sets[lengths(sets) >= 2]
  sc <- gsva(ssgseaParam(exprData = as.matrix(expr), geneSets = sets, minSize = 2))
  mz <- t(scale(t(as.matrix(sc))))
  grp <- grp[match(colnames(mz), grp$sample), ]
  ord <- order(grp$group)
  long <- expand.grid(cell = rownames(mz), sample = colnames(mz)[ord], stringsAsFactors = FALSE)
  long$z <- as.vector(mz[, ord])
  long$cell <- factor(long$cell, levels = rev(rownames(mz)))
  long$sample <- factor(long$sample, levels = colnames(mz)[ord])
  long$group <- grp$group[ord][match(long$sample, colnames(mz)[ord])]
  br <- cumsum(table(grp$group[ord]))
  p <- ggplot(long, aes(sample, cell, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = "#3C5488FF", mid = "white", high = "#E64B35FF",
                         midpoint = 0, name = "z-score") +
    labs(title = "GSE69528 ssGSEA immune-cell scores (row z-score)", x = NULL, y = NULL) +
    theme_ajrccm() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.key.height = unit(4, "mm"), legend.title = element_text(size = 6))
  save_pdf(p, file.path(OUT_SUPP, "FigS4_ssGSEA_heatmap_GSE69528.pdf"), 150, 70)
})

## ---- S6: scTenifoldKnk KO top-10 differentially regulated genes ----
safe("S6", {
  ko_bar <- function(gene, letter) {
    x <- read.csv(file.path(SC, "knockout", paste0("KO_", gene, "_diffRegulation.csv")))
    cn <- colnames(x); gcol <- cn[1]
    zcol <- cn[grepl("Z|diff", cn, ignore.case = TRUE)][1]
    x <- x[order(-abs(x[[zcol]])), ][1:10, ]
    x <- x[order(x[[zcol]]), ]
    x[[gcol]] <- factor(x[[gcol]], levels = x[[gcol]])
    ggplot(x, aes_string(gcol, zcol)) +
      geom_col(fill = "#3C5488FF", width = 0.65) + coord_flip() +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      labs(title = substitute(paste(l, "  ", italic(g), " KO top 10 (scTenifoldKnk)",
                                    list(l = letter, g = gene))),
           x = NULL, y = "Differential regulation (Z)") +
      theme_ajrccm()
  }
  save_pdf(ko_bar("UPP1", "A"), file.path(OUT_SUPP, "FigS6A_KO_UPP1_top10.pdf"), 80, 62)
  save_pdf(ko_bar("EXOSC4", "B"), file.path(OUT_SUPP, "FigS6B_KO_EXOSC4_top10.pdf"), 80, 62)
})

## ---- S7: HPA infection logFC ----
safe("S7", {
  hpa <- read.csv(file.path(B, "10_HPA_validation", "HPA_CD177_VNN2_infection_vs_healthy.csv"))
  hpa$Disease <- factor(hpa$Disease, levels = hpa$Disease[order(hpa$logFC[hpa$Gene == "CD177"])])
  p <- ggplot(hpa, aes(Disease, logFC, fill = Gene)) +
    geom_col(position = position_dodge(0.75), width = 0.7) +
    geom_text(aes(label = sig, group = Gene), vjust = -0.4, size = 2.2,
              position = position_dodge(0.75)) +
    scale_fill_manual(values = pal_npg(2)) +
    labs(title = "HPA Blood Atlas: infection vs healthy (log2 fold change)",
         x = NULL, y = "log2FC", fill = NULL) +
    theme_ajrccm() + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 5.5),
                           legend.position = "top")
  save_pdf(p, file.path(OUT_SUPP, "FigS7_HPA_infection_logFC.pdf"), 150, 68)
})

## ---- S9: GSE167363 atlas - hub mean-expression heatmap by cell type ----
safe("S9", {
  m <- read.csv(file.path(SC, "GSE167363", "GSE167363_hub_celltype_mean.csv"), check.names = FALSE)
  long <- pivot_longer(m, cols = all_of(HUB6), names_to = "gene", values_to = "mean_expr")
  long$gene <- factor(long$gene, levels = rev(HUB6))
  long$celltype <- factor(long$celltype, levels = sort(unique(long$celltype)))
  p <- ggplot(long, aes(celltype, gene, fill = log1p(mean_expr))) +
    geom_tile(color = "white", linewidth = 0.2) +
    facet_wrap(~ group, nrow = 1) +
    scale_fill_gradient(low = "white", high = "#E64B35FF", name = "log(1+expr)") +
    scale_y_discrete(labels = function(x) parse(text = gene_italic(x))) +
    labs(title = "GSE167363 hub-gene mean expression by cell type", x = NULL, y = NULL) +
    theme_ajrccm() + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 5.5),
                           legend.key.height = unit(4, "mm"), legend.title = element_text(size = 6))
  save_pdf(p, file.path(OUT_SUPP, "FigS9_hub_expr_heatmap_GSE167363.pdf"), 130, 62)
})

## ---- S10: AUCell cell-type heatmap + GSE167363 CellOracle TF bars ----
safe("S10", {
  ac <- read.csv(file.path(SC, "GSE167363", "aucell", "GSE167363_AUCell_scores.csv"),
                 check.names = FALSE)
  progs <- c(paste0(HUB6, "_program"), "HubSignature")
  agg <- ac %>% group_by(celltype, group) %>%
    summarise(across(all_of(progs), mean), .groups = "drop")
  long <- pivot_longer(agg, cols = all_of(progs), names_to = "program", values_to = "score")
  long$program <- gsub("_program", "", long$program)
  long$program <- factor(long$program, levels = rev(c(HUB6, "HubSignature")))
  p1 <- ggplot(long, aes(celltype, program, fill = score)) +
    geom_tile(color = "white", linewidth = 0.2) + facet_wrap(~ group, nrow = 1) +
    scale_fill_gradient(low = "white", high = "#3C5488FF", name = "AUCell") +
    labs(title = "A  GSE167363 hub-program activity by cell type", x = NULL, y = NULL) +
    theme_ajrccm() + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 5.5),
                           legend.key.height = unit(4, "mm"), legend.title = element_text(size = 6))
  save_pdf(p1, file.path(OUT_SUPP, "FigS10A_AUCell_heatmap_GSE167363.pdf"), 140, 66)

  ko <- read.csv(file.path(SC, "GSE167363", "celloracle", "CellOracle_TF_KO_hub_effect.csv"))
  tgs <- c("UPP1", "CDK5RAP2", "EXOSC4", "VNN1")
  ko <- ko[ko$target %in% tgs, ]
  top_tf <- ko %>% group_by(target) %>% slice_min(delta, n = 8, with_ties = FALSE) %>%
    ungroup()
  top_tf$KO_TF <- factor(top_tf$KO_TF, levels = rev(unique(top_tf$KO_TF[order(top_tf$delta)])))
  p2 <- ggplot(top_tf, aes(KO_TF, delta)) +
    geom_col(fill = "#3C5488FF", width = 0.65) + coord_flip() +
    facet_wrap(~ target, scales = "free_y", nrow = 1,
               labeller = labeller(target = function(x) paste0("italic('", x, "')"))) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    labs(title = "B  GSE167363 CellOracle upstream TF knockout", x = NULL,
         y = "Simulated expression change (delta)") +
    theme_ajrccm() + theme(strip.text = element_text(size = 6),
                           axis.text.y = element_text(size = 5))
  save_pdf(p2, file.path(OUT_SUPP, "FigS10B_CellOracle_TF_GSE167363.pdf"), 178, 62)
})

## ---- S12: metabolic pathway group differences ----
safe("S12", {
  m1 <- read.csv(file.path(SC, "metabolism", "metabolism_group_stats_GSE151263.csv"))
  m1$dataset <- "GSE151263 (PBMC)"
  m2 <- read.csv(file.path(SC, "metabolism", "metabolism_group_stats_GSE167363.csv"))
  m2$dataset <- "GSE167363 (CD14+ mono)"
  d <- rbind(m1, m2)
  d$pathway <- factor(d$pathway, levels = rev(sort(unique(d$pathway))))
  p <- ggplot(d, aes(pathway, delta, fill = dataset)) +
    geom_col(position = position_dodge(0.75), width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    coord_flip() + scale_fill_manual(values = pal_npg(2)) +
    labs(title = "Metabolic pathway scores, disease vs control (delta)", x = NULL,
         y = "Mean score difference", fill = NULL) +
    theme_ajrccm() + theme(legend.position = "top", legend.text = element_text(size = 6))
  save_pdf(p, file.path(OUT_SUPP, "FigS12_metabolic_delta.pdf"), 120, 72)
})

## ---- S13: stratified knockout comparison ----
safe("S13", {
  cmp <- read.csv(file.path(SC, "knockout_stratified", "KO_stratified_group_comparison.csv"))
  cmp$gene <- factor(cmp$gene, levels = HUB6)
  p1 <- ggplot(cmp, aes(gene, spearman_rho, fill = dataset)) +
    geom_col(position = position_dodge(0.75), width = 0.7) +
    scale_fill_manual(values = pal_npg(2)) +
    scale_x_discrete(labels = function(x) parse(text = gene_italic(x))) +
    labs(title = "A  KO effect concordance between groups (Spearman rho)",
         x = NULL, y = "Spearman rho", fill = NULL) +
    theme_ajrccm() + theme(legend.position = "top", legend.text = element_text(size = 6))
  save_pdf(p1, file.path(OUT_SUPP, "FigS13A_KO_concordance.pdf"), 110, 62)

  co <- bind_rows(read.csv(file.path(SC, "celloracle_stratified",
                                 "CellOracle_TF_KO_hub_effect_GSE167363_HC.csv")),
              read.csv(file.path(SC, "celloracle_stratified",
                                 "CellOracle_TF_KO_hub_effect_GSE167363_Sepsis.csv")),
              read.csv(file.path(SC, "celloracle_stratified",
                                 "CellOracle_TF_KO_hub_effect_GSE151263_Sepsis.csv")),
              read.csv(file.path(SC, "celloracle_stratified",
                                 "CellOracle_TF_KO_hub_effect_GSE151263_ARDS.csv")))
  cebp <- co[co$KO_TF == "CEBPB" & co$target == "UPP1", ]
  cebp$grp <- factor(paste0(cebp$dataset, "_", cebp$group),
                     levels = c("GSE167363_HC", "GSE167363_Sepsis",
                                "GSE151263_Sepsis", "GSE151263_ARDS"))
  p2 <- ggplot(cebp, aes(grp, delta)) +
    geom_col(fill = "#3C5488FF", width = 0.55) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", delta)), vjust = -0.4, size = 2.2) +
    labs(title = "B  CEBPB knockout effect on UPP1 by disease group", x = NULL,
         y = "Simulated delta") +
    theme_ajrccm() + theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 5.5))
  save_pdf(p2, file.path(OUT_SUPP, "FigS13B_CEBPB_UPP1_gradient.pdf"), 90, 62)
})

## ---- S14: donor-level sensitivity ----
safe("S14", {
  d <- read.csv(file.path(SC, "donor_sensitivity", "donor_level_summary.csv"))
  pt <- d[d$metric == "pseudotime", ]
  pt$group <- factor(pt$group, levels = c("HC", "Sepsis", "ARDS"))
  p1 <- ggplot(pt, aes(group, value, color = group)) +
    geom_point(size = 2, alpha = 0.85) +
    facet_wrap(~ dataset, scales = "free", nrow = 1) +
    scale_color_manual(values = pal_npg(3)) +
    labs(title = "A  Donor-level pseudotime", x = NULL, y = "Median pseudotime") +
    theme_ajrccm() + theme(legend.position = "none",
                           axis.text.x = element_text(angle = 25, hjust = 1))
  up <- d[d$metric == "UPP1_program", ]
  up$group <- factor(up$group, levels = c("HC", "Sepsis", "ARDS"))
  p2 <- ggplot(up, aes(group, value, color = group)) +
    geom_point(size = 2, alpha = 0.85) +
    facet_wrap(~ dataset, scales = "free", nrow = 1) +
    scale_color_manual(values = pal_npg(3)) +
    labs(title = substitute(paste("B  Donor-level ", italic("UPP1"), " program")),
         x = NULL, y = "Mean AUCell score") +
    theme_ajrccm() + theme(legend.position = "none",
                           axis.text.x = element_text(angle = 25, hjust = 1))
  gly <- d[d$metric == "metabolism" & d$pathway == "Glycolysis", ]
  gly$group <- factor(gly$group, levels = c("HC", "Sepsis", "ARDS"))
  p3 <- ggplot(gly, aes(group, value, color = group)) +
    geom_point(size = 2, alpha = 0.85) +
    facet_wrap(~ dataset, scales = "free", nrow = 1) +
    scale_color_manual(values = pal_npg(3)) +
    labs(title = "C  Donor-level glycolysis score", x = NULL, y = "Mean score") +
    theme_ajrccm() + theme(legend.position = "none",
                           axis.text.x = element_text(angle = 25, hjust = 1))
  ox <- d[d$metric == "metabolism" & d$pathway == "Oxidative phosphorylation", ]
  ox$group <- factor(ox$group, levels = c("HC", "Sepsis", "ARDS"))
  p4 <- ggplot(ox, aes(group, value, color = group)) +
    geom_point(size = 2, alpha = 0.85) +
    facet_wrap(~ dataset, scales = "free", nrow = 1) +
    scale_color_manual(values = pal_npg(3)) +
    labs(title = "D  Donor-level oxidative phosphorylation", x = NULL, y = "Mean score") +
    theme_ajrccm() + theme(legend.position = "none",
                           axis.text.x = element_text(angle = 25, hjust = 1))
  save_pdf((p1 + p2) / (p3 + p4), file.path(OUT_SUPP, "FigS14_donor_sensitivity.pdf"), 178, 120)
})

cat("Supplementary panels done.\n")
