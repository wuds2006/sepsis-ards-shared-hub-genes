# ---------------------------------------------------------------------------
# 02_fig2.R -- Figure 2, core results
#   A  cross-validated performance of the machine-learning classifiers
#   B  SHAP importance, top 12 of the 179 shared candidates
#   C  external validation of the six hub genes in the two validation cohorts
#   D  M2 macrophage infiltration score in the discovery cohorts
#   E  shift in immune-cell scores, both discovery cohorts
# Inputs : data/frozen/05_ML_SHAP, 06_external_validation, 07_hub_immune
# Output : Figures/main/Figure2.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

FZ <- function(...) file.path("data", "frozen", ...)

# --- panel A ---------------------------------------------------------------

ml_a <- function() {
  auc <- read.csv(FZ("05_ML_SHAP", "ML_AUC_results.csv"), stringsAsFactors = FALSE)
  est <- auc[!is.na(auc$mean_AUC), ]
  est <- est[order(est$mean_AUC), ]
  gap <- auc[is.na(auc$mean_AUC), ]

  lv <- c(paste0(gap$model, "\u2020"), est$model)
  est$model <- factor(est$model, levels = lv)
  gap$model <- factor(paste0(gap$model, "\u2020"), levels = lv)

  p <- ggplot() +
    geom_col(data = est, aes(mean_AUC, model), width = 0.66,
             position = "identity", fill = pal[["steel"]], colour = NA) +
    geom_errorbarh(data = est, aes(y = model, xmin = mean_AUC - sd_AUC,
                                   xmax = mean_AUC + sd_AUC),
                   height = 0.22, linewidth = 0.35, colour = pal[["ink"]]) +
    geom_text(data = est, aes(mean_AUC + sd_AUC + 0.035, model,
                              label = sprintf("%.2f", mean_AUC)),
              hjust = 0, family = FONT_BODY, size = 7 / .pt,
              colour = pal[["ink"]]) +
    geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    scale_x_continuous(breaks = seq(0, 1, 0.25),
                       expand = expansion(mult = c(0, 0.14))) +
    scale_y_discrete(expand = expansion(add = 0.42)) +
    labs(x = "Mean cross-validated AUC", y = NULL,
         subtitle = sprintf("%d estimable classifiers, 179 shared candidates",
                            nrow(est)),
         caption = "\u2020 GBM not estimable at 19 training samples per fold") +
    theme_ajr() +
    theme(axis.text.y = element_text(size = PT_BODY),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["slate"]]),
          plot.caption = element_text(size = PT_BODY, hjust = 0,
                                      colour = pal[["slate"]]),
          plot.caption.position = "plot")

  if (nrow(gap)) {
    p <- p + annotate("rect", xmin = 0, xmax = 0.045, ymin = 0.62,
                      ymax = 1.38, fill = NA, colour = pal[["slate"]],
                      linewidth = 0.3)
  }
  p
}

# --- panel B ---------------------------------------------------------------

shap_b <- function() {
  s <- read.csv(FZ("05_ML_SHAP", "SHAP_importance.csv"), stringsAsFactors = FALSE)
  s <- head(s[order(-s$mean_abs_shap), ], 12)
  s$hub <- s$gene %in% HUB6
  s$gene <- factor(s$gene, levels = rev(s$gene[order(s$mean_abs_shap)]))

  ggplot(s, aes(mean_abs_shap, gene)) +
    geom_col(aes(fill = hub), width = 0.68, position = "identity", colour = NA) +
    scale_fill_manual(values = c(`TRUE` = pal[["brick"]], `FALSE` = "#CFCBC4"),
                      guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Mean absolute SHAP value", y = NULL,
         subtitle = "Top 12 genes; hub genes in red") +
    theme_ajr() +
    theme(axis.text.y = element_text(size = PT_BODY, face = "italic"),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["slate"]]))
}

# --- panel C ---------------------------------------------------------------

ext_c <- function() {
  ci <- read.csv(FZ("06_external_validation", "AUC_CI.csv"), stringsAsFactors = FALSE)
  ord <- c(HUB6, "Combined_signature")
  ci <- ci[ci$gene %in% ord, ]
  ci$gene <- factor(ci$gene, levels = ord,
                    labels = c(HUB6, "6-gene\nsignature"))
  ci$cohort <- factor(ci$cohort,
                      levels = c("GSE100159_Sepsis", "GSE243066_ARDS"),
                      labels = c("Sepsis (GSE100159)", "ARDS (GSE243066)"))

  ggplot(ci, aes(gene, AUC, fill = cohort)) +
    geom_col(position = position_dodge(width = 0.74), width = 0.7,
             colour = NA) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                  position = position_dodge(width = 0.74), width = 0.16,
                  linewidth = 0.35, colour = pal[["ink"]]) +
    geom_text(aes(label = sprintf("%.2f", AUC), y = pmin(CI_upper + 0.05, 1.17)),
              position = position_dodge(width = 0.74), vjust = 0,
              family = FONT_BODY, size = 7 / .pt, colour = pal[["ink"]]) +
    geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    scale_fill_manual(values = c(pal[["steel"]], pal[["brick"]]), name = NULL) +
    scale_y_continuous(limits = c(0, 1.24), breaks = seq(0, 1, 0.25),
                       expand = expansion(mult = c(0, 0.01))) +
    labs(x = NULL, y = "AUC (95% CI)",
         subtitle = "Single-gene discrimination, no model fitting") +
    theme_ajr() +
    theme(legend.position = "top", legend.direction = "horizontal",
          legend.margin = margin(0, 0, -2, 0),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["slate"]]))
}

# --- panel D ---------------------------------------------------------------

m2_d <- function() {
  rd <- function(f) {
    d <- read.csv(FZ("07_hub_immune", "ssgsea_immune", f), stringsAsFactors = FALSE)
    d[, c("group", "M2_Macro")]
  }
  a <- rd("GSE76293_ssgsea_scores.csv"); a$cohort <- "ARDS (GSE76293)"
  s <- rd("GSE69528_ssgsea_scores.csv"); s$cohort <- "Sepsis (GSE69528)"
  d <- rbind(a, s)
  d$grp <- factor(paste(d$cohort, d$group, sep = "|"),
                  levels = c("ARDS (GSE76293)|Control", "ARDS (GSE76293)|ARDS",
                             "Sepsis (GSE69528)|Control", "Sepsis (GSE69528)|Sepsis"))
  d$fill <- c("ARDS (GSE76293)|Control" = pal[["sage"]],
              "ARDS (GSE76293)|ARDS" = pal[["brick"]],
              "Sepsis (GSE69528)|Control" = pal[["sage"]],
              "Sepsis (GSE69528)|Sepsis" = pal[["steel"]])[as.character(d$grp)]

  lab <- data.frame(grp = levels(d$grp),
                    x = seq_len(4),
                    txt = c("Control", "ARDS", "Control", "Sepsis"))
  n <- as.data.frame(table(d$grp))
  lab$txt <- sprintf("%s\nn = %d", lab$txt, n$Freq[match(lab$grp, n$Var1)])

  rng <- range(d$M2_Macro); span <- diff(rng)
  br <- data.frame(x1 = c(1, 3), x2 = c(2, 4),
                   y = c(max(d$M2_Macro[d$grp == levels(d$grp)[2]]) + span * 0.10,
                         max(d$M2_Macro[d$grp == levels(d$grp)[4]]) + span * 0.10),
                   p = c("P = 7.4e-7", "P = 5.8e-14"))

  ggplot(d, aes(as.integer(grp), M2_Macro)) +
    geom_boxplot(aes(fill = fill, group = grp), width = 0.58, outlier.shape = NA,
                 linewidth = 0.35, colour = "#4A4A4A") +
    geom_point(position = position_jitter(width = 0.11, height = 0, seed = 7),
               size = 0.55, colour = "#4A4A4A", stroke = 0) +
    scale_fill_identity() +
    geom_segment(data = br, aes(x = x1, xend = x1, y = y - span * 0.01,
                                yend = y), inherit.aes = FALSE,
                 linewidth = 0.35, colour = pal[["ink"]]) +
    geom_segment(data = br, aes(x = x2, xend = x2, y = y - span * 0.01,
                                yend = y), inherit.aes = FALSE,
                 linewidth = 0.35, colour = pal[["ink"]]) +
    geom_segment(data = br, aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.35, colour = pal[["ink"]]) +
    geom_text(data = br, aes(x = (x1 + x2) / 2, y = y + span * 0.02, label = p),
              inherit.aes = FALSE, family = FONT_BODY, size = 7 / .pt,
              vjust = 0, colour = pal[["ink"]]) +
    scale_x_continuous(breaks = 1:4, labels = lab$txt, expand = expansion(add = 0.55)) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.10))) +
    labs(x = NULL, y = "M2 macrophage ssGSEA score",
         subtitle = "Two-sided Wilcoxon rank-sum, discovery cohorts") +
    theme_ajr() +
    theme(plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["slate"]]))
}

# --- panel E ---------------------------------------------------------------

immune_e <- function() {
  frames <- lapply(c("GSE76293", "GSE69528"), function(tag) {
    d <- read.csv(FZ("07_hub_immune", paste0(tag, "_immune_diff.csv")),
                  stringsAsFactors = FALSE)
    d$padj <- p.adjust(d$p, "BH")
    d$cohort <- ifelse(tag == "GSE76293", "ARDS (GSE76293)", "Sepsis (GSE69528)")
    d
  })
  imm <- do.call(rbind, frames)
  ord <- names(sort(tapply(abs(imm$log2FC), imm$cell, max)))
  imm$cell <- factor(imm$cell, levels = ord)
  imm$cohort <- factor(imm$cohort, levels = c("ARDS (GSE76293)", "Sepsis (GSE69528)"))
  imm$sig <- imm$padj < 0.05
  imm$nlp <- pmax(-log10(imm$padj), 0)

  ggplot(imm, aes(log2FC, cell)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    geom_point(aes(colour = cohort, size = nlp), shape = 16) +
    scale_colour_manual(values = c(pal[["brick"]], pal[["steel"]]), name = NULL) +
    scale_size_area(max_size = 2.2, breaks = c(2, 6, 10, 14),
                    name = "-log10 adjusted P") +
    labs(x = "log2 fold change (disease vs control)", y = NULL) +
    theme_ajr() +
    theme(legend.position = "top", legend.box = "vertical",
          legend.margin = margin(0, 0, -4, 0),
          legend.spacing.x = unit(3, "mm"),
          legend.key.size = unit(6, "pt"),
          legend.box.just = "left",
          legend.justification = "left")
}

# --- assemble --------------------------------------------------------------

# Panel letters go through plot_annotation() rather than tag_plot(): tag_plot()
# pins the tag at a fixed npc position above the panel, and for the top row
# (panels with a subtitle and no headroom) that position falls outside the
# canvas, so A, B and C were clipped in the exported files. plot_annotation()
# reserves tag space inside each panel margin instead.
fig2 <- (ml_a() | shap_b()) / ext_c() / (m2_d() | immune_e()) +
  plot_layout(heights = c(1, 0.95, 1.12)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.tag = element_text(family = FONT_TAG, size = PT_TAG, face = "bold",
                                colour = pal[["ink"]]))

save_fig(fig2, "Figure2", width_mm = FIG_W_MM, height_mm = 214)
msg("Figure2 written")
