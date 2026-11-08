# ---------------------------------------------------------------------------
# 06_fig6.R -- Figure 6, external confirmation of the six-gene signature
#   A  single-gene AUC with DeLong confidence intervals in the two validation
#      cohorts
#   B  forest plot of the standardised mean difference (Hedges g) per gene,
#      with the fixed-effect pooled estimate
#   C  ROC curves in the ARDS validation cohort (GSE243066)
#   D  ROC curves in the sepsis validation cohort (GSE100159)
# Input  : data/derived/hub6_expression_long.csv,
#          data/frozen/06_external_validation/AUC_CI.csv
# Output : Figures/main/Figure6.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")
suppressPackageStartupMessages(library(pROC))

long <- read.csv(file.path("data", "derived", "hub6_expression_long.csv"),
                 stringsAsFactors = FALSE)
ci <- read.csv(file.path("data", "frozen", "06_external_validation", "AUC_CI.csv"),
               stringsAsFactors = FALSE)
ci <- ci[ci$gene %in% HUB6, ]
ci$cohort <- factor(ci$cohort, levels = c("GSE100159_Sepsis", "GSE243066_ARDS"),
                    labels = c("Sepsis (GSE100159)", "ARDS (GSE243066)"))
ci$gene <- factor(ci$gene, levels = rev(HUB6))

# --- panel A ---------------------------------------------------------------

panel_a <- function() {
  ggplot(ci, aes(AUC, gene, colour = cohort)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.18,
                   linewidth = 0.4, position = position_dodge(width = 0.45)) +
    geom_point(size = 1.1, position = position_dodge(width = 0.45)) +
    scale_colour_manual(values = c(pal[["steel"]], pal[["brick"]]), name = NULL) +
    scale_x_continuous(limits = c(0.4, 1.02), breaks = seq(0.5, 1, 0.1)) +
    labs(x = "AUC (95% CI, DeLong)", y = NULL) +
    theme_ajr() +
    theme(axis.text.y = element_text(face = "italic"),
          legend.position = "top", legend.direction = "horizontal",
          legend.margin = margin(0, 0, -3, 0))
}

# --- panel B ---------------------------------------------------------------
# Pooled values are taken from the frozen meta-analysis table
# (06_external_validation/reviewer_w5/metaDE_single_gene_validation.csv), so the
# figure reports the same fixed-effect estimates and I2 values as the text
# (g 1.14 to 2.32; I2 = 0.60 to 0.91 for VNN1, UPP1 and CD177). Per-cohort
# standard errors use the Hedges variance approximation that reproduces the
# frozen pooled confidence intervals to two decimals.

panel_b <- function() {
  m <- read.csv(file.path("data", "frozen", "06_external_validation", "reviewer_w5",
                          "metaDE_single_gene_validation.csv"), stringsAsFactors = FALSE)
  m <- m[m$gene %in% HUB6, ]

  smp <- unique(long[long$cohort %in% c("GSE243066", "GSE100159"),
                     c("sample", "cohort", "group")])
  nn <- table(smp$cohort, smp$group)
  se_g <- function(g, n_dis, n_ctl) {
    sqrt((n_dis + n_ctl) / (n_dis * n_ctl) + g^2 / (2 * (n_dis + n_ctl - 4)))
  }

  per <- rbind(
    data.frame(gene = m$gene, cohort = "ARDS (GSE243066)", g = m$g_GSE243066,
               se = se_g(m$g_GSE243066, nn["GSE243066", "Disease"],
                         nn["GSE243066", "Control"])),
    data.frame(gene = m$gene, cohort = "Sepsis (GSE100159)", g = m$g_GSE100159,
               se = se_g(m$g_GSE100159, nn["GSE100159", "Disease"],
                         nn["GSE100159", "Control"]))
  )
  split_ci <- strsplit(as.character(m$pooled_g_FE_ci), "-")
  pooled <- data.frame(gene = m$gene, g = m$pooled_g_FE,
                       lo = as.numeric(vapply(split_ci, "[", "", 1)),
                       hi = as.numeric(vapply(split_ci, "[", "", 2)),
                       i2 = m$I2)

  per$gene <- factor(per$gene, levels = rev(HUB6))
  per$cohort <- factor(per$cohort, levels = c("Sepsis (GSE100159)", "ARDS (GSE243066)"))
  pooled$gene <- factor(pooled$gene, levels = rev(HUB6))
  pooled$y <- as.numeric(pooled$gene) + 0.26
  pooled$lab <- sprintf("I2 = %.2f", pooled$i2)
  xr <- max(pooled$hi)

  ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    geom_errorbarh(data = per, aes(y = as.numeric(gene) - 0.14, xmin = g - 1.96 * se,
                                   xmax = g + 1.96 * se, colour = cohort),
                   height = 0, linewidth = 0.35) +
    geom_point(data = per, aes(g, as.numeric(gene) - 0.14, colour = cohort),
               size = 0.9) +
    geom_errorbarh(data = pooled, aes(y = y, xmin = lo, xmax = hi),
                   height = 0, colour = pal[["ink"]], linewidth = 0.35) +
    geom_point(data = pooled, aes(g, y), shape = 23, size = 1.5,
               fill = pal[["ink"]], colour = pal[["ink"]]) +
    geom_text(data = pooled, aes(x = xr + 0.55, y = y, label = lab), hjust = 0,
              family = FONT_BODY, size = 7 / .pt, colour = pal[["ink"]]) +
    annotate("text", x = xr + 0.55, y = length(HUB6) + 0.62, label = "Heterogeneity",
             hjust = 0, family = FONT_BODY, size = 7 / .pt, colour = pal[["ink"]]) +
    scale_colour_manual(values = c(pal[["steel"]], pal[["brick"]]), name = NULL) +
    scale_y_continuous(breaks = seq_along(HUB6), labels = rev(HUB6),
                       limits = c(0.45, length(HUB6) + 0.75)) +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.26))) +
    labs(x = "Standardised mean difference (Hedges g)", y = NULL) +
    theme_ajr() +
    theme(axis.text.y = element_text(face = "italic"),
          legend.position = "top", legend.direction = "horizontal",
          legend.margin = margin(0, 0, -3, 0))
}

# --- panels C, D -----------------------------------------------------------

panel_roc <- function(cd) {
  d <- long[long$cohort == cd, ]
  smp <- unique(d[, c("sample", "group")])
  wide <- tapply(d$expr, list(d$sample, d$gene), mean)
  rownames(wide) <- rownames(wide)
  y <- as.integer(smp$group[match(rownames(wide), smp$sample)] == "Disease")

  curves <- lapply(HUB6, function(g) {
    r <- roc(y, wide[, g], quiet = TRUE, direction = "<")
    data.frame(fpr = 1 - r$specificities, tpr = r$sensitivities,
               grp = sprintf("%s (%.3f)", g, as.numeric(auc(r))))
  })
  curves <- do.call(rbind, curves)
  curves$grp <- factor(curves$grp, levels = unique(curves$grp))
  cols <- setNames(hub_col, paste0(HUB6, " (", c("1", "1", "1", "1", "1", "1"), ")"))
  cols <- setNames(hub_col[seq_along(unique(curves$grp))], unique(curves$grp))

  n_dis <- sum(y == 1); n_ctl <- sum(y == 0)
  ggplot(curves, aes(fpr, tpr, colour = grp)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3,
                colour = pal[["slate"]]) +
    geom_step(linewidth = 0.4) +
    scale_colour_manual(values = unname(hub_col[seq_along(levels(curves$grp))]),
                        name = NULL) +
    scale_x_continuous(expand = expansion(mult = 0.01)) +
    scale_y_continuous(expand = expansion(mult = 0.01)) +
    labs(x = "1 - specificity", y = "Sensitivity",
         subtitle = sprintf("%s validation (n = %d vs %d)", cd, n_dis, n_ctl)) +
    theme_ajr() +
    theme(legend.position = "right",
          legend.key.height = unit(6, "pt"),
          legend.text = element_text(size = PT_BODY),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

# --- assemble --------------------------------------------------------------

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid,
                c(list(panel_a(), panel_b(), ncol = 2, labels = c("A", "B")), lab_args))
row2 <- do.call(cowplot::plot_grid,
                c(list(panel_roc("GSE243066"), panel_roc("GSE100159"),
                       ncol = 2, labels = c("C", "D")), lab_args))

fig6 <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1.05))

save_fig(fig6, "Figure6", width_mm = FIG_W_MM, height_mm = 200)
msg("Figure6 written")
