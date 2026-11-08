# ---------------------------------------------------------------------------
# 42_supp_S22.R -- Figure S22. Classifier performance and feature attribution in
# both discovery cohorts
#   A  cross-validated AUC of the eleven estimable classifiers
#   B  SHAP importance, ARDS discovery cohort
#   C  SHAP importance, sepsis discovery cohort
# Input  : data/frozen/05_ML_SHAP
# Output : Figures/supplementary/FigureS22.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

ML <- function(...) file.path("data", "frozen", "05_ML_SHAP", ...)

# --- panel A ---------------------------------------------------------------

panel_a <- function() {
  a <- read.csv(ML("ML_AUC_results.csv"), stringsAsFactors = FALSE)
  b <- read.csv(ML("sepsis", "ML_AUC_results_sepsis.csv"), stringsAsFactors = FALSE)
  d <- rbind(data.frame(model = a$model, auc = a$mean_AUC, sd = a$sd_AUC,
                        cohort = "ARDS discovery"),
             data.frame(model = b$model, auc = b$mean_AUC, sd = b$sd_AUC,
                        cohort = "Sepsis discovery"))
  d$cohort <- factor(d$cohort, levels = c("ARDS discovery", "Sepsis discovery"))
  ord <- names(sort(tapply(d$auc, d$model, mean)))
  d$model <- factor(d$model, levels = ord)

  ggplot(d, aes(auc, model, fill = cohort)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.62, colour = NA) +
    scale_fill_manual(values = c(`ARDS discovery` = pal[["brick"]],
                                 `Sepsis discovery` = pal[["steel"]]), name = NULL) +
    scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25),
                       expand = expansion(mult = c(0, 0.03))) +
    labs(x = "Cross-validated AUC", y = NULL) +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"))
}

# --- panel B, C ------------------------------------------------------------

shap_panel <- function(file, sub) {
  d <- read.csv(ML(file), stringsAsFactors = FALSE)
  d <- d[order(-d$mean_abs_shap), ][1:12, ]
  d$hub <- d$gene %in% HUB6
  d$gene <- factor(d$gene, levels = rev(d$gene))

  ggplot(d, aes(mean_abs_shap, gene, fill = hub)) +
    geom_col(width = 0.64, colour = NA) +
    scale_fill_manual(values = c(`TRUE` = pal[["brick"]], `FALSE` = pal[["slate"]]),
                      labels = c(`TRUE` = "six-gene hub", `FALSE` = "other feature"),
                      name = NULL) +
    scale_y_discrete(labels = function(v) v) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Mean absolute SHAP value", y = NULL, subtitle = sub) +
    theme_ajr() +
    theme(axis.text.y = element_text(face = "italic"),
          legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid, c(list(panel_a(), labels = "A"), lab_args))
row2 <- do.call(cowplot::plot_grid,
                c(list(shap_panel("SHAP_importance.csv", "ARDS discovery"),
                       shap_panel(file.path("sepsis", "SHAP_importance_sepsis.csv"),
                                  "Sepsis discovery"),
                       ncol = 2, labels = c("B", "C")), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1.1))
save_fig(fig, "FigureS22", dir = out_dir("supplementary"), height_mm = 206)
msg("FigureS22 written")
