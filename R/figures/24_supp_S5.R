# ---------------------------------------------------------------------------
# 24_supp_S5.R -- Figure S5. Individual machine-learning classifier performance
#   A  out-of-fold ROC curves, ARDS discovery cohort
#   B  out-of-fold ROC curves, sepsis discovery cohort
#   C  out-of-fold predicted probability by group, best-performing classifier
# Input  : data/frozen/05_ML_SHAP/GSE76293_OOF_predictions.csv, GSE69528_...
# Output : Figures/supplementary/FigureS5.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")
suppressPackageStartupMessages(library(pROC))

ML <- function(...) file.path("data", "frozen", "05_ML_SHAP", ...)
MODELS <- c("LASSO", "Ridge", "ENet", "RF", "SVM", "KNN", "NB", "NN", "LR", "XGB", "DT")

curve_panel <- function(file, sub) {
  d <- read.csv(ML(file), stringsAsFactors = FALSE)
  y <- as.integer(d$group != "Control")
  fits <- lapply(MODELS, function(m) roc(y, d[[m]], quiet = TRUE, direction = "<"))
  aucs <- vapply(fits, function(r) as.numeric(auc(r)), 1)
  ord <- order(-aucs)
  labs <- sprintf("%s (%.3f)", MODELS, aucs)

  out <- do.call(rbind, lapply(seq_along(MODELS), function(i) {
    r <- fits[[i]]
    data.frame(fpr = 1 - r$specificities, tpr = r$sensitivities, lab = labs[i])
  }))
  out$lab <- factor(out$lab, levels = labs[ord])
  cols <- setNames(c(pal[["brick"]], pal[["steel"]], pal[["sage"]], pal[["ochre"]],
                     pal[["plum"]], pal[["olive"]], pal[["slate"]], "#8C6D5B",
                     "#5F7A6B", "#A97C83", "#7A7F8C")[seq_along(MODELS)], labs[ord])

  ggplot(out, aes(fpr, tpr, colour = lab)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3,
                colour = pal[["slate"]]) +
    geom_line(linewidth = 0.35) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_x_continuous(breaks = c(0, 0.5, 1), expand = c(0.01, 0.01)) +
    scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0.01, 0.01)) +
    guides(colour = guide_legend(ncol = 2, override.aes = list(linewidth = 0.5))) +
    labs(x = "1 - specificity", y = "Sensitivity", subtitle = sub) +
    theme_ajr() +
    theme(aspect.ratio = 1,
          legend.position = "right", legend.key.size = unit(2.4, "mm"),
          legend.spacing.y = unit(0, "mm"),
          legend.text = element_text(size = PT_BODY),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

panel_c <- function() {
  a <- read.csv(ML("GSE76293_OOF_predictions.csv"), stringsAsFactors = FALSE)
  b <- read.csv(ML("GSE69528_OOF_predictions.csv"), stringsAsFactors = FALSE)
  d <- rbind(data.frame(cohort = "ARDS discovery", group = a$group, p = a$RF),
             data.frame(cohort = "Sepsis discovery", group = b$group, p = b$RF))
  d$grp <- factor(ifelse(d$group == "Control", "Control", "Disease"),
                  levels = c("Control", "Disease"))

  ggplot(d, aes(p, grp, fill = grp)) +
    geom_violin(width = 0.7, linewidth = 0.25, colour = "#4A4A4A", scale = "width") +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.25) +
    facet_wrap(~ cohort, nrow = 1) +
    scale_fill_manual(values = c(Control = pal[["sage"]], Disease = pal[["brick"]]),
                      guide = "none") +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    labs(x = "Out-of-fold predicted probability (random forest)", y = NULL) +
    theme_ajr() +
    theme(strip.text = element_text(size = PT_BODY, hjust = 0),
          panel.spacing.x = unit(6, "mm"))
}

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid,
                c(list(curve_panel("GSE76293_OOF_predictions.csv", "ARDS discovery (GSE76293)"),
                       curve_panel("GSE69528_OOF_predictions.csv", "Sepsis discovery (GSE69528)"),
                       ncol = 2, labels = c("A", "B")), lab_args))
row2 <- do.call(cowplot::plot_grid, c(list(panel_c(), labels = "C"), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(1.02, 0.98))
save_fig(fig, "FigureS5", dir = out_dir("supplementary"), height_mm = 200)
msg("FigureS5 written")
