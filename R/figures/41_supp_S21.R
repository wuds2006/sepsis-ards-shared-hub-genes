# ---------------------------------------------------------------------------
# 41_supp_S21.R -- Figure S21. Extended prognostic modelling over the 179
# shared candidate genes in GSE65682
#   A  model discrimination
#   B  calibration
#   C  decision-curve analysis
#   D  reclassification and discrimination improvement
#   E  penalised regression coefficients
#   F  random-survival-forest importance
# Input  : data/frozen/13_survival_prognosis/tables
# Output : Figures/supplementary/FigureS21.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

SP <- function(...) file.path("data", "frozen", "13_survival_prognosis", "tables", ...)

# model names differ between the performance and the decision-curve tables, so
# every label that can appear in either table is mapped to its colour up front
SHORT <- c("Clinical" = "Clinical",
           "Hub6" = "Hub6",
           "LASSO-179" = "LASSO-179",
           "Clinical + Hub6" = "Clinical + Hub6",
           "Clinical + LASSO-179" = "Clinical + LASSO",
           "RSF" = "RSF",
           "RSF (OOB)" = "RSF",
           "Treat all" = "Treat all",
           "Treat none" = "Treat none")
LBL_COL <- c("Clinical" = unname(pal[["slate"]]),
             "Hub6" = unname(pal[["brick"]]),
             "LASSO-179" = unname(pal[["olive"]]),
             "Clinical + Hub6" = unname(pal[["ochre"]]),
             "Clinical + LASSO" = unname(pal[["sage"]]),
             "RSF" = unname(pal[["plum"]]),
             "Treat all" = "#9AA0A6",
             "Treat none" = "#C4C9CE")

relabel <- function(x) {
  lab <- ifelse(x %in% names(SHORT), SHORT[x], x)
  factor(lab, levels = unique(lab))
}

lbl_values <- function(lab) {
  v <- LBL_COL[levels(lab)]
  if (anyNA(v)) stop("no colour assigned to: ", paste(levels(lab)[is.na(v)], collapse = ", "))
  v
}

# --- panel A ---------------------------------------------------------------

panel_a <- function() {
  d <- read.csv(SP("ext_model_performance.csv"), stringsAsFactors = FALSE)
  d <- d[order(d$C_apparent), ]
  d$model <- factor(d$model, levels = d$model)

  ggplot(d, aes(C_apparent, model)) +
    geom_col(width = 0.62, fill = pal[["steel"]], colour = NA) +
    geom_text(aes(label = sprintf("%.3f", C_apparent)), hjust = -0.2,
              family = FONT_BODY, size = 7 / .pt, colour = pal[["ink"]]) +
    scale_x_continuous(limits = c(0, 0.88), breaks = seq(0, 0.6, 0.2),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = "C-index", y = NULL, subtitle = "Discrimination") +
    theme_ajr() +
    theme(plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

# --- panel B ---------------------------------------------------------------

panel_b <- function() {
  d <- read.csv(SP("ext_calibration.csv"), stringsAsFactors = FALSE)
  d$lab <- relabel(d$model)

  ggplot(d, aes(pred, obs, colour = lab)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3,
                colour = pal[["slate"]]) +
    geom_line(linewidth = 0.35) +
    geom_point(size = 1.0) +
    scale_colour_manual(values = lbl_values(d$lab), name = NULL) +
    guides(colour = guide_legend(ncol = 2)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Predicted 28-day mortality", y = "Observed", subtitle = "Calibration") +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0),
          legend.key.size = unit(2.4, "mm"),
          legend.text = element_text(size = PT_BODY),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

# --- panel C ---------------------------------------------------------------

panel_c <- function() {
  d <- read.csv(SP("ext_dca.csv"), stringsAsFactors = FALSE)
  d$lab <- relabel(d$model)
  lt <- ifelse(levels(d$lab) %in% c("Treat all", "Treat none"), "dashed", "solid")

  ggplot(d, aes(threshold, net_benefit, colour = lab)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = pal[["slate"]]) +
    geom_line(aes(linetype = lab), linewidth = 0.4) +
    scale_colour_manual(values = lbl_values(d$lab), name = NULL) +
    scale_linetype_manual(values = setNames(lt, levels(d$lab)), guide = "none") +
    guides(colour = guide_legend(ncol = 2)) +
    scale_x_continuous(limits = c(0.05, 0.6), breaks = seq(0.1, 0.6, 0.1)) +
    labs(x = "Threshold probability", y = "Net benefit", subtitle = "Decision curves") +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0),
          legend.key.size = unit(2.4, "mm"),
          legend.text = element_text(size = PT_BODY),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

# --- panel D ---------------------------------------------------------------

panel_d <- function() {
  d <- read.csv(SP("ext_nri_idi.csv"), stringsAsFactors = FALSE)
  d$lab <- sub(" vs Clinical", "", d$comparison)
  d$star <- ifelse(d$NRI_P < 0.05, "*", "")
  m <- rbind(data.frame(model = d$lab, metric = "NRI", value = d$NRI, star = d$star),
             data.frame(model = d$lab, metric = "IDI", value = d$IDI, star = d$star))

  ggplot(m, aes(value, model, fill = metric)) +
    geom_col(position = position_dodge(width = 0.68), width = 0.6, colour = NA) +
    geom_text(aes(label = star), position = position_dodge(width = 0.68),
              hjust = -0.5, family = FONT_BODY, size = 7 / .pt,
              colour = pal[["ink"]]) +
    scale_fill_manual(values = c(NRI = pal[["steel"]], IDI = pal[["ochre"]]), name = NULL) +
    scale_x_continuous(limits = c(0, 0.185), breaks = seq(0, 0.16, 0.04),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "Improvement vs clinical", y = NULL,
         subtitle = "Reclassification") +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

# --- panel E ---------------------------------------------------------------

panel_e <- function() {
  d <- read.csv(SP("ext_lasso_coefficients.csv"), stringsAsFactors = FALSE)
  d <- d[order(-abs(d$beta_standardised)), ][1:15, ]
  d$gene <- factor(d$gene, levels = rev(d$gene))
  d$dir <- ifelse(d$beta_standardised > 0, "Higher risk", "Lower risk")

  ggplot(d, aes(beta_standardised, gene, fill = dir)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = pal[["ink"]]) +
    geom_col(width = 0.62, colour = NA) +
    scale_fill_manual(values = c(`Higher risk` = pal[["brick"]],
                                 `Lower risk` = pal[["steel"]]), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
    labs(x = "Standardised coefficient", y = NULL, subtitle = "Penalised regression") +
    theme_ajr() +
    theme(axis.text.y = element_text(face = "italic"),
          legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

# --- panel F ---------------------------------------------------------------

panel_f <- function() {
  d <- read.csv(SP("ext_rsf_importance.csv"), stringsAsFactors = FALSE)
  d$var <- sub("_imp$", "", d$variable)
  d <- d[order(-d$importance), ][1:15, ]
  d$var <- ifelse(d$var == "age", "Age", d$var)
  d$var <- factor(d$var, levels = rev(d$var))

  ggplot(d, aes(importance, var)) +
    geom_col(width = 0.62, fill = pal[["plum"]], colour = NA) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(x = "Variable importance", y = NULL, subtitle = "Random survival forest") +
    theme_ajr() +
    theme(axis.text.y = element_text(face = "italic"),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid,
                c(list(panel_a(), panel_b(), panel_c(), ncol = 3,
                       labels = c("A", "B", "C")), lab_args))
row2 <- do.call(cowplot::plot_grid,
                c(list(panel_d(), panel_e(), panel_f(), ncol = 3,
                       labels = c("D", "E", "F")), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1)
save_fig(fig, "FigureS21", dir = out_dir("supplementary"), height_mm = 208)
msg("FigureS21 written")
