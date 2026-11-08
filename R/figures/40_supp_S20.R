# ---------------------------------------------------------------------------
# 40_supp_S20.R -- Figure S20. Six-hub-gene prognostic performance in GSE65682
#   A  survival by risk group
#   B  log-rank P value of each hub gene
#   C  time-dependent AUC of the risk score
#   D  calibration of the clinical and combined models
#   E  risk score adjusted for age and sex
# Input  : data/frozen/13_survival_prognosis/tables
# Output : Figures/supplementary/FigureS20.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")
suppressPackageStartupMessages(library(survival))

SP <- function(...) file.path("data", "frozen", "13_survival_prognosis", "tables", ...)

risk <- read.csv(SP("GSE65682_risk_scores.csv"), stringsAsFactors = FALSE)
risk$riskGroup <- factor(risk$riskGroup, levels = c("Low", "High"))

# --- panel A ---------------------------------------------------------------

panel_a <- function() {
  fit <- survfit(Surv(tte, event) ~ riskGroup, data = risk)
  sd <- survdiff(Surv(tte, event) ~ riskGroup, data = risk)
  p <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  s <- summary(fit)
  d <- data.frame(time = s$time, surv = s$surv,
                  grp = sub("riskGroup=", "", as.character(s$strata)))
  d <- rbind(data.frame(time = 0, surv = 1, grp = unique(d$grp)), d)

  ggplot(d, aes(time, surv, colour = grp)) +
    geom_step(linewidth = 0.5) +
    annotate("text", x = 1, y = 0.08, hjust = 0,
             label = sprintf("log-rank P = %s", formatC(p, format = "e", digits = 1)),
             family = FONT_BODY, size = 7 / .pt, colour = pal[["ink"]]) +
    scale_colour_manual(values = c(Low = pal[["sage"]], High = pal[["brick"]]),
                        name = NULL) +
    scale_x_continuous(breaks = seq(0, 28, 7), expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       expand = c(0.01, 0.01)) +
    labs(x = "Days from admission", y = "Survival probability",
         subtitle = "28-day survival by risk group") +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

# --- panel B ---------------------------------------------------------------

panel_b <- function() {
  d <- read.csv(SP("pergene_KM_logrank.csv"), stringsAsFactors = FALSE)
  d <- d[d$gene %in% HUB6, ]
  d$gene <- factor(d$gene, levels = rev(HUB6))
  d$nlp <- -log10(d$logrank_P)

  ggplot(d, aes(nlp, gene)) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    geom_col(width = 0.6, fill = pal[["slate"]], colour = NA) +
    geom_text(aes(label = sprintf("%.2f", logrank_P)), hjust = -0.2,
              family = FONT_BODY, size = 7 / .pt, colour = pal[["ink"]]) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.26))) +
    labs(x = "-log10 log-rank P", y = NULL, subtitle = "Single-gene survival") +
    theme_ajr() +
    theme(axis.text.y = element_text(face = "italic"),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

# --- panel C ---------------------------------------------------------------

panel_c <- function() {
  d <- read.csv(SP("timeROC_riskScore.csv"), stringsAsFactors = FALSE)
  ci <- do.call(rbind, strsplit(as.character(d$CI), "-"))
  d$lo <- as.numeric(ci[, 1]); d$hi <- as.numeric(ci[, 2])

  ggplot(d, aes(time, AUC)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 1.2, linewidth = 0.35,
                  colour = pal[["steel"]]) +
    geom_point(size = 1.3, colour = pal[["steel"]]) +
    scale_x_continuous(breaks = c(7, 14, 28)) +
    scale_y_continuous(limits = c(0.5, 0.75), breaks = seq(0.5, 0.75, 0.05)) +
    labs(x = "Days", y = "Time-dependent AUC", subtitle = "Risk-score discrimination") +
    theme_ajr() +
    theme(plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

# --- panel D ---------------------------------------------------------------

panel_d <- function() {
  d <- read.csv(SP("ext_calibration.csv"), stringsAsFactors = FALSE)
  ggplot(d, aes(pred, obs, colour = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3,
                colour = pal[["slate"]]) +
    geom_line(linewidth = 0.4) +
    geom_point(size = 1.2) +
    scale_colour_manual(values = c(Clinical = pal[["slate"]],
                                   `Clinical + Hub6` = pal[["brick"]]), name = NULL) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Predicted 28-day mortality", y = "Observed mortality",
         subtitle = "Calibration") +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

# --- panel E ---------------------------------------------------------------

panel_e <- function() {
  d <- read.csv(SP("adjusted_cox_forest.csv"), stringsAsFactors = FALSE)
  d$var <- factor(d$var, levels = rev(d$var))
  d$lab <- sprintf("%.2f (%.2f-%.2f)", d$HR, d$lo, d$hi)

  ggplot(d, aes(HR, var)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.4,
                   colour = pal[["ink"]]) +
    geom_point(shape = 23, size = 1.5, fill = pal[["ink"]], colour = pal[["ink"]]) +
    geom_text(aes(x = hi, label = lab), hjust = -0.12, family = FONT_BODY,
              size = 7 / .pt, colour = pal[["ink"]]) +
    scale_x_log10(breaks = c(0.5, 1, 2, 5),
                  expand = expansion(mult = c(0.10, 0.42))) +
    labs(x = "Hazard ratio (95% CI)", y = NULL, subtitle = "Adjusted Cox model") +
    theme_ajr() +
    theme(plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid,
                c(list(panel_a(), panel_b(), ncol = 2, rel_widths = c(0.55, 0.45),
                       labels = c("A", "B")), lab_args))
row2 <- do.call(cowplot::plot_grid,
                c(list(panel_c(), panel_d(), panel_e(), ncol = 3,
                       labels = c("C", "D", "E")), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(1.05, 1))
save_fig(fig, "FigureS20", dir = out_dir("supplementary"), height_mm = 210)
msg("FigureS20 written")
