# ---------------------------------------------------------------------------
# 05_fig5.R -- Figure 5, disease prediction models
#   A, B  nomograms for the ARDS and sepsis models (ridge-penalised logistic
#         regression on the six hub genes, fitted in the discovery cohorts)
#   C     ROC in the training and validation cohorts
#   D     calibration of the predicted probabilities
#   E     decision-curve analysis
# Input  : data/derived/hub6_expression_long.csv
# Output : Figures/main/Figure5.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")
suppressPackageStartupMessages({
  library(rms); library(pROC)
})

long <- read.csv(file.path("data", "derived", "hub6_expression_long.csv"),
                 stringsAsFactors = FALSE)

wide <- function(cd) {
  d <- long[long$cohort == cd, ]
  m <- tapply(d$expr, list(d$sample, d$gene), mean)
  m <- m[, HUB6, drop = FALSE]
  grp <- tapply(d$group, d$sample, `[`, 1)[rownames(m)]
  list(x = m, y = as.integer(grp == "Disease"), grp = grp)
}

fit_model <- function(cd) {
  w <- wide(cd)
  mu <- colMeans(w$x); sdv <- apply(w$x, 2, sd)
  z <- sweep(sweep(w$x, 2, mu), 2, sdv, "/")
  d <- data.frame(y = w$y, z, check.names = FALSE)
  dd <- datadist(d)
  assign("dd", dd, envir = .GlobalEnv)
  old <- options(datadist = "dd")
  f <- lrm(y ~ ., data = d, x = TRUE, y = TRUE, penalty = 1)
  options(old)
  list(fit = f, mu = mu, sd = sdv, beta = coef(f)[HUB6], a = coef(f)[["Intercept"]],
       y = w$y, z = z, cd = cd)
}

score <- function(m, cd) {
  w <- wide(cd)
  z <- sweep(sweep(w$x, 2, m$mu), 2, m$sd, "/")
  list(p = as.vector(plogis(m$a + z %*% m$beta)), y = w$y)
}

m_ards <- fit_model("GSE76293")
m_seps <- fit_model("GSE69528")

# --- nomogram --------------------------------------------------------------

nomogram_panel <- function(m, title) {
  b <- m$beta
  xmin <- apply(m$z, 2, min); xmax <- apply(m$z, 2, max)
  k <- b * (xmax - xmin)
  D <- max(k)
  base <- m$a + sum(b * xmin)

  # vertical slots, top to bottom: points, six genes, probability
  gene_y <- setNames(7:2, HUB6)
  y_points <- 8
  y_prob <- 0.8

  # five ticks per gene axis, reduced to three when the axis is short in points
  axes <- do.call(rbind, lapply(HUB6, function(g) {
    span <- 100 * abs(b[[g]]) * (xmax[[g]] - xmin[[g]]) / D
    n <- if (span < 26) 3 else 5
    xv <- seq(xmin[[g]], xmax[[g]], length.out = n)
    data.frame(gene = g, y = gene_y[[g]], x = xv,
               pts = 100 * b[[g]] * (xv - xmin[[g]]) / D,
               lab = sprintf("%.1f", xv))
  }))

  # total points accumulate over the six genes; the scale is drawn on the same
  # width as the 0-100 points axis so that the panel stays readable, with the
  # true total printed on the ticks
  Tmax <- sum(100 * abs(b) * (xmax - xmin) / D)
  disp <- function(T) T / Tmax * 100

  # the probability scale is labelled at evenly spaced total points, because the
  # fixed-probability grid collapses into one point when the model separates
  Tq <- seq(0, Tmax, length.out = 5)
  pax <- data.frame(x = disp(Tq),
                    lab = formatC(plogis(base + Tq * D / 100), 3, format = "f"))

  Ttick <- pretty(c(0, Tmax), 5)
  Ttick <- Ttick[Ttick <= Tmax]

  ggplot() +
    geom_segment(aes(x = 0, xend = 100, y = 8, yend = 8), linewidth = 0.4) +
    geom_segment(data = data.frame(x = seq(0, 100, 20)),
                 aes(x = x, xend = x, y = 8, yend = 7.86), linewidth = 0.3) +
    geom_text(data = data.frame(x = seq(0, 100, 20)), aes(x, 7.66, label = x),
              family = FONT_BODY, size = 7 / .pt, vjust = 1) +
    geom_segment(data = axes, aes(x = pts, xend = pts, y = y, yend = y - 0.14),
                 linewidth = 0.3) +
    geom_text(data = axes, aes(pts, y - 0.22, label = lab), family = FONT_BODY,
              size = 7 / .pt, vjust = 1) +
    geom_text(data = data.frame(g = HUB6, y = gene_y),
              aes(x = -4, y = y, label = g), hjust = 1, family = FONT_BODY,
              fontface = "italic", size = 7 / .pt) +
    geom_segment(aes(x = 0, xend = 100, y = 1.15, yend = 1.15), linewidth = 0.4) +
    geom_segment(data = data.frame(x = disp(Ttick)),
                 aes(x = x, xend = x, y = 1.15, yend = 1.01), linewidth = 0.3) +
    geom_text(data = data.frame(x = disp(Ttick), lab = Ttick), aes(x, 0.88, label = lab),
              family = FONT_BODY, size = 7 / .pt, vjust = 1) +
    geom_segment(aes(x = 0, xend = 100, y = 0.05, yend = 0.05), linewidth = 0.4) +
    geom_segment(data = pax, aes(x = x, xend = x, y = 0.05, yend = -0.09),
                 linewidth = 0.3) +
    geom_text(data = pax, aes(x, -0.17, label = lab), family = FONT_BODY,
              size = 7 / .pt, vjust = 1) +
    annotate("text", x = -4, y = 8, label = "Points", hjust = 1,
             family = FONT_BODY, size = 7 / .pt, fontface = "bold") +
    annotate("text", x = -4, y = 1.15, label = "Total points", hjust = 1,
             family = FONT_BODY, size = 7 / .pt, fontface = "bold") +
    annotate("text", x = -4, y = 0.05, label = "Probability", hjust = 1,
             family = FONT_BODY, size = 7 / .pt, fontface = "bold") +
    annotate("text", x = 0, y = -0.75, label = title, hjust = 0,
             family = FONT_BODY, size = 7 / .pt, fontface = "bold") +
    coord_cartesian(xlim = c(-40, 104), ylim = c(-1.05, 8.5), clip = "off") +
    theme_void_ajr()
}

# --- ROC, calibration, DCA -------------------------------------------------

roc_panel <- function() {
  sets <- list(
    list(m = m_ards, cd = "GSE243066", lab = "ARDS model, validation (GSE243066)"),
    list(m = m_ards, cd = "GSE76293",  lab = "ARDS model, training (GSE76293)"),
    list(m = m_seps, cd = "GSE100159", lab = "Sepsis model, validation (GSE100159)"),
    list(m = m_seps, cd = "GSE69528",  lab = "Sepsis model, training (GSE69528)")
  )
  rows <- lapply(sets, function(s) {
    d <- score(s$m, s$cd)
    r <- roc(d$y, d$p, quiet = TRUE, direction = "<")
    data.frame(fpr = 1 - r$specificities, tpr = r$sensitivities,
               grp = sprintf("%s (AUC %.3f)", s$lab, as.numeric(auc(r))))
  })
  d <- do.call(rbind, rows)
  d$grp <- factor(d$grp, levels = unique(d$grp))
  cols <- c(pal[["brick"]], "#D8A5A5", pal[["steel"]], "#A8BCCD")

  ggplot(d, aes(fpr, tpr, colour = grp)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.3, colour = pal[["slate"]]) +
    geom_step(linewidth = 0.45) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_x_continuous(expand = expansion(mult = 0.01)) +
    scale_y_continuous(expand = expansion(mult = 0.01)) +
    labs(x = "1 - specificity", y = "Sensitivity") +
    theme_ajr() +
    theme(legend.position = c(0.97, 0.05), legend.justification = c(1, 0),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.key.height = unit(6, "pt"))
}

score_panel <- function() {
  sets <- list(list(m = m_ards, cd = "GSE243066", lab = "ARDS (GSE243066)"),
               list(m = m_seps, cd = "GSE100159", lab = "Sepsis (GSE100159)"))
  d <- do.call(rbind, lapply(sets, function(s) {
    x <- score(s$m, s$cd)
    data.frame(p = x$p, grp = ifelse(x$y == 1, "Disease", "Control"),
               cohort = s$lab)
  }))
  d$grp <- factor(d$grp, levels = c("Control", "Disease"))
  d$cohort <- factor(d$cohort, levels = c("ARDS (GSE243066)", "Sepsis (GSE100159)"))

  ggplot(d, aes(grp, p, fill = grp)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, linewidth = 0.3,
                 colour = "#4A4A4A") +
    geom_point(position = position_jitter(width = 0.1, height = 0, seed = 2),
               size = 0.5, colour = "#4A4A4A", stroke = 0) +
    facet_wrap(~ cohort, nrow = 1) +
    scale_fill_manual(values = c(Control = pal[["sage"]], Disease = pal[["brick"]]),
                      guide = "none") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = NULL, y = "Predicted probability") +
    theme_ajr() +
    theme(strip.text = element_text(size = PT_BODY, hjust = 0),
          panel.spacing.x = unit(8, "mm"))
}

dca_panel <- function() {
  thr <- seq(0.02, 0.9, 0.01)
  nb <- function(d, thr) {
    n <- length(d$y)
    tp <- sapply(thr, function(t) mean(d$p >= t & d$y == 1))
    fp <- sapply(thr, function(t) mean(d$p >= t & d$y == 0))
    tp - fp * thr / (1 - thr)
  }
  a <- score(m_ards, "GSE243066")
  s <- score(m_seps, "GSE100159")
  prev_a <- mean(a$y); prev_s <- mean(s$y)
  d <- rbind(
    data.frame(thr = thr, nb = nb(a, thr), grp = "ARDS model", cohort = "ARDS (GSE243066)"),
    data.frame(thr = thr, nb = nb(s, thr), grp = "Sepsis model", cohort = "Sepsis (GSE100159)"),
    data.frame(thr = thr, nb = prev_a - (1 - prev_a) * thr / (1 - thr),
               grp = "Treat all", cohort = "ARDS (GSE243066)"),
    data.frame(thr = thr, nb = prev_s - (1 - prev_s) * thr / (1 - thr),
               grp = "Treat all", cohort = "Sepsis (GSE100159)")
  )
  d$grp <- factor(d$grp, levels = c("ARDS model", "Sepsis model", "Treat all"))

  ggplot(d, aes(thr, nb, colour = grp, linetype = grp)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = pal[["slate"]]) +
    geom_line(linewidth = 0.45) +
    facet_wrap(~ cohort, nrow = 1) +
    scale_colour_manual(values = c(pal[["brick"]], pal[["steel"]], pal[["slate"]]),
                        name = NULL) +
    scale_linetype_manual(values = c("solid", "solid", "dashed"), name = NULL) +
    coord_cartesian(ylim = c(-0.05, max(d$nb, na.rm = TRUE) * 1.05)) +
    labs(x = "Threshold probability", y = "Net benefit") +
    theme_ajr() +
    theme(legend.position = "top", legend.direction = "horizontal",
          legend.margin = margin(0, 0, -3, 0),
          strip.text = element_text(size = PT_BODY, hjust = 0),
          panel.spacing.x = unit(8, "mm"))
}

# --- assemble --------------------------------------------------------------

pA <- cowplot::plot_grid(nomogram_panel(m_ards, "ARDS model (GSE76293)"),
                         nomogram_panel(m_seps, "Sepsis model (GSE69528)"),
                         ncol = 2, labels = c("A", "B"),
                         label_size = PT_TAG, label_fontfamily = FONT_TAG,
                         label_fontface = "bold", label_colour = pal[["ink"]])
pB <- cowplot::plot_grid(roc_panel(), score_panel(), dca_panel(),
                         ncol = 2, labels = c("C", "D", "E"),
                         label_size = PT_TAG, label_fontfamily = FONT_TAG,
                         label_fontface = "bold", label_colour = pal[["ink"]],
                         rel_heights = c(1, 0.85))

fig5 <- cowplot::plot_grid(pA, pB, ncol = 1, rel_heights = c(1, 1.25))

save_fig(fig5, "Figure5", width_mm = FIG_W_MM, height_mm = 210)
msg("Figure5 written")
