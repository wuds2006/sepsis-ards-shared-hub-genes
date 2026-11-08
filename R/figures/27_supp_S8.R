# ---------------------------------------------------------------------------
# 27_supp_S8.R -- Figure S8. Sensitivity analyses for the null Mendelian
# randomisation estimates: pleiotropy, heterogeneity and direction of effect
#   A  forward analysis (sepsis to ARDS)
#   B  reverse analysis (ARDS to sepsis)
#   C  Cochran Q by estimator
# Input  : data/frozen/04_MR_SMR
# Output : Figures/supplementary/FigureS8.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

MR <- function(...) file.path("data", "frozen", "04_MR_SMR", ...)

read_dir <- function(prefix) {
  pl <- read.csv(MR(sprintf("MR_%spleiotropy.csv", prefix)), stringsAsFactors = FALSE)
  he <- read.csv(MR(sprintf("MR_%sheterogeneity.csv", prefix)), stringsAsFactors = FALSE)
  st <- read.csv(MR(sprintf("MR_%ssteiger.csv", prefix)), stringsAsFactors = FALSE)
  ivw <- he[he$method == "Inverse variance weighted", ][1, ]
  list(pl = pl[1, ], he = he, st = st[1, ], ivw = ivw)
}

fw <- read_dir("sepsis4980_")
rv <- read_dir("reverse_")

sens_panel <- function(x, sub) {
  d <- data.frame(
    test = factor(c("Pleiotropy (Egger intercept)", "Heterogeneity (Cochran Q)",
                    "Direction (Steiger)"),
                  levels = rev(c("Pleiotropy (Egger intercept)", "Heterogeneity (Cochran Q)",
                                 "Direction (Steiger)"))),
    p = c(x$pl$pval, x$ivw$Q_pval, x$st$steiger_pval)
  )
  d$nlp <- -log10(d$p)
  d$lab <- formatC(d$p, format = "e", digits = 1)

  ggplot(d, aes(nlp, test)) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    geom_segment(aes(x = 0, xend = nlp, yend = test), linewidth = 0.4,
                 colour = "#C4C9CE") +
    geom_point(size = 1.6, colour = pal[["steel"]]) +
    geom_text(aes(label = sprintf("P = %s", lab)), hjust = -0.25,
              family = FONT_BODY, size = 7 / .pt, colour = pal[["ink"]]) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.55))) +
    labs(x = "-log10 P", y = NULL, subtitle = sub) +
    theme_ajr() +
    theme(plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

# --- panel C ---------------------------------------------------------------

panel_c <- function() {
  d <- rbind(
    data.frame(dir = "Forward", fw$he[, c("method", "Q", "Q_pval")]),
    data.frame(dir = "Reverse", rv$he[, c("method", "Q", "Q_pval")])
  )
  d$dir <- factor(d$dir, levels = c("Forward", "Reverse"))
  d$lab <- ifelse(d$Q_pval < 0.001, "<0.001", sprintf("%.3f", d$Q_pval))

  ggplot(d, aes(Q, method, fill = dir)) +
    geom_col(position = position_dodge(width = 0.66), width = 0.6, colour = NA) +
    geom_text(aes(label = sprintf("P = %s", lab)),
              position = position_dodge(width = 0.66), hjust = -0.18,
              family = FONT_BODY, size = 7 / .pt, colour = pal[["ink"]]) +
    scale_fill_manual(values = c(Forward = pal[["brick"]], Reverse = pal[["steel"]]),
                      name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(x = "Cochran Q", y = NULL) +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"))
}

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid,
                c(list(sens_panel(fw, "Forward: sepsis to ARDS"),
                       sens_panel(rv, "Reverse: ARDS to sepsis"),
                       ncol = 2, labels = c("A", "B")), lab_args))
row2 <- do.call(cowplot::plot_grid, c(list(panel_c(), labels = "C"), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 0.7))
save_fig(fig, "FigureS8", dir = out_dir("supplementary"), height_mm = 200)
msg("FigureS8 written")
