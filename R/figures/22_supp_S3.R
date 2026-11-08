# ---------------------------------------------------------------------------
# 22_supp_S3.R -- Figure S3. CMap glucocorticoid signal against random gene-set nulls
#   A  distribution of the best glucocorticoid P value in DEG-matched null sets
#   B  frequency of a glucocorticoid term in the top ten, by background
# Input  : data/frozen/09_docking_drug/cmap_deg_matched_null.csv,
#          cmap_negative_control_summary.csv
# Output : Figures/supplementary/FigureS3.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

DK <- function(...) file.path("data", "frozen", "09_docking_drug", ...)

nul <- read.csv(DK("cmap_deg_matched_null.csv"), stringsAsFactors = FALSE)
smy <- read.csv(DK("cmap_negative_control_summary.csv"), stringsAsFactors = FALSE)
REAL_P <- 0.0003664321680484968   # hub six-gene set, Old_CMAP_up, dexamethasone

# --- panel A ---------------------------------------------------------------

panel_a <- function() {
  x <- nul[!is.na(nul$best_gluco_p), ]
  br <- sort(unique(x$best_gluco_p))
  ecdf_d <- data.frame(
    p = c(min(br) / 2, br),
    frac = c(0, vapply(br, function(t) mean(x$best_gluco_p <= t), 1))
  )
  frac_hub <- mean(x$best_gluco_p <= REAL_P)

  ggplot(ecdf_d, aes(p, 100 * frac)) +
    geom_step(direction = "hv", linewidth = 0.5, colour = pal[["ink"]]) +
    geom_point(size = 0.9, colour = pal[["ink"]]) +
    geom_vline(xintercept = REAL_P, linetype = "dashed", linewidth = 0.35,
               colour = pal[["brick"]]) +
    annotate("text", x = REAL_P, y = 96, hjust = -0.05,
             label = sprintf("hub six-gene set\nP = 3.7e-04, %.0f%% of nulls", 100 * frac_hub),
             family = FONT_BODY, size = 7 / .pt, colour = pal[["brick"]],
             lineheight = 0.9) +
    scale_x_log10(breaks = c(3e-6, 1e-4, 1e-3, 1e-2, 3e-2),
                  labels = c("3e-06", "1e-04", "1e-03", "1e-02", "3e-02")) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Best glucocorticoid P value per null gene set",
         y = "Null gene sets at or below (%)") +
    theme_ajr()
}

# --- panel B ---------------------------------------------------------------

panel_b <- function() {
  x <- smy[!is.na(smy$gluco_top10_freq_pct), ]
  x <- x[order(x$library, x$gluco_top10_freq_pct), ]
  x$lab <- ifelse(x$background == "real", "Hub six-gene set",
                  gsub("_", " ", sub("^background_", "", x$background)))
  x$lib <- factor(x$library, levels = c("Old_CMAP_up", "Old_CMAP_down"),
                  labels = c("CMap up-regulated", "CMap down-regulated"))
  x$lab <- factor(x$lab, levels = unique(x$lab[order(x$library, x$gluco_top10_freq_pct)]))
  x$is_real <- x$background == "real"

  ggplot(x, aes(gluco_top10_freq_pct, lab, fill = is_real)) +
    geom_col(width = 0.66, colour = NA) +
    geom_text(aes(label = sprintf("%.0f", gluco_top10_freq_pct)), hjust = -0.2,
              family = FONT_BODY, size = 7 / .pt, colour = pal[["ink"]]) +
    facet_wrap(~ lib, nrow = 1) +
    scale_fill_manual(values = c(`TRUE` = pal[["brick"]], `FALSE` = pal[["slate"]]),
                      guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(x = "Null sets with a glucocorticoid in the top ten (%)", y = NULL) +
    theme_ajr() +
    theme(strip.text = element_text(size = PT_BODY, hjust = 0),
          panel.spacing.x = unit(5, "mm"))
}

# --- assemble --------------------------------------------------------------

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid, c(list(panel_a(), labels = "A"), lab_args))
row2 <- do.call(cowplot::plot_grid, c(list(panel_b(), labels = "B"), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1.15))
save_fig(fig, "FigureS3", dir = out_dir("supplementary"), height_mm = 200)
msg("FigureS3 written")
