# ---------------------------------------------------------------------------
# 30_supp_S10.R -- Figure S10. ssGSEA immune infiltration across four cohorts
# Input  : data/frozen/07_hub_immune/ssgsea_immune/<cohort>_ssgsea_wilcox.csv
# Output : Figures/supplementary/FigureS10.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

HI <- function(...) file.path("data", "frozen", "07_hub_immune", ...)

coh <- c(GSE76293 = "ARDS discovery", GSE69528 = "Sepsis discovery",
         GSE243066 = "ARDS validation", GSE100159 = "Sepsis validation")
# the four cohorts are summarised with the same statistic: the difference in
# mean ssGSEA score between disease and control samples, as written by the
# per-cohort Wilcoxon tables
d <- do.call(rbind, lapply(names(coh), function(cd) {
  x <- read.csv(HI("ssgsea_immune", paste0(cd, "_ssgsea_wilcox.csv")),
                stringsAsFactors = FALSE)
  x$cohort <- coh[[cd]]
  x
}))
d <- d[!is.na(d$delta), ]
d$sig <- ifelse(d$padj < 0.05, "*", "")
ord <- names(sort(tapply(d$delta, d$cell, function(v) mean(abs(v)))))
d$cell <- factor(d$cell, levels = ord)
d$cohort <- factor(d$cohort, levels = rev(coh))

panel_a <- function() {
  ggplot(d, aes(cell, cohort, fill = delta)) +
    geom_tile(colour = "white", linewidth = 0.45) +
    geom_text(aes(label = sig), family = FONT_BODY, size = 9 / .pt,
              colour = pal[["ink"]], vjust = 0.78) +
    scale_fill_gradientn(colours = ramp_div, limits = c(-0.2, 0.2),
                         oob = scales::squish, name = "ssGSEA score\ndifference") +
    scale_x_discrete(labels = function(v) stringr::str_wrap(v, 14)) +
    labs(x = NULL, y = NULL) +
    theme_ajr() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1,
                                     lineheight = 0.85),
          axis.line = element_blank(), axis.ticks = element_blank(),
          legend.key.height = unit(8, "mm"), legend.key.width = unit(2.2, "mm"),
            legend.title = element_text(size = PT_BODY))
}

# --- panel B ---------------------------------------------------------------

panel_b <- function() {
  s <- read.csv(HI("ssgsea_immune", "ssgsea_direction_summary.csv"),
                stringsAsFactors = FALSE)
  s$cohort <- factor(coh[s$cohort], levels = rev(coh))
  d2 <- rbind(data.frame(cohort = s$cohort, dir = "Higher in disease", n = s$n_up),
              data.frame(cohort = s$cohort, dir = "Lower in disease", n = s$n_down))
  d2$dir <- factor(d2$dir, levels = c("Higher in disease", "Lower in disease"))

  ggplot(d2, aes(n, cohort, fill = dir)) +
    geom_col(position = position_dodge(width = 0.68), width = 0.62, colour = NA) +
    geom_text(aes(label = n), position = position_dodge(width = 0.68),
              hjust = -0.35, family = FONT_BODY, size = 7 / .pt,
              colour = pal[["ink"]]) +
    scale_fill_manual(values = c(`Higher in disease` = pal[["brick"]],
                                 `Lower in disease` = pal[["steel"]]), name = NULL) +
    scale_x_continuous(breaks = seq(0, 12, 3),
                       expand = expansion(mult = c(0, 0.14))) +
    labs(x = "Immune cell types (of 19)", y = NULL) +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"))
}

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid, c(list(panel_a(), labels = "A"), lab_args))
row2 <- do.call(cowplot::plot_grid, c(list(panel_b(), labels = "B"), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(1.5, 0.8))
save_fig(fig, "FigureS10", dir = out_dir("supplementary"), height_mm = 200)
msg("FigureS10 written")
