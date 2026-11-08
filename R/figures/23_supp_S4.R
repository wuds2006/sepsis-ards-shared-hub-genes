# ---------------------------------------------------------------------------
# 23_supp_S4.R -- Figure S4. Redocking validation of the docking protocol
#   A  per-mode heavy-atom RMSD against affinity for the two co-crystal cases
#   B  top-ranked versus best affinity of the redocked pose
# Input  : data/frozen/09_docking_drug/redocking
# Output : Figures/supplementary/FigureS4.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

DK <- function(...) file.path("data", "frozen", "09_docking_drug", "redocking", ...)

mode <- read.csv(DK("redocking_per_mode.csv"), stringsAsFactors = FALSE)
top <- read.csv(DK("redocking_results.csv"), stringsAsFactors = FALSE)

case_lab <- c(UPP1 = "UPP1 (3EUF, BAU)", VNN1 = "VNN1 (7SLV, 9TC)")
mode$top_ranked <- as.logical(mode$top_ranked)
top$target <- sub("_.*$", "", top$case)
mode$case <- factor(case_lab[mode$target], levels = case_lab)
top$case <- factor(case_lab[top$target], levels = case_lab)

# --- panel A ---------------------------------------------------------------

panel_a <- function() {
  ggplot(mode, aes(rmsd_A, affinity_kcal_mol)) +
    geom_vline(xintercept = 2, linetype = "dashed", linewidth = 0.3,
               colour = pal[["brick"]]) +
    geom_point(aes(fill = top_ranked), shape = 21, size = 1.6, stroke = 0.3,
               colour = pal[["ink"]]) +
    geom_text(aes(label = mode), vjust = -1.0, family = FONT_BODY,
              size = 7 / .pt, colour = pal[["ink"]]) +
    facet_wrap(~ case, nrow = 1) +
    scale_fill_manual(values = c(`TRUE` = pal[["brick"]], `FALSE` = "white"),
                      labels = c(`TRUE` = "top-ranked pose", `FALSE` = "other poses"),
                      name = NULL) +
    scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
    labs(x = "Heavy-atom RMSD to the co-crystal ligand (A)",
         y = "Predicted affinity (kcal/mol)") +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"),
          strip.text = element_text(size = PT_BODY, hjust = 0),
          panel.spacing.x = unit(6, "mm"))
}

# --- panel B ---------------------------------------------------------------

panel_b <- function() {
  d <- rbind(
    data.frame(case = top$case, what = "top-ranked pose",
               aff = top$top1_affinity_kcal_mol),
    data.frame(case = top$case, what = "best of all poses",
               aff = top$best_affinity_kcal_mol)
  )
  d$what <- factor(d$what, levels = c("best of all poses", "top-ranked pose"))
  d$lab <- sprintf("%.2f", d$aff)

  ggplot(d, aes(aff, case)) +
    geom_line(aes(group = case), linewidth = 0.4, colour = "#9AA0A6") +
    geom_point(aes(fill = what), shape = 21, size = 1.8, stroke = 0.35,
               colour = pal[["ink"]]) +
    geom_text(aes(label = lab), hjust = -0.35, family = FONT_BODY,
              size = 7 / .pt, colour = pal[["ink"]]) +
    scale_fill_manual(values = c(`best of all poses` = palette()[1],
                                 `top-ranked pose` = pal[["brick"]]), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.12))) +
    labs(x = "Predicted affinity (kcal/mol)", y = NULL) +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"))
}

# --- assemble --------------------------------------------------------------

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid, c(list(panel_a(), labels = "A"), lab_args))
row2 <- do.call(cowplot::plot_grid, c(list(panel_b(), labels = "B"), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(1.25, 0.75))
save_fig(fig, "FigureS4", dir = out_dir("supplementary"), height_mm = 200)
msg("FigureS4 written")
