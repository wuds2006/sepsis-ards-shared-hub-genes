# ---------------------------------------------------------------------------
# 44_supp_S24.R -- Figure S24. CIBERSORT immune-cell composition of the four
# bulk cohorts (the stacked-bar view retained for completeness)
# Input  : data/frozen/07_hub_immune/cibersort/*_CIBERSORT_results.csv
# Output : Figures/supplementary/FigureS24.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

CI <- function(...) file.path("data", "frozen", "07_hub_immune", "cibersort", ...)

coh <- c(GSE76293 = "ARDS discovery (GSE76293, n = 24)",
         GSE69528 = "Sepsis discovery (GSE69528, n = 111)",
         GSE243066 = "ARDS validation (GSE243066, n = 49)",
         GSE100159 = "Sepsis validation (GSE100159, n = 45)")

sub <- c("B cells naive", "B cells memory", "Plasma cells", "T cells CD8",
         "T cells CD4 naive", "T cells CD4 memory resting",
         "T cells CD4 memory activated", "T cells follicular helper",
         "T cells regulatory (Tregs)", "T cells gamma delta", "NK cells resting",
         "NK cells activated", "Monocytes", "Macrophages M0", "Macrophages M1",
         "Macrophages M2", "Dendritic cells resting", "Dendritic cells activated",
         "Mast cells resting", "Mast cells activated", "Eosinophils", "Neutrophils")

d <- do.call(rbind, lapply(names(coh), function(cd) {
  x <- read.csv(CI(paste0(cd, "_CIBERSORT_results.csv")), stringsAsFactors = FALSE,
                check.names = FALSE)
  x <- x[!is.na(x[[sub[1]]]), ]
  long <- do.call(rbind, lapply(sub, function(s) {
    data.frame(sample = x$sample, group = x$group, cell = s, frac = x[[s]])
  }))
  long$cohort <- coh[[cd]]
  long
}))
d$cohort <- factor(d$cohort, levels = coh)
d$cell <- factor(d$cell, levels = sub)
d$grp <- factor(ifelse(d$group == "Control", "Control", "Disease"),
                levels = c("Control", "Disease"))
ord <- unique(d$sample[order(d$cohort, d$grp, d$sample)])
d$sample <- factor(d$sample, levels = ord)

# 22 low-chroma hues: five families (B/plasma, T, NK, myeloid, granulocyte/mast)
fam_col <- c(rep(c("#8A5F5F", "#B07A6E", "#C99A8B"), 3),      # B, plasma, T CD4
             "#4C6E91", "#6E8CA8", "#9BB0C4",                  # CD8, Tfh, Treg
             "#7E9E95", "#A6BDB4",                             # gamma delta, NK
             "#C99A5B", "#D9BC8C", "#8A7391", "#A994AE",       # mono, macro
             "#7D8C5C", "#A3B084",                             # DC
             "#6F7B85", "#98A2AA", "#B9C0C6",                  # mast
             "#B05B5B", "#C98B7B", "#D6C3AE", "#8C8C8C", "#C4C9CE")
names(fam_col) <- sub

# short display labels: the full LM22 names do not fit a 45 mm legend column
disp <- c("B cells naive" = "B naive",
          "B cells memory" = "B memory",
          "Plasma cells" = "Plasma cells",
          "T cells CD8" = "CD8 T",
          "T cells CD4 naive" = "CD4 T naive",
          "T cells CD4 memory resting" = "CD4 T memory rest.",
          "T cells CD4 memory activated" = "CD4 T memory act.",
          "T cells follicular helper" = "Tfh",
          "T cells regulatory (Tregs)" = "Tregs",
          "T cells gamma delta" = "T gamma-delta",
          "NK cells resting" = "NK resting",
          "NK cells activated" = "NK activated",
          "Monocytes" = "Monocytes",
          "Macrophages M0" = "Mac M0",
          "Macrophages M1" = "Mac M1",
          "Macrophages M2" = "Mac M2",
          "Dendritic cells resting" = "DC resting",
          "Dendritic cells activated" = "DC activated",
          "Mast cells resting" = "Mast resting",
          "Mast cells activated" = "Mast activated",
          "Eosinophils" = "Eosinophils",
          "Neutrophils" = "Neutrophils")
stopifnot(all(sub %in% names(disp)))

ggplot(d, aes(sample, frac, fill = cell)) +
  geom_col(width = 1, colour = NA) +
  facet_wrap(~ cohort, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = fam_col, labels = disp, name = NULL) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.25)) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  labs(x = "Samples ordered by group", y = "Estimated cell fraction") +
  theme_ajr() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text = element_text(size = PT_BODY, hjust = 0),
        legend.position = "right",
        legend.text = element_text(size = PT_BODY),
        legend.key.size = unit(2.2, "mm"),
        legend.spacing.y = unit(0.6, "mm"),
        legend.margin = margin(0, 0, 0, 1),
        panel.spacing.y = unit(2, "mm"),
        plot.margin = margin(2, 1, 2, 2))

fig <- ggplot2::last_plot()
save_fig(fig, "FigureS24", dir = out_dir("supplementary"), height_mm = 208)
msg("FigureS24 written")
