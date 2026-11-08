# ---------------------------------------------------------------------------
# 28_supp_S9.R -- Figure S9. Hub-gene expression against CIBERSORT cell fractions
#   Spearman correlation matrices for the four bulk cohorts
# Input  : data/frozen/07_hub_immune/cibersort/*_hub_ciber_cor.csv
# Output : Figures/supplementary/FigureS9.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

CI <- function(...) file.path("data", "frozen", "07_hub_immune", "cibersort", ...)

coh <- c(GSE76293 = "ARDS discovery (GSE76293)",
         GSE69528 = "Sepsis discovery (GSE69528)",
         GSE243066 = "ARDS validation (GSE243066)",
         GSE100159 = "Sepsis validation (GSE100159)")

long <- do.call(rbind, lapply(names(coh), function(cd) {
  m <- read.csv(CI(paste0(cd, "_hub_ciber_cor.csv")), row.names = 1,
                check.names = FALSE)
  m <- as.matrix(m[HUB6, , drop = FALSE])
  d <- as.data.frame(as.table(m))
  names(d) <- c("gene", "cell", "rho")
  d$cohort <- coh[[cd]]
  d
}))
long$cohort <- factor(long$cohort, levels = coh)
long$gene <- factor(long$gene, levels = rev(HUB6))
# keep cell types with at least one non-missing estimate in every cohort
# cell types that are absent from a cohort (all-zero fractions) carry no
# estimate there; drop only those that are missing in every cohort
keep <- names(which(tapply(!is.na(long$rho), long$cell, sum) > 0))
long <- long[long$cell %in% keep, ]
ord <- names(sort(tapply(long$rho, long$cell, function(v) mean(abs(v), na.rm = TRUE))))
long$cell <- factor(long$cell, levels = ord)

ggplot(long, aes(cell, gene, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  facet_wrap(~ cohort, ncol = 1) +
  scale_fill_gradientn(colours = ramp_div, limits = c(-0.6, 0.6),
                       oob = scales::squish, na.value = "#E8E8E8",
                       name = "Spearman rho") +
  scale_x_discrete(labels = function(v) stringr::str_wrap(v, 14)) +
  labs(x = NULL, y = NULL) +
  theme_ajr() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1,
                                   lineheight = 0.85),
          axis.text.y = element_text(face = "italic"),
          strip.text = element_text(size = PT_BODY, hjust = 0),
          axis.line = element_blank(), axis.ticks = element_blank(),
          panel.spacing.y = unit(2.5, "mm"),
          legend.key.height = unit(8, "mm"), legend.key.width = unit(2.2, "mm"),
          legend.title = element_text(size = PT_BODY))

fig <- ggplot2::last_plot()
save_fig(fig, "FigureS9", dir = out_dir("supplementary"), height_mm = 210)
msg("FigureS9 written")
