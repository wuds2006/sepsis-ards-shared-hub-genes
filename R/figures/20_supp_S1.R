# ---------------------------------------------------------------------------
# 20_supp_S1.R -- Figure S1. WGCNA soft-threshold selection and module structure
#   A  scale-free topology fit, ARDS discovery cohort
#   B  scale-free topology fit, sepsis discovery cohort
#   C  module sizes in the ARDS discovery network
#   D  module preservation of the ARDS modules in the sepsis cohort
# Input  : data/frozen/01_DEG/*_gene_expr.rds, data/frozen/02_WGCNA
# Output : Figures/supplementary/FigureS1.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")
suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)

EX <- function(...) file.path("data", "frozen", "01_DEG", ...)
WG <- function(...) file.path("data", "frozen", "02_WGCNA", ...)

# pickSoftThreshold is re-run on the frozen expression matrices so that the
# curve, the chosen power and the connectivity are traceable to the inputs.
# The fit is evaluated on the 5,000 most variable genes: the full 22,880-gene
# adjacency matrix does not fit in memory on this machine, and the fit is a
# property of the bulk of the transcriptome rather than of single genes.
N_TOP <- 5000

soft_fit <- function(file, chosen) {
  x <- readRDS(EX(file))
  madv <- apply(x, 1, mad, na.rm = TRUE)
  x <- x[order(-madv)[seq_len(min(N_TOP, nrow(x)))], , drop = FALSE]
  x <- x[apply(x, 1, var, na.rm = TRUE) > 0, , drop = FALSE]
  sft <- pickSoftThreshold(t(x), powerVector = c(1:10, seq(12, 20, 2)),
                           networkType = "unsigned", verbose = 0)
  f <- sft$fitIndices
  data.frame(power = f$Power, rsq = -sign(f$slope) * f$SFT.R.sq,
             conn = f$mean.k., chosen = f$Power == chosen)
}

pw_a <- soft_fit("GSE76293_gene_expr.rds", 10)
pw_b <- soft_fit("GSE69528_gene_expr.rds", 10)

fit_panel <- function(d, sub) {
  ggplot(d, aes(power, rsq)) +
    geom_hline(yintercept = 0.85, linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    geom_vline(xintercept = d$power[d$chosen], linetype = "dotted",
               linewidth = 0.3, colour = pal[["brick"]]) +
    geom_line(linewidth = 0.4, colour = pal[["ink"]]) +
    geom_point(size = 0.9, colour = pal[["ink"]]) +
    scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Soft-thresholding power", y = "Scale-free topology fit (R2)",
         subtitle = sub) +
    theme_ajr() +
    theme(plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

# --- panel C ---------------------------------------------------------------

panel_c <- function() {
  gm <- readRDS(WG("GSE76293_wgcna.rds"))$gene_module
  tb <- as.data.frame(table(gm$module))
  names(tb) <- c("module", "n")
  tb <- tb[tb$module != "grey", ]
  tb <- tb[order(-tb$n), ][1:12, ]
  tb$module <- factor(tb$module, levels = rev(tb$module))

  ggplot(tb, aes(n, module)) +
    geom_col(width = 0.68, fill = pal[["steel"]], colour = NA) +
    geom_text(aes(label = n), hjust = -0.18, family = FONT_BODY,
              size = 7 / .pt, colour = pal[["ink"]]) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = "Genes in module", y = NULL, subtitle = "ARDS discovery network") +
    theme_ajr() +
    theme(plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

# --- panel D ---------------------------------------------------------------

panel_d <- function() {
  d <- read.csv(WG("preservation", "modulePreservation_GSE76293_in_GSE69528.csv"))
  d <- d[!is.na(d$Zsummary), ]
  ggplot(d, aes(moduleSize, Zsummary)) +
    geom_hline(yintercept = 10, linetype = "dashed", linewidth = 0.3,
               colour = pal[["sage"]]) +
    geom_hline(yintercept = 2, linetype = "dotted", linewidth = 0.3,
               colour = pal[["slate"]]) +
    geom_point(size = 1.2, colour = pal[["plum"]]) +
    ggrepel::geom_text_repel(aes(label = module), size = 7 / .pt,
                             family = FONT_BODY, colour = pal[["ink"]],
                             min.segment.length = 0, segment.size = 0.2,
                             box.padding = 0.18, max.overlaps = Inf, seed = 5) +
    scale_x_log10() +
    labs(x = "Module size (genes)", y = "Preservation Zsummary",
         subtitle = "ARDS modules in sepsis") +
    theme_ajr() +
    theme(plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

# --- assemble --------------------------------------------------------------

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid,
                c(list(fit_panel(pw_a, "ARDS discovery (GSE76293)"),
                       fit_panel(pw_b, "Sepsis discovery (GSE69528)"),
                       ncol = 2, labels = c("A", "B")), lab_args))
row2 <- do.call(cowplot::plot_grid,
                c(list(panel_c(), panel_d(), ncol = 2, labels = c("C", "D")), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1))
save_fig(fig, "FigureS1", dir = out_dir("supplementary"), height_mm = 200)
msg("FigureS1 written")
