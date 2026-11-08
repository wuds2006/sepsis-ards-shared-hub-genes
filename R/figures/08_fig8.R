# ---------------------------------------------------------------------------
# 08_fig8.R -- Figure 8, shared immune microenvironment
#   A  six-gene signature score in the two discovery cohorts
#   B  CIBERSORT cell-fraction shift, disease versus control, three cohorts
#   C, D  top immune pathways from GSVA in the two discovery cohorts
# Input  : data/frozen/07_hub_immune
# Output : Figures/main/Figure8.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

FZ <- function(...) file.path("data", "frozen", "07_hub_immune", ...)

# --- panel A ---------------------------------------------------------------

panel_a <- function() {
  st <- read.csv(FZ("hub6_GSVA_stats.csv"), stringsAsFactors = FALSE)
  coh <- c(GSE76293 = "ARDS (GSE76293)", GSE69528 = "Sepsis (GSE69528)")
  d <- do.call(rbind, lapply(names(coh), function(cd) {
    x <- read.csv(FZ(paste0(cd, "_hub6_GSVA.csv")), stringsAsFactors = FALSE)
    x$cohort <- coh[[cd]]
    x
  }))
  d$grp <- factor(ifelse(d$group == "Control", "Control", "Disease"),
                  levels = c("Control", "Disease"))
  d$cohort <- factor(d$cohort, levels = coh)

  n <- table(d$cohort, d$grp)
  xlab <- setNames(sprintf("%s\nn = %d / %d", levels(d$cohort),
                           n[levels(d$cohort), "Disease"],
                           n[levels(d$cohort), "Control"]), levels(d$cohort))

  ann <- data.frame(
    cohort = factor(coh[st$dataset], levels = coh),
    lab = sprintf("P = %s", formatC(st$P, format = "e", digits = 1))
  )
  ann$y <- as.numeric(tapply(d$score, d$cohort, max)[as.character(ann$cohort)]) + 0.30

  ggplot(d, aes(cohort, score)) +
    geom_violin(aes(fill = grp), width = 0.72, linewidth = 0.3,
                colour = "#4A4A4A", scale = "width",
                position = position_dodge(0.82)) +
    geom_boxplot(aes(fill = grp), width = 0.10, linewidth = 0.28,
                 outlier.shape = NA, fill = "white",
                 position = position_dodge(0.82)) +
    geom_point(aes(fill = grp), shape = 21, size = 0.45, stroke = 0.1,
               colour = "#6E6E6E", alpha = 0.55,
               position = position_jitterdodge(jitter.width = 0.13,
                                               dodge.width = 0.82, seed = 3)) +
    geom_text(data = ann, aes(x = cohort, y = y, label = lab), inherit.aes = FALSE,
              family = FONT_BODY, size = 7 / .pt, colour = pal[["ink"]]) +
    scale_fill_manual(values = c(Control = pal[["sage"]], Disease = pal[["brick"]]),
                      name = NULL) +
    scale_x_discrete(labels = xlab) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.16))) +
    labs(x = NULL, y = "Six-gene GSVA score") +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -3, 0),
          axis.text.x = element_text(lineheight = 0.95))
}

# --- panel B ---------------------------------------------------------------

panel_b <- function() {
  coh <- c(GSE76293 = "ARDS (GSE76293)", GSE69528 = "Sepsis (GSE69528)",
           GSE100159 = "Sepsis (GSE100159)")
  d <- do.call(rbind, lapply(names(coh), function(cd) {
    x <- read.csv(FZ("cibersort", paste0(cd, "_celltype_wilcox.csv")),
                  stringsAsFactors = FALSE)
    x$cohort <- coh[[cd]]
    x
  }))
  d <- d[!is.na(d$log2FC) & !is.na(d$padj), ]
  d$sig <- factor(ifelse(d$padj < 0.05, "adjusted P < 0.05", "not significant"),
                  levels = c("adjusted P < 0.05", "not significant"))
  d$nlp <- pmax(-log10(d$padj), 0)
  ord <- names(sort(tapply(abs(d$log2FC), d$cell, max)))
  d$cell <- factor(d$cell, levels = ord)
  d$cohort <- factor(d$cohort, levels = coh)

  ggplot(d, aes(log2FC, cell)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    geom_point(aes(size = nlp, colour = cohort, shape = sig), stroke = 0.35,
               position = position_dodge(width = 0.72)) +
    scale_colour_manual(values = c(pal[["brick"]], pal[["steel"]], pal[["olive"]]),
                        labels = c("GSE76293", "GSE69528", "GSE100159"),
                        name = NULL) +
    scale_shape_manual(values = c(16, 1), name = NULL,
                       labels = c("adj. P < 0.05", "not significant")) +
    scale_size_area(max_size = 2.1, breaks = c(1, 6),
                    name = "-log10 adj. P") +
    labs(x = "log2 fold change in cell fraction (disease vs control)", y = NULL) +
    theme_ajr() +
    theme(legend.position = "top", legend.box = "horizontal",
          legend.margin = margin(0, 0, -5, 0),
          legend.spacing.x = unit(2.5, "mm"),
          legend.key.size = unit(2.6, "mm"))
}

# --- panels C, D -----------------------------------------------------------

panel_cd <- function(cd, sub) {
  d <- read.csv(FZ("gsva_top20", paste0(cd, "_GSVA_top20.csv")),
                stringsAsFactors = FALSE)
  d <- d[order(d$logFC), ]
  d$pathway <- gsub("_", " ", sub("^GOBP_|^HALLMARK_|^REACTOME_", "", d$pathway))
  d$pathway <- tolower(d$pathway)
  d$pathway <- paste0(toupper(substr(d$pathway, 1, 1)), substr(d$pathway, 2, nchar(d$pathway)))
  d$sig <- ifelse(d$adj.P.Val < 0.05, "*", "")

  ggplot(d, aes(logFC, reorder(pathway, logFC))) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = pal[["ink"]]) +
    geom_col(aes(fill = logFC), width = 0.68, colour = NA) +
    geom_text(aes(x = ifelse(logFC > 0, logFC + 0.03, logFC - 0.03), label = sig),
              hjust = ifelse(d$logFC > 0, 0, 1), family = FONT_BODY,
              size = 7 / .pt, colour = pal[["ink"]]) +
    scale_fill_gradientn(colours = ramp_div, guide = "none") +
    scale_x_continuous(breaks = scales::pretty_breaks(4),
                       expand = expansion(mult = c(0.18, 0.18))) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, 34)) +
    labs(x = "GSVA log2 fold change (disease vs control)", y = NULL,
         subtitle = sub) +
    theme_ajr() +
    theme(axis.text.y = element_text(lineheight = 0.88, size = PT_BODY),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold",
                                       margin = margin(0, 0, 1.5, 0)))
}

# --- assemble --------------------------------------------------------------

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid, c(list(panel_a(), labels = "A"), lab_args))
row2 <- do.call(cowplot::plot_grid, c(list(panel_b(), labels = "B"), lab_args))
row3 <- do.call(cowplot::plot_grid,
                c(list(panel_cd("GSE76293", "ARDS (GSE76293)"),
                       panel_cd("GSE69528", "Sepsis (GSE69528)"),
                       ncol = 2, labels = c("C", "D")), lab_args))

fig8 <- cowplot::plot_grid(row1, row2, row3, ncol = 1,
                           rel_heights = c(0.60, 0.80, 1.62))

save_fig(fig8, "Figure8", width_mm = FIG_W_MM, height_mm = 210)
msg("Figure8 written")
