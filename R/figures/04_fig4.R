# ---------------------------------------------------------------------------
# 04_fig4.R -- Figure 4, expression and single-gene enrichment of the hub genes
#   A-D  per-sample expression of the six hub genes in the four cohorts
#   E    single-gene GSEA: normalised enrichment scores of the leading
#        Hallmark pathways, six genes by four cohorts
# Inputs : data/derived/hub6_expression_long.csv (built by 01_prep_tables.R),
#          data/frozen/12_single_gene_gsea/single_gene_gsea_all_results_4cohorts.csv
# Output : Figures/main/Figure4.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")
suppressPackageStartupMessages(library(forcats))

FZ <- function(...) file.path("data", "frozen", ...)

expr <- read.csv(file.path("data", "derived", "hub6_expression_long.csv"),
                 stringsAsFactors = FALSE)
expr$gene <- factor(expr$gene, levels = HUB6)

cohort_lab <- c(GSE76293 = "ARDS (GSE76293)", GSE69528 = "Sepsis (GSE69528)",
                GSE100159 = "Sepsis (GSE100159)", GSE243066 = "ARDS (GSE243066)")

# --- panels A-D ------------------------------------------------------------

box_panel <- function(cd) {
  d <- expr[expr$cohort == cd, ]
  d$group <- factor(d$group, levels = c("Control", "Disease"))
  n <- tapply(d$sample, d$group, function(x) length(unique(x)))

  # two-sided Wilcoxon per gene, BH across the six genes of that cohort
  pv <- sapply(HUB6, function(g) {
    x <- d$expr[d$gene == g & d$group == "Control"]
    y <- d$expr[d$gene == g & d$group == "Disease"]
    if (length(x) < 3 || length(y) < 3) return(NA_real_)
    wilcox.test(y, x)$p.value
  })
  pv <- p.adjust(pv, "BH")
  st <- data.frame(gene = factor(HUB6, levels = HUB6),
                   p = pv,
                   y = tapply(d$expr, d$gene, max)[HUB6] +
                       diff(range(d$expr)) * 0.04)
  st$lab <- ifelse(is.na(st$p), "", as.character(star_p(st$p)))

  p <- ggplot(d, aes(gene, expr, fill = group)) +
    geom_boxplot(width = 0.66, position = position_dodge(width = 0.72),
                 outlier.shape = NA, linewidth = 0.3, colour = "#4A4A4A") +
    geom_point(position = position_jitterdodge(jitter.width = 0.12,
                                               dodge.width = 0.72, seed = 5),
               size = 0.4, colour = "#4A4A4A", stroke = 0) +
    scale_fill_manual(values = c(Control = pal[["sage"]], Disease = pal[["brick"]]),
                      labels = c(sprintf("Control (n = %d)", n[["Control"]]),
                                 sprintf("Disease (n = %d)", n[["Disease"]])),
                      name = NULL) +
    geom_text(data = st, aes(gene, y, label = lab), inherit.aes = FALSE,
              family = FONT_BODY, size = 7 / .pt, colour = pal[["ink"]],
              vjust = 0) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.10))) +
    labs(x = NULL, y = "log2 expression", subtitle = cohort_lab[[cd]]) +
    theme_ajr() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, face = "italic"),
          legend.position = "top", legend.direction = "horizontal",
          legend.margin = margin(0, 0, -3, 0),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
  p
}

# --- panel E ---------------------------------------------------------------

gsea_e <- function() {
  g <- read.csv(FZ("12_single_gene_gsea", "single_gene_gsea_all_results_4cohorts.csv"),
                stringsAsFactors = FALSE)

  # leading pathways: the two strongest per gene in each cohort, pooled, then
  # kept only if they reach the top two in at least two cohorts
  # pathways must be scored in all four cohorts, otherwise the grid would
  # carry empty cells that are not biological missingness
  complete <- names(which(tapply(g$NES, g$pathway,
                                 function(x) sum(!is.na(x))) == length(HUB6) * 4))
  top <- g[g$pathway %in% complete, ] %>%
    group_by(cohort, hub_gene) %>%
    slice_max(order_by = NES, n = 2) %>%
    ungroup()
  rank_tab <- sort(table(top$pathway), decreasing = TRUE)
  keep <- names(rank_tab)[seq_len(min(10, length(rank_tab)))]

  d <- g[g$pathway %in% keep, ]
  d$pathway_lab <- gsub("^HALLMARK_", "", d$pathway)
  d$pathway_lab <- gsub("_", " ", d$pathway_lab)
  d$pathway_lab <- tolower(d$pathway_lab)
  d$pathway_lab <- paste0(toupper(substr(d$pathway_lab, 1, 1)),
                          substr(d$pathway_lab, 2, nchar(d$pathway_lab)))
  d$pathway_lab <- factor(d$pathway_lab, levels = rev(sort(unique(d$pathway_lab))))
  d$hub_gene <- factor(d$hub_gene, levels = HUB6)
  d$cohort <- factor(cohort_lab[d$cohort],
                     levels = cohort_lab[c("GSE76293", "GSE69528",
                                           "GSE100159", "GSE243066")])
  d$sig <- ifelse(d$padj < 0.05, "*", "")

  lim <- max(abs(d$NES), na.rm = TRUE)

  ggplot(d, aes(hub_gene, pathway_lab, fill = NES)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sig), family = FONT_BODY, size = 7 / .pt,
              colour = pal[["ink"]], vjust = 0.78) +
    facet_wrap(~ cohort, nrow = 1) +
    scale_fill_gradientn(colours = ramp_div, limits = c(-lim, lim),
                         name = "NES") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_ajr() +
    theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5,
                                        face = "italic"),
          axis.line = element_blank(), axis.ticks = element_blank(),
          strip.text = element_text(size = PT_BODY, hjust = 0),
          panel.spacing.x = unit(4, "mm"),
          legend.key.height = unit(12, "pt"),
          legend.key.width = unit(4, "pt"))
}

# --- assemble --------------------------------------------------------------

pA <- box_panel("GSE76293")
pB <- box_panel("GSE69528")
pC <- box_panel("GSE100159")
pD <- box_panel("GSE243066")
pE <- gsea_e()

# cowplot places the 9 pt bold tags at one fixed corner of every subplot, which
# keeps them aligned even though the panels have very different axis furniture
lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid, c(list(pA, pB, ncol = 2, labels = c("A", "B"),
                                           align = "h", axis = "tb"), lab_args))
row2 <- do.call(cowplot::plot_grid, c(list(pC, pD, ncol = 2, labels = c("C", "D"),
                                           align = "h", axis = "tb"), lab_args))
row3 <- do.call(cowplot::plot_grid, c(list(pE, ncol = 1, labels = "E"), lab_args))

fig4 <- cowplot::plot_grid(row1, row2, row3, ncol = 1,
                           rel_heights = c(1, 1, 1.05))

save_fig(fig4, "Figure4", width_mm = FIG_W_MM, height_mm = 212)
msg("Figure4 written")
