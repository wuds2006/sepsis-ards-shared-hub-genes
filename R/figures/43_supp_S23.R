# ---------------------------------------------------------------------------
# 43_supp_S23.R -- Figure S23. Functional enrichment of the 179 shared
# candidate genes
#   A  GO biological process
#   B  KEGG pathways
#   C  Reactome pathways
# Input  : data/frozen/03_intersection, data/frozen/11_network_pharma
# Output : Figures/supplementary/FigureS23.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

IN <- function(...) file.path("data", "frozen", "03_intersection", ...)
NP <- function(...) file.path("data", "frozen", "11_network_pharma", ...)

enrich_panel <- function(d, term_col, p_col, n, sub, wrap = 34,
                         xlab = "-log10 adjusted P") {
  d <- d[order(d[[p_col]]), ][seq_len(min(n, nrow(d))), ]
  d$term <- d[[term_col]]
  d$term <- paste0(toupper(substr(d$term, 1, 1)), substr(d$term, 2, nchar(d$term)))
  d$nlp <- -log10(d[[p_col]])
  d$term <- factor(d$term, levels = rev(d$term))

  ggplot(d, aes(nlp, term)) +
    geom_col(width = 0.64, fill = pal[["steel"]], colour = NA) +
    geom_text(aes(label = Count), hjust = -0.25, family = FONT_BODY,
              size = 7 / .pt, colour = pal[["ink"]]) +
    scale_y_discrete(labels = function(v) stringr::str_wrap(v, wrap)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(x = xlab, y = NULL, subtitle = sub) +
    theme_ajr() +
    theme(axis.text.y = element_text(lineheight = 0.85),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

go <- read.csv(IN("Shared_CGs_GO_BP.csv"), stringsAsFactors = FALSE)
kg <- read.csv(IN("Shared_CGs_KEGG.csv"), stringsAsFactors = FALSE)
re <- read.csv(NP("Reactome_2022_enrichment.csv"), stringsAsFactors = FALSE)
re$Count <- NA_integer_
re$Description <- sub(" R-HSA-\\d+$", "", re$term)

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid,
                c(list(enrich_panel(go, "Description", "p.adjust", 12,
                                    "GO biological process"),
                       labels = "A"), lab_args))
row2 <- do.call(cowplot::plot_grid,
                c(list(enrich_panel(kg, "Description", "p.adjust", 10,
                                    "KEGG pathway"),
                       enrich_panel(re, "Description", "pvalue", 10,
                                    "Reactome pathway", xlab = "-log10 P"),
                       ncol = 2, labels = c("B", "C")), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1))
save_fig(fig, "FigureS23", dir = out_dir("supplementary"), height_mm = 204)
msg("FigureS23 written")
