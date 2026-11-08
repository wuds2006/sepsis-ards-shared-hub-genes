# ---------------------------------------------------------------------------
# 38_supp_S18.R -- Figure S18. Group-stratified virtual knockout
#   A  number of significantly re-regulated genes per knockout and group
#   B  regulatory shift of the most affected genes across knockout settings
# Input  : data/frozen/08_singlecell/knockout_stratified
# Output : Figures/supplementary/FigureS18.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

KS <- function(...) file.path("data", "frozen", "08_singlecell", "knockout_stratified", ...)

files <- list.files(KS(), pattern = "_diffRegulation\\.csv$")
meta <- do.call(rbind, lapply(files, function(f) {
  p <- strsplit(sub("^KO_|\\.csv$", "", f), "_")[[1]]
  data.frame(file = f, cohort = p[1], group = p[2], ko = p[3])
}))
meta$lab <- sprintf("%s %s %s", meta$cohort, meta$group, meta$ko)

d <- do.call(rbind, lapply(seq_len(nrow(meta)), function(i) {
  x <- read.csv(KS(meta$file[i]), stringsAsFactors = FALSE)
  x$lab <- meta$lab[i]
  x$grp <- meta$group[i]
  x$ko <- meta$ko[i]
  x
}))
d$sig <- !is.na(d$p.adj) & d$p.adj < 0.05

# --- panel A ---------------------------------------------------------------

panel_a <- function() {
  a <- data.frame(lab = names(table(d$lab)), n = as.numeric(tapply(d$sig, d$lab, sum)),
                  grp = tapply(d$grp, d$lab, unique))
  a$grp <- factor(a$grp, levels = c("ARDS", "Sepsis", "HC"))
  a$lab <- factor(a$lab, levels = a$lab[order(a$n, a$lab)])

  ggplot(a, aes(n, lab, fill = grp)) +
    geom_col(width = 0.64, colour = NA) +
    geom_text(aes(label = n), hjust = -0.4, family = FONT_BODY,
              size = 7 / .pt, colour = pal[["ink"]]) +
    scale_fill_manual(values = c(ARDS = pal[["steel"]], Sepsis = pal[["brick"]],
                                 HC = pal[["sage"]]), name = NULL) +
    scale_x_continuous(breaks = 0:3, expand = expansion(mult = c(0, 0.14))) +
    labs(x = "Genes with adjusted P < 0.05", y = NULL) +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.6, "mm"),
          axis.text.y = element_text(size = PT_BODY))
}

# --- panel B ---------------------------------------------------------------

panel_b <- function() {
  # twelve genes with the largest regulatory shift in any knockout setting
  z <- aggregate(abs(Z) ~ gene, data = d, FUN = max)
  keep <- as.character(z$gene[order(-z$`abs(Z)`)][1:12])
  m <- d[d$gene %in% keep, c("gene", "lab", "Z")]
  m$Z <- pmax(pmin(m$Z, 8), -8)
  m$gene <- factor(m$gene, levels = rev(keep))

  ggplot(m, aes(lab, gene, fill = Z)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_gradientn(colours = ramp_div, limits = c(-8, 8),
                         oob = scales::squish, na.value = "#F2F2F2", name = "Z") +
    labs(x = NULL, y = NULL) +
    theme_ajr() +
    theme(axis.text.x = element_text(angle = 32, hjust = 1, vjust = 1,
                                     lineheight = 0.85),
          axis.text.y = element_text(face = "italic"),
          axis.line = element_blank(), axis.ticks = element_blank(),
          plot.margin = margin(2, 3, 2, 6),
          legend.key.height = unit(6, "mm"), legend.key.width = unit(2.2, "mm"),
          legend.title = element_text(size = PT_BODY))
}

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid, c(list(panel_a(), labels = "A"), lab_args))
row2 <- do.call(cowplot::plot_grid, c(list(panel_b(), labels = "B"), lab_args))

fig <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(0.95, 1.05))
save_fig(fig, "FigureS18", dir = out_dir("supplementary"), height_mm = 206)
msg("FigureS18 written")
