# ---------------------------------------------------------------------------
# 01_fig1.R -- Figure 1
#   A  study design and analysis flow
#   B  volcano, ARDS discovery cohort (GSE76293, 12 vs 12)
#   C  volcano, sepsis discovery cohort (GSE69528, 83 vs 28)
#   D  WGCNA module-trait correlation, both discovery cohorts
#   E  four-way intersection of DEG and disease-module gene sets
# Inputs : data/frozen/01_DEG, 02_WGCNA, 03_intersection
# Output : Figures/main/Figure1.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")
suppressPackageStartupMessages({
  library(ggrepel); library(VennDiagram); library(forcats)
})

FZ <- function(...) file.path("data", "frozen", ...)

# --- panel A: flow diagram -------------------------------------------------

flow_a <- function() {
  box <- function(x, y, w, h, lab, fill = "#EFF2F4", col = pal[["slate"]],
                  lwd = 0.4, fontface = "plain") {
    list(x = x, y = y, w = w, h = h, lab = lab, fill = fill, col = col,
         lwd = lwd, fontface = fontface)
  }
  bx <- list(
    box( 24, 28, 44, 12, "Discovery cohorts\nARDS  GSE76293 (12 vs 12)\nsepsis  GSE69528 (83 vs 28)"),
    box( 70, 28, 36, 12, "Differential expression\n557 ARDS DEGs\n1,458 sepsis DEGs"),
    box(112, 28, 36, 12, "WGCNA\n13 ARDS modules\n11 sepsis modules"),
    box(157, 28, 42, 12, "DEG and module intersection\n179 shared\ncandidate genes",
        fill = "#F4E7E7", col = pal[["brick"]], lwd = 0.6, fontface = "bold"),
    box( 34, 8, 60, 11, "Machine learning\n12 classifiers per cohort (11 estimable)\nSHAP ranking"),
    box( 99, 8, 40, 11, "Six hub genes\nCDK5RAP2  CD177  GRB10\nUPP1  EXOSC4  VNN1",
        fill = "#F4E7E7", col = pal[["brick"]], lwd = 0.6, fontface = "bold"),
    box(154, 8, 44, 11, "Validation and interpretation\nGSE243066 · GSE100159\nGSE65682 · single-cell atlases")
  )
  boxes <- data.frame(
    xmin = sapply(bx, `[[`, "x") - sapply(bx, `[[`, "w") / 2,
    xmax = sapply(bx, `[[`, "x") + sapply(bx, `[[`, "w") / 2,
    ymin = sapply(bx, `[[`, "y") - sapply(bx, `[[`, "h") / 2,
    ymax = sapply(bx, `[[`, "y") + sapply(bx, `[[`, "h") / 2,
    fill = sapply(bx, `[[`, "fill"),
    col  = sapply(bx, `[[`, "col"),
    lwd  = sapply(bx, `[[`, "lwd"),
    lab  = sapply(bx, `[[`, "lab"),
    fontface = sapply(bx, `[[`, "fontface"),
    stringsAsFactors = FALSE
  )

  # straight arrows between neighbours in a row, elbow from row 1 to row 2
  arr <- data.frame(
    x    = c(46.6, 88.6, 130.6, 64.6, 119.6, 34),
    y    = c(28,   28,   28,    8,    8,     18.4),
    xend = c(51.4, 93.4, 135.4, 78.6, 124.4, 34),
    yend = c(28,   28,   28,    8,    8,     13.8)
  )
  elbow <- data.frame(x = c(157, 157, 34), y = c(22, 18.4, 18.4))

  ggplot() +
    geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin,
                                ymax = ymax, fill = fill, colour = col,
                                linewidth = lwd)) +
    scale_fill_identity() + scale_colour_identity() + scale_linewidth_identity() +
    geom_text(data = boxes, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
                                label = lab, fontface = fontface),
              family = FONT_BODY, size = 7 / .pt, lineheight = 1.18) +
    geom_path(data = elbow, aes(x, y), colour = pal[["slate"]], linewidth = 0.35) +
    geom_segment(data = arr, aes(x = x, y = y, xend = xend, yend = yend),
                 colour = pal[["slate"]], linewidth = 0.35,
                 arrow = arrow(type = "closed", length = unit(1.2, "mm"))) +
    coord_cartesian(xlim = c(0, 182), ylim = c(0, 36), expand = FALSE, clip = "off") +
    theme_void_ajr() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

# --- panels B, C: volcano --------------------------------------------------

volcano_bc <- function(deg_file, subtitle) {
  d <- read.csv(deg_file, stringsAsFactors = FALSE)
  d$grp <- ifelse(d$adj.P.Val < 0.05 & d$logFC > 1, "Up",
           ifelse(d$adj.P.Val < 0.05 & d$logFC < -1, "Down", "NS"))
  d$grp <- factor(d$grp, levels = c("Up", "Down", "NS"))
  d$y <- -log10(d$adj.P.Val)
  n_up <- sum(d$grp == "Up"); n_dn <- sum(d$grp == "Down")
  lab <- d[d$gene %in% HUB6, ]

  ggplot(d, aes(logFC, y)) +
    geom_point(aes(colour = grp), size = 0.55, stroke = 0) +
    scale_colour_manual(values = c(Up = pal[["brick"]], Down = pal[["steel"]],
                                   NS = "#C9C9C9"),
                        labels = c(paste0("Up, n = ", n_up),
                                   paste0("Down, n = ", n_dn), "Not significant"),
                        name = NULL) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    geom_point(data = lab, shape = 21, size = 1.1, stroke = 0.4,
               colour = pal[["ink"]], fill = NA) +
    geom_text_repel(data = lab, aes(label = gene), family = FONT_BODY,
                    fontface = "italic", size = 7 / .pt, colour = pal[["ink"]],
                    min.segment.length = 0, segment.size = 0.2,
                    box.padding = 0.25, point.padding = 0.15, max.overlaps = 20,
                    seed = 1) +
    scale_x_continuous(limits = c(-6.5, 6.5), breaks = seq(-6, 6, 2)) +
    labs(subtitle = subtitle, x = "log2 fold change",
         y = "-log10 adjusted P") +
    theme_ajr() +
    theme(legend.position = "top", legend.direction = "horizontal",
          legend.margin = margin(0, 0, -2, 0),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"),
          plot.margin = margin(2, 4, 2, 3))
}

# --- panel D: module-trait correlation -------------------------------------

module_trait_d <- function() {
  one <- function(rds, cohort) {
    m <- readRDS(rds)$mt
    m$cohort <- cohort
    m
  }
  mt <- rbind(
    one(FZ("02_WGCNA", "GSE76293_wgcna.rds"), "ARDS (GSE76293)"),
    one(FZ("02_WGCNA", "GSE69528_wgcna.rds"), "Sepsis (GSE69528)")
  )
  mt <- mt %>%
    mutate(cohort = factor(cohort, levels = c("ARDS (GSE76293)", "Sepsis (GSE69528)")),
           module_lab = fct_reorder(paste0("ME", module), cor),
           p_lab = star_p(p),
           lab = sprintf("%.2f%s", cor, ifelse(p_lab == "ns", "", p_lab)),
           hjust = ifelse(cor > 0, -0.12, 1.12),
           dir = ifelse(cor > 0, "Positive", "Negative"))

  ggplot(mt, aes(x = cor, y = module_lab)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = pal[["ink"]]) +
    geom_col(aes(fill = dir), width = 0.62, colour = NA) +
    geom_text(aes(label = lab, hjust = hjust), family = FONT_BODY,
              size = 7 / .pt, colour = pal[["ink"]]) +
    facet_wrap(~ cohort, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = c(Positive = pal[["brick"]],
                                 Negative = pal[["steel"]]), guide = "none") +
    scale_x_continuous(limits = c(-1.25, 1.25), breaks = seq(-1, 1, 0.5)) +
    labs(x = "Module-trait correlation with disease status", y = NULL) +
    theme_ajr() +
    theme(panel.spacing.x = unit(10, "mm"),
          strip.text = element_text(hjust = 0, size = PT_BODY),
          axis.text.y = element_text(size = PT_BODY, colour = pal[["slate"]]),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = margin(2, 5, 2, 3))
}

# --- panel E: four-set intersection ----------------------------------------

venn_e <- function() {
  rd <- function(f) unique(trimws(readLines(FZ(f), warn = FALSE)))
  clean <- function(x) unique(na.omit(trimws(x))[trimws(x) != ""])
  sets <- list(
    "ARDS DEG"      = clean(read.csv(FZ("01_DEG", "GSE76293_DEG.csv"))$gene),
    "ARDS module"   = clean(rd("02_WGCNA/GSE76293_disease_module_genes.txt")),
    "Sepsis DEG"    = clean(read.csv(FZ("01_DEG", "GSE69528_DEG.csv"))$gene),
    "Sepsis module" = clean(rd("02_WGCNA/GSE69528_disease_module_genes.txt"))
  )
  stopifnot(length(Reduce(intersect, sets)) ==
              nrow(read.csv(FZ("03_intersection", "Shared_CGs.csv"))))

  cols <- c(pal[["brick"]], pal[["ochre"]], pal[["steel"]], pal[["sage"]])
  v <- venn.diagram(sets, filename = NULL, disable.logging = TRUE,
                    fill = cols, alpha = 0.38, col = cols, lwd = 0.6,
                    cex = 7/12, cat.cex = 7/12, cat.fontface = "plain",
                    cat.pos = c(-25, 28, -38, 38),
                    cat.dist = c(0.22, 0.24, 0.15, 0.15),
                    cat.col = rep(pal[["ink"]], 4),
                    fontfamily = FONT_BODY, cat.fontfamily = FONT_BODY,
                    margin = 0.16, euler.d = FALSE, scaled = FALSE)
  wrap_elements(full = grobTree(v))
}

# --- assemble --------------------------------------------------------------

# Panel letters go through plot_annotation(), not tag_plot(). tag_plot() pins the
# letter at a fixed npc position inside each panel's own gtable, and that gtable
# includes the y-axis label and tick text on the left. The offset therefore
# depends on how wide a panel's left decoration happens to be: A, D and E landed
# at the outer figure edge while B and C landed a further ~10 % in, and C ended
# up crowding the right margin of panel B. plot_annotation() reserves one tag
# column across the whole layout, so all five letters sit on a single vertical
# rule at the same offset above their panel.
pA <- flow_a()
pB <- volcano_bc(FZ("01_DEG", "GSE76293_DEG_all.csv"),
                 "ARDS vs control (GSE76293)")
pC <- volcano_bc(FZ("01_DEG", "GSE69528_DEG_all.csv"),
                 "Sepsis vs control (GSE69528)")
pD <- module_trait_d()
pE <- venn_e()

# The outer margin goes through plot_annotation(theme =), which wraps the whole
# layout; applying the same theme with `&` only reaches the sub-plots, and the
# tags then land 0.5 mm from the trimmed canvas edge.
fig1 <- (pA / (pB | pC) / pD / pE) +
  plot_layout(heights = c(0.76, 1, 0.72, 0.80)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.tag = element_text(family = FONT_TAG, size = PT_TAG, face = "bold",
                              colour = pal[["ink"]]),
      plot.tag.position = c(0, 1),
      plot.margin = margin(5, 3.5, 3.5, 5)
    )
  ) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

save_fig(fig1, "Figure1", width_mm = FIG_W_MM, height_mm = 216)
msg("Figure1 written: %s", out_dir("main"))
