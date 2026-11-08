# ---------------------------------------------------------------------------
# 03_fig3.R -- Figure 3, discovery-internal screening and interaction network
#   A  cross-disease testing of the nominated candidates within the discovery
#      cohorts (ARDS list scored in sepsis, sepsis list scored in ARDS)
#   B  STRING protein-interaction network of the shared candidates
# Inputs : data/frozen/05_ML_SHAP/discovery_crossdisease_AUC.csv,
#          data/frozen/11_network_pharma/PPI_{nodes,edges}.csv
# Output : Figures/main/Figure3.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")
suppressPackageStartupMessages({ library(igraph); library(ggrepel) })

FZ <- function(...) file.path("data", "frozen", ...)

# --- panel A ---------------------------------------------------------------

cross_a <- function() {
  d <- read.csv(FZ("05_ML_SHAP", "discovery_crossdisease_AUC.csv"),
                stringsAsFactors = FALSE)
  names(d) <- c("gene", "list", "auc_ards", "p_ards", "auc_sepsis", "p_sepsis")
  d$list <- trimws(d$list)

  ards_list <- d[grepl("ARDS_candidate", d$list), ]
  ards_list$auc <- ards_list$auc_sepsis          # ARDS list scored in sepsis
  sepsis_list <- d[grepl("Sepsis_candidate", d$list), ]
  sepsis_list$auc <- sepsis_list$auc_ards        # sepsis list scored in ARDS

  classify <- function(x) {
    ifelse(x$gene %in% HUB6, "Confirmed hub gene",
           ifelse(x$auc < 0.5, "Direction-reversed (excluded)",
                  "Same direction (not retained)"))
  }

  one <- function(x, xlab) {
    x$class <- factor(classify(x), levels = c("Confirmed hub gene",
                                              "Direction-reversed (excluded)",
                                              "Same direction (not retained)"))
    x$gene <- factor(x$gene, levels = x$gene[order(x$auc)])
    ggplot(x, aes(auc, gene)) +
      geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3,
                 colour = pal[["slate"]]) +
      geom_segment(aes(x = 0.5, xend = auc, y = gene, yend = gene),
                   linewidth = 0.3, colour = pal[["slate"]]) +
      geom_point(aes(colour = class, shape = class, fill = class), size = 1.5,
                 stroke = 0.4) +
      scale_colour_manual(values = c("Confirmed hub gene" = pal[["brick"]],
                                     "Direction-reversed (excluded)" = pal[["steel"]],
                                     "Same direction (not retained)" = pal[["slate"]]),
                          name = NULL, drop = FALSE) +
      scale_fill_manual(values = c("Confirmed hub gene" = pal[["brick"]],
                                   "Direction-reversed (excluded)" = pal[["steel"]],
                                   "Same direction (not retained)" = "white"),
                        name = NULL, drop = FALSE) +
      scale_shape_manual(values = c("Confirmed hub gene" = 21,
                                    "Direction-reversed (excluded)" = 24,
                                    "Same direction (not retained)" = 21),
                         name = NULL, drop = FALSE) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                         expand = expansion(mult = c(0.02, 0.06))) +
      labs(x = xlab, y = NULL) +
      theme_ajr() +
      theme(legend.position = "right",
            legend.text = element_text(size = PT_BODY),
            legend.key.height = unit(6, "pt"),
            legend.key.width = unit(8, "pt"),
            axis.text.y = element_text(size = PT_BODY),
            plot.margin = margin(2, 2, 2, 2))
  }

  p1 <- one(ards_list, "Single-gene AUC in the sepsis cohort") +
    labs(subtitle = "ARDS-nominated candidates scored in sepsis") +
    theme(plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
  p2 <- one(sepsis_list, "Single-gene AUC in the ARDS cohort") +
    labs(subtitle = "Sepsis-nominated candidates scored in ARDS") +
    theme(plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))

  # One key for the whole panel, set under both sub-panels: the three classes
  # apply to either list, and no single sub-panel carries all three, so a key
  # inside one of them would document the wrong symbol set. The key is pulled
  # out as a grob and placed on its own row, which keeps both panels full
  # width and the axis labels legible.
  key <- cowplot::get_legend(
    p1 + theme(legend.position = "bottom",
               legend.direction = "horizontal",
               legend.spacing.x = unit(5, "pt"),
               legend.margin = margin(0, 0, 0, 0),
               legend.box.margin = margin(0, 0, 0, 0)))

  row <- cowplot::plot_grid(p1 + theme(legend.position = "none"),
                            p2 + theme(legend.position = "none"),
                            nrow = 1, align = "h", axis = "tb")

  cowplot::plot_grid(row, key, ncol = 1, rel_heights = c(1, 0.075))
}

# --- panel B ---------------------------------------------------------------

ppi_b <- function() {
  nodes <- read.csv(FZ("11_network_pharma", "PPI_nodes.csv"), stringsAsFactors = FALSE)
  edges <- read.csv(FZ("11_network_pharma", "PPI_edges.csv"), stringsAsFactors = FALSE)

  g <- graph_from_data_frame(edges[, c("preferredName_A", "preferredName_B")],
                             directed = FALSE, vertices = nodes$gene)
  set.seed(11)
  lay <- layout_with_fr(g, niter = 2000)
  xy <- data.frame(gene = V(g)$name, x = lay[, 1], y = lay[, 2])
  xy$degree <- nodes$degree[match(xy$gene, nodes$gene)]
  xy$is_hub <- xy$gene %in% HUB6

  module <- function(gene) {
    if (gene %in% c("UPP1", "UMPS", "CMPK1", "NME1", "NME2", "DTYMK"))
      return("Pyrimidine metabolism")
    if (grepl("^EXOSC|^DIS3", gene)) return("RNA exosome")
    if (gene %in% c("VNN1", "VNN2", "VNN3")) return("Pantetheinase")
    "Other"
  }
  xy$module <- ifelse(xy$is_hub, "Hub gene", vapply(xy$gene, module, ""))
  xy$module <- factor(xy$module, levels = c("Hub gene", "Pyrimidine metabolism",
                                            "RNA exosome", "Pantetheinase", "Other"))

  seg <- data.frame(
    x    = xy$x[match(edges$preferredName_A, xy$gene)],
    y    = xy$y[match(edges$preferredName_A, xy$gene)],
    xend = xy$x[match(edges$preferredName_B, xy$gene)],
    yend = xy$y[match(edges$preferredName_B, xy$gene)],
    w    = edges$score
  )

  # Labels are placed by rule, not by collision relaxation: each hub name is
  # pushed outward along the ray from the network centroid, which is where the
  # free space is, and joined to its node by a thin leader. Repulsion alone
  # kept dropping names back onto the dense pyrimidine module.
  hub_lab <- xy[xy$is_hub, ]
  ctr <- c(mean(xy$x), mean(xy$y))
  vx <- hub_lab$x - ctr[1]
  vy <- hub_lab$y - ctr[2]
  vn <- sqrt(vx^2 + vy^2)
  hub_lab$vx <- vx / vn
  hub_lab$vy <- vy / vn
  # how far along that ray each name has to travel to clear its own module;
  # UPP1 sits in the largest cluster and needs the longest reach
  clear <- c(UPP1 = 2.15, EXOSC4 = 1.70, CDK5RAP2 = 1.50, GRB10 = 1.45,
             CD177 = 1.40, VNN1 = 1.35)
  hub_lab$dist <- ifelse(is.na(clear[hub_lab$gene]), 1.35, clear[hub_lab$gene])
  hub_lab$lx <- hub_lab$x + hub_lab$dist * hub_lab$vx
  hub_lab$ly <- hub_lab$y + hub_lab$dist * hub_lab$vy
  hub_lab$hj <- ifelse(hub_lab$vx < -0.25, 1,
                       ifelse(hub_lab$vx > 0.25, 0, 0.5))

  mod_col <- c("Hub gene" = pal[["brick"]], "Pyrimidine metabolism" = pal[["steel"]],
               "RNA exosome" = pal[["sage"]], "Pantetheinase" = pal[["ochre"]],
               "Other" = "#C9CDD2")

  ggplot() +
    geom_segment(data = seg, aes(x, y, xend = xend, yend = yend,
                                 linewidth = w), colour = "#B9BEC4",
                 show.legend = FALSE) +
    scale_linewidth_continuous(range = c(0.12, 0.7)) +
    geom_point(data = xy, aes(x, y, size = degree, fill = module),
               shape = 21, colour = "white", stroke = 0.35) +
    geom_segment(data = hub_lab,
                 aes(x = x + 0.30 * vx, y = y + 0.30 * vy,
                     xend = lx, yend = ly),
                 linewidth = 0.2, colour = pal[["slate"]]) +
    geom_text(data = hub_lab, aes(lx, ly, label = gene, hjust = hj),
              family = FONT_BODY, fontface = "italic", size = 7 / .pt,
              colour = pal[["ink"]]) +
    scale_fill_manual(values = mod_col, name = NULL,
                      guide = guide_legend(override.aes = list(size = 1.6),
                                           order = 1)) +
    # a size key needs visible symbols: fill-driven points leave the key blank
    # on a white page, so the key is drawn with a solid grey glyph
    scale_size_continuous(range = c(0.6, 3.2), name = "STRING degree",
                          breaks = c(2, 6, 12),
                          guide = guide_legend(
                            override.aes = list(shape = 16, colour = "grey45"),
                            order = 2)) +
    coord_equal(clip = "off") +
    theme_void_ajr() +
    theme(legend.position = "right",
          legend.key.size = unit(6, "pt"),
          legend.spacing.y = unit(1, "mm"),
          plot.margin = margin(2, 4, 2, 2))
}

# --- assemble --------------------------------------------------------------

# Panel A is a row of two sub-panels under one key, so its letter is attached
# to the whole block; tagging the left sub-panel alone would tie the letter to
# half a panel and let the key drift away from it.
panel_a <- cowplot::ggdraw(cross_a()) +
  cowplot::draw_label("A", x = 0.004, y = 0.995, hjust = 0, vjust = 1,
                      fontfamily = FONT_TAG, size = PT_TAG, fontface = "bold",
                      colour = pal[["ink"]])

fig3 <- cowplot::plot_grid(panel_a,
                           tag_plot(ppi_b(), "B", x = 0.005, y = 1.005),
                           ncol = 1, rel_heights = c(0.92, 1.18))

fig3 <- fig3 + theme(plot.background = element_rect(fill = "white", colour = NA))

save_fig(fig3, "Figure3", width_mm = FIG_W_MM, height_mm = 200)
msg("Figure3 written")
