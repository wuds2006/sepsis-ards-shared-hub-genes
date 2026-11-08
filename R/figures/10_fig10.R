# ---------------------------------------------------------------------------
# 10_fig10.R -- Figure 10, protein-level validation and therapeutic relevance
#   A  HPA Blood Atlas protein change across ten infectious diseases
#   B  CMap connectivity results, two compound libraries
#   C  site-specific docking affinity, six targets x ten ligands
#   D  AlphaFold confidence over the whole model and over the docking box
#   E  PLIP contacts of the VNN1-luteolin complex
# Input  : data/frozen/10_HPA_validation, data/frozen/09_docking_drug
# Output : Figures/main/Figure10.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

HP <- function(...) file.path("data", "frozen", "10_HPA_validation", ...)
DK <- function(...) file.path("data", "frozen", "09_docking_drug", ...)

GENE6 <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")

# --- panel A ---------------------------------------------------------------

panel_a <- function() {
  d <- read.csv(HP("HPA_CD177_VNN2_infection_vs_healthy.csv"),
                stringsAsFactors = FALSE)
  d$Gene <- factor(d$Gene, levels = c("CD177", "VNN2"))
  ord <- d$Disease[d$Gene == "CD177"][order(d$logFC[d$Gene == "CD177"])]
  d$Disease <- factor(d$Disease, levels = ord)

  ggplot(d, aes(logFC, Disease, colour = Gene)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    geom_point(aes(size = neglogp), shape = 16,
               position = position_dodge(width = 0.55)) +
    scale_colour_manual(values = c(CD177 = pal[["brick"]], VNN2 = pal[["steel"]]),
                        name = NULL) +
    scale_size_area(max_size = 2.4, breaks = c(1, 3),
                    name = "-log10 adj. P") +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, 26)) +
    labs(x = "Protein log2 fold change (infection vs healthy)", y = NULL) +
    theme_ajr() +
    theme(legend.position = "top", legend.box = "horizontal",
          legend.margin = margin(0, 0, -5, 0), legend.spacing.x = unit(2.5, "mm"),
          legend.key.size = unit(2.6, "mm"),
          axis.text.y = element_text(lineheight = 0.86),
          plot.margin = margin(1, 2, 1, 4))
}

# --- panel B ---------------------------------------------------------------

panel_b <- function() {
  cm <- read.csv(DK("cmap_hub6_realtime_query.csv"), stringsAsFactors = FALSE)
  cm <- cm[cm$library == "Old_CMAP_up", ]
  cm <- cm[order(cm$rank), c("term", "pvalue")]
  cm$set <- "CMap"
  ge <- read.csv(DK("Drug_Perturbations_from_GEO_2014_hub_prediction.csv"),
                 stringsAsFactors = FALSE)
  ge <- ge[order(ge$pvalue), c("term", "pvalue")]
  ge$set <- "GEO"
  d <- rbind(cm, ge)
  d$lab <- sub("-\\d+$", "", d$term)
  d$lab <- sub(" homo sapiens.*$| mus musculus.*$| rattus norvegicus.*$", "", d$lab)
  d$lab <- paste0(toupper(substr(d$lab, 1, 1)), substr(d$lab, 2, nchar(d$lab)))
  d$nlp <- -log10(d$pvalue)
  d <- d[order(d$set, d$pvalue), ]
  d <- d[!duplicated(paste(d$set, d$lab)), ]
  d <- do.call(rbind, lapply(split(d, d$set), function(x) head(x, 8)))
  d$set <- factor(d$set, levels = c("CMap", "GEO"))

  # Eight CMap compounds are exactly tied at the top (P = 3.664e-4, z = 101.51,
  # combined score 803.12), so "rank" is an arbitrary ordering inside the tie.
  # Two consequences are handled here: (i) within the CMap facet the highest
  # ranked compound must sit at the TOP, otherwise rank 1 is drawn at the bottom
  # and reads as last place; (ii) the eight identical -log10 P labels are
  # collapsed to a single "(8 tied)" annotation instead of eight repeats of 3.4.
  d$idx <- stats::ave(seq_len(nrow(d)), d$set, FUN = seq_along)
  lv <- unlist(lapply(c("CMap", "GEO"), function(s) {
    x <- d[d$set == s, ]
    x$lab[if (s == "CMap") order(-x$idx) else order(x$nlp)]
  }))
  d$lab <- factor(d$lab, levels = lv)

  d$tie <- stats::ave(d$nlp, d$set, FUN = function(v)
    vapply(v, function(z) sum(abs(v - z) < 1e-9), numeric(1)))
  # One label per distinct value: eight identical "3.4"s would read as noise.
  # The tie itself is stated in the legend, not repeated eight times.
  d$labtxt <- ifelse(duplicated(paste(d$set, round(d$nlp, 3))), "",
                     sprintf("%.1f", d$nlp))

  ggplot(d, aes(nlp, lab, fill = set)) +
    geom_col(width = 0.66, colour = NA) +
    geom_text(data = d[d$labtxt != "", ], aes(label = labtxt), hjust = -0.15,
              family = FONT_BODY, size = 7 / .pt, colour = pal[["ink"]]) +
    facet_wrap(~ set, nrow = 1, scales = "free") +
    scale_fill_manual(values = c(CMap = pal[["slate"]],
                                 GEO = pal[["plum"]]),
                      guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.30))) +
    labs(x = "-log10 P", y = NULL) +
    theme_ajr() +
    theme(strip.text = element_text(size = PT_BODY, hjust = 0),
          panel.spacing.x = unit(5, "mm"),
          axis.text.y = element_text(lineheight = 0.86))
}

# --- panel C ---------------------------------------------------------------

panel_c <- function() {
  d <- read.csv(DK("docking_site_results.csv"), stringsAsFactors = FALSE)
  d$gene <- factor(d$gene, levels = rev(GENE6))
  lig <- unique(d$ligand)
  d$ligand <- factor(d$ligand, levels = lig)
  best <- d[which.min(d$affinity), ]

  ggplot(d, aes(ligand, gene, fill = affinity)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_tile(data = best, colour = pal[["ink"]], linewidth = 0.5, fill = NA) +
    scale_fill_gradientn(colours = rev(ramp_seq), name = "Vina\n(kcal/mol)") +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, 12)) +
    scale_y_discrete(labels = function(x) ifelse(x %in% HUB6, x, x)) +
    labs(x = NULL, y = NULL) +
    theme_ajr() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1,
                                     lineheight = 0.86),
          axis.text.y = element_text(face = "italic"),
          axis.line = element_blank(), axis.ticks = element_blank(),
          legend.key.height = unit(4, "mm"), legend.key.width = unit(2.2, "mm"),
          legend.title = element_text(size = PT_BODY))
}

# --- panel D ---------------------------------------------------------------

panel_d <- function() {
  d <- read.csv(DK("alphafold_plddt_assessment.csv"), stringsAsFactors = FALSE)
  d <- d[d$gene %in% GENE6, ]
  d$gene <- factor(d$gene, levels = rev(GENE6[GENE6 %in% d$gene]))

  ggplot(d, aes(y = gene)) +
    geom_vline(xintercept = 70, linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    geom_segment(aes(x = mean_pLDDT, xend = box_mean_pLDDT, yend = gene),
                 linewidth = 0.5, colour = "#9AA0A6") +
    geom_point(aes(x = mean_pLDDT), shape = 21, size = 1.5, stroke = 0.35,
               fill = "white", colour = pal[["ink"]]) +
    geom_point(aes(x = box_mean_pLDDT), shape = 16, size = 1.6,
               colour = pal[["brick"]]) +
    scale_x_continuous(limits = c(50, 100), breaks = seq(50, 100, 10)) +
    labs(x = "AlphaFold pLDDT", y = NULL) +
    theme_ajr() +
    theme(axis.text.y = element_text(face = "italic"))
}

# --- panel E ---------------------------------------------------------------

panel_e <- function() {
  d <- read.csv(DK("plip_vnn1_luteolin.csv"), stringsAsFactors = FALSE)
  d$type <- factor(d$type, levels = c("hydrogen bond", "hydrophobic", "pi-stacking"))
  d$res <- factor(d$residue, levels = unique(d$residue[order(d$distance)]))
  d$side[d$side == ""] <- "contact"

  ggplot(d, aes(distance, res)) +
    geom_segment(aes(x = 2.4, xend = distance, yend = res), linewidth = 0.3,
                 colour = "#C4C9CE") +
    geom_point(aes(colour = type, shape = type), size = 1.7) +
    scale_colour_manual(values = c(`hydrogen bond` = pal[["brick"]],
                                   hydrophobic = pal[["olive"]],
                                   `pi-stacking` = pal[["steel"]]), name = NULL) +
    scale_shape_manual(values = c(16, 15, 17), name = NULL) +
    scale_x_continuous(limits = c(2.4, 5.0), breaks = seq(2.5, 5, 0.5),
                       expand = expansion(mult = c(0.01, 0.04))) +
    labs(x = "Contact distance (A)", y = NULL) +
    theme_ajr() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -5, 0),
          legend.key.size = unit(2.8, "mm"), legend.spacing.x = unit(2.5, "mm"))
}

# --- assemble --------------------------------------------------------------

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

row1 <- do.call(cowplot::plot_grid,
                c(list(panel_a(), panel_b(), ncol = 2,
                       rel_widths = c(0.48, 0.52), labels = c("A", "B")), lab_args))
row2 <- do.call(cowplot::plot_grid,
                c(list(panel_c(), panel_d(), ncol = 2,
                       rel_widths = c(0.62, 0.38), labels = c("C", "D")), lab_args))
row3 <- do.call(cowplot::plot_grid,
                c(list(panel_e(), labels = "E"), lab_args))

fig10 <- cowplot::plot_grid(row1, row2, row3, ncol = 1,
                            rel_heights = c(1.10, 0.95, 0.60))

save_fig(fig10, "Figure10", width_mm = FIG_W_MM, height_mm = 206)
msg("Figure10 written")
