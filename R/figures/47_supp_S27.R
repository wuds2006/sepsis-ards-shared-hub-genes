# ---------------------------------------------------------------------------
# 47_supp_S27.R -- Figure S27. Monocyte co-expression modules tested by hdWGCNA
#   A  scale-free topology fit across soft-threshold powers (GSE167363)
#   B  module assignment and module membership (kME) of the six hub genes
#   C  module eigengene versus sepsis status across metacells
#   D  module eigengene versus ETosis activity, donor level (descriptive)
# Input  : data/frozen/08_singlecell/hdwgcna
# Output : Figures/supplementary/FigureS27.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

HD <- function(...) file.path("data", "frozen", "08_singlecell", "hdwgcna", ...)

sp <- read.csv(HD("GSE167363_mono_softPowers.csv"), stringsAsFactors = FALSE)
mb <- read.csv(HD("GSE167363_mono_hub6_module_membership.csv"), stringsAsFactors = FALSE)
mt <- read.csv(HD("GSE167363_mono_module_trait.csv"), stringsAsFactors = FALSE)
ec <- read.csv(HD("GSE167363_mono_module_etosis_cor.csv"), stringsAsFactors = FALSE)

SEL_POWER <- 4          # first power with SFT.R.sq >= 0.80

# --- panel A: soft-threshold selection --------------------------------------
panel_a <- function() {
  sp$flag <- sp$Power == SEL_POWER
  r2 <- sp$SFT.R.sq[sp$flag]
  ggplot(sp, aes(Power, SFT.R.sq)) +
    geom_hline(yintercept = 0.80, linetype = "22", linewidth = 0.3,
               colour = pal[["slate"]]) +
    geom_vline(xintercept = SEL_POWER, linetype = "22", linewidth = 0.3,
               colour = pal[["brick"]]) +
    geom_line(linewidth = 0.35, colour = pal[["steel"]]) +
    geom_point(size = 1.5, colour = pal[["steel"]]) +
    annotate("text", x = SEL_POWER + 1.2, y = 0.42, hjust = 0,
             label = sprintf("selected power = %d (R\u00b2 = %.2f)", SEL_POWER, r2),
             size = PT_BODY / .pt, family = FONT_BODY, colour = pal[["ink"]]) +
    scale_x_continuous(breaks = c(1, seq(5, 30, by = 5))) +
    labs(x = "Soft-threshold power", y = "Scale-free fit (R\u00b2)") +
    theme_ajr()
}

# --- panel B: hub-gene module assignment ------------------------------------
panel_b <- function() {
  d <- mb
  d$cat <- ifelse(d$module == "not_variable", "not variable",
                  ifelse(d$module == "grey", "unassigned", d$module))
  d$disp <- ifelse(is.na(d$kME_own), 0, d$kME_own)
  d$gene <- factor(d$gene, levels = rev(d$gene))
  col <- c(turquoise = pal[["sage"]], unassigned = "#B8B8B8",
           `not variable` = "#E2E2E2")
  ggplot(d, aes(disp, gene, fill = cat)) +
    geom_col(width = 0.62, colour = pal[["ink"]], linewidth = 0.25) +
    geom_text(aes(label = ifelse(cat == "turquoise", sprintf("kME = %.2f", kME_own),
                         ifelse(cat == "unassigned", "grey", "n/a")),
                  x = disp + 0.02),
              hjust = 0, size = PT_BODY / .pt, family = FONT_BODY,
              colour = pal[["ink"]]) +
    scale_fill_manual(values = col, name = NULL,
                      breaks = c("not variable", "turquoise", "unassigned")) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_x_continuous(limits = c(0, 0.86), breaks = c(0, 0.25, 0.50, 0.75),
                       expand = expansion(mult = c(0, 0.01))) +
    labs(x = "Module membership (kME)", y = NULL) +
    theme_ajr() +
    theme(legend.position = "bottom", legend.margin = margin(-3, 0, 0, 0),
          legend.key.size = unit(2.4, "mm"),
          legend.spacing.x = unit(1.2, "mm"))
}

# --- panel C: module eigengene versus sepsis status -------------------------
# All five modules reach BH P < 0.05, so a significance fill map would carry no
# information; a single fill is used and significance is reported in the legend.
panel_c <- function() {
  d <- mt
  d$module <- factor(d$module, levels = d$module[order(d$cor_group)])
  ggplot(d, aes(cor_group, module)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = pal[["ink"]]) +
    geom_col(width = 0.62, fill = pal[["brick"]], colour = pal[["ink"]],
             linewidth = 0.25) +
    geom_text(aes(x = ifelse(cor_group >= 0, cor_group + 0.03, cor_group - 0.03),
                  label = sprintf("r = %.2f", cor_group),
                  hjust = ifelse(cor_group >= 0, 0, 1)),
              size = PT_BODY / .pt, family = FONT_BODY, colour = pal[["ink"]]) +
    scale_x_continuous(limits = c(min(d$cor_group) - 0.22,
                                  max(d$cor_group) + 0.18)) +
    labs(x = "Pearson r (module eigengene vs sepsis status)", y = NULL) +
    theme_ajr()
}

# --- panel D: module eigengene versus ETosis activity, donor level ----------
panel_d <- function() {
  d <- ec
  d$module <- factor(d$module, levels = d$module[order(d$rho_etosis)])
  ggplot(d, aes(rho_etosis, module)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = pal[["ink"]]) +
    geom_col(width = 0.62, fill = pal[["plum"]], colour = pal[["ink"]],
             linewidth = 0.25) +
    geom_text(aes(x = ifelse(rho_etosis >= 0, rho_etosis + 0.06, rho_etosis - 0.06),
                  label = sprintf("rho = %.1f", rho_etosis),
                  hjust = ifelse(rho_etosis >= 0, 0, 1)),
              size = PT_BODY / .pt, family = FONT_BODY, colour = pal[["ink"]]) +
    scale_x_continuous(limits = c(-1.45, 1.22),
                       breaks = c(-1.0, -0.5, 0, 0.5, 1.0)) +
    labs(x = "Spearman rho vs ETosis activity (donor level, n = 5)",
         y = NULL) +
    theme_ajr()
}

# --- assemble ----------------------------------------------------------------
fig <- cowplot::plot_grid(panel_a(), panel_b(), panel_c(), panel_d(),
                          ncol = 2, labels = c("A", "B", "C", "D"),
                          label_size = PT_TAG, label_fontfamily = FONT_TAG,
                          label_fontface = "bold", label_colour = pal[["ink"]])

save_fig(fig, "FigureS27", width_mm = FIG_W_MM, height_mm = 150,
         dir = out_dir("supplementary"))
msg("FigureS27 written")
