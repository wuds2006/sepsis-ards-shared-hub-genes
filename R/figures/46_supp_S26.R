# ---------------------------------------------------------------------------
# 46_supp_S26.R -- Figure S26. ETosis and NADPH-oxidase program activity
#   A  ETosis AUCell score by monocyte subset, healthy control vs sepsis
#   B  Spearman correlation of the hub-gene programs with ETosis activity
#   C  donor S2 before/after treatment, by monocyte subset
#   D  donor-level ETosis score in PBMC monocytes, ARDS vs sepsis
# Input  : data/frozen/08_singlecell/etosis_aucell
# Output : Figures/supplementary/FigureS26.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

ET <- function(...) file.path("data", "frozen", "08_singlecell", "etosis_aucell", ...)

sc  <- read.csv(ET("GSE167363_ET_AUCell_scores.csv"), stringsAsFactors = FALSE,
                check.names = FALSE)
st  <- read.csv(ET("GSE167363_ET_stats.csv"), stringsAsFactors = FALSE)
cr  <- read.csv(ET("GSE167363_ET_correlations.csv"), stringsAsFactors = FALSE)
dm  <- read.csv(ET("GSE167363_ET_donor_means.csv"), stringsAsFactors = FALSE)
d2  <- read.csv(ET("GSE151263_ET_donor_means.csv"), stringsAsFactors = FALSE)
st2 <- read.csv(ET("GSE151263_ET_stats.csv"), stringsAsFactors = FALSE)

MONO <- c("ClassicalMono", "InflammatoryMono", "IntermediateMono",
          "NonclassicalMono", "ISGMono")
MONO_LAB <- c("Classical", "Inflammatory", "Intermediate", "Non-classical", "ISG+")

GRP_COL <- c(HC = pal[["sage"]], Sepsis = pal[["brick"]],
             ARDS = pal[["steel"]])

# subset colours: five lightness steps of the land/sea palette
CT_COL <- c(Classical = "#3E5C7E", Inflammatory = "#C99A5B",
            Intermediate = "#4C8A7E", `Non-classical` = "#8A7391",
            `ISG+` = "#A9BFD0")

sig_lab <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 2.2e-16, "P<2.2e-16", sprintf("P=%.1e", p)))
}

# --- panel A: ETosis by monocyte subset, HC vs sepsis ------------------------
panel_a <- function() {
  d <- sc[sc$celltype %in% MONO, ]
  d$celltype <- factor(d$celltype, levels = MONO, labels = MONO_LAB)
  set.seed(42)
  ds <- d %>% group_by(celltype, group) %>% slice_sample(n = 400) %>% ungroup()

  s <- st[st$geneset == "ETosis" & st$celltype %in% MONO, ]
  s$celltype <- factor(s$celltype, levels = MONO, labels = MONO_LAB)
  s$lab <- sig_lab(s$padj)
  s$y <- max(d$ETosis) * 1.02

  ggplot(d, aes(group, ETosis, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.62, linewidth = 0.3,
                 colour = pal[["ink"]], alpha = 0.85) +
    geom_point(data = ds, position = position_jitter(width = 0.18),
               size = 0.22, alpha = 0.22, colour = pal[["slate"]]) +
    geom_text(data = s, aes(x = 1.5, y = y, label = lab), inherit.aes = FALSE,
              size = PT_BODY / .pt, family = FONT_BODY, colour = pal[["ink"]]) +
    facet_wrap(~ celltype, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = GRP_COL, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = NULL, y = "ETosis program activity (AUCell)") +
    theme_ajr() +
    theme(strip.text = element_text(size = PT_BODY, face = "bold"),
          axis.text.x = element_text(size = PT_BODY),
          legend.position = "top", legend.margin = margin(0, 0, -3, 0),
          legend.key.size = unit(2.6, "mm"))
}

# --- panel B: hub programs versus ETosis activity ---------------------------
panel_b <- function() {
  b <- data.frame(program = c("Six-gene signature",
                              sub("_program$", "", cr$program)),
                  rho = c(cr$rho_hub[1], cr$rho_prog),
                  stringsAsFactors = FALSE)
  b$program <- factor(b$program, levels = b$program[order(b$rho)])
  b$grp <- ifelse(b$program == "Six-gene signature",
                  "Six-gene signature", "Per-hub program")

  ggplot(b, aes(rho, program, fill = grp)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = pal[["ink"]]) +
    geom_col(width = 0.64, colour = pal[["ink"]], linewidth = 0.25) +
    scale_fill_manual(values = c(`Six-gene signature` = pal[["brick"]],
                                 `Per-hub program` = pal[["steel"]]), name = NULL) +
    scale_x_continuous(limits = c(-0.45, 0.90),
                       breaks = c(-0.25, 0, 0.25, 0.50, 0.75)) +
    labs(x = "Spearman rho vs ETosis activity", y = NULL) +
    theme_ajr() +
    theme(legend.position = "none")
}

# --- panel C: donor S2 before and after treatment, by subset ----------------
# The two time points differ strongly in cell number, so a pooled weighted mean
# would be driven by composition rather than by a program shift. Panels are
# therefore drawn per subset with the cell count annotated. Positions are
# numeric so that the extra geom layer cannot demote the discrete axis.
panel_c <- function() {
  d <- dm[grepl("S2_T", dm$sample) & dm$celltype %in% MONO, ]
  d$tp <- ifelse(grepl("T0", d$sample), "T0", "T6")
  d$celltype <- factor(d$celltype, levels = MONO, labels = MONO_LAB)
  ord <- c("T0", "T6")
  d$xi <- match(d$tp, ord)

  ggplot(d, aes(xi, ETosis, group = celltype, colour = celltype)) +
    geom_line(linewidth = 0.5, alpha = 0.9) +
    geom_point(aes(size = sqrt(n_cells)), alpha = 0.9) +
    geom_text(aes(label = n_cells), vjust = -1.0, size = PT_BODY / .pt,
              family = FONT_BODY, show.legend = FALSE) +
    scale_x_continuous(breaks = 1:2, labels = ord, limits = c(0.85, 2.15)) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.16))) +
    scale_size_continuous(range = c(0.7, 3.0), guide = "none") +
    scale_colour_manual(values = CT_COL, name = NULL) +
    labs(x = "Donor S2 time point", y = "ETosis program activity (AUCell)") +
    theme_ajr() +
    theme(legend.position = "bottom", legend.margin = margin(-4, 0, 0, 0),
          legend.key.size = unit(2.4, "mm")) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE))
}

# --- panel D: PBMC monocytes, ARDS vs sepsis --------------------------------
panel_d <- function() {
  d <- d2[d2$celltype == "Mono", ]
  d$group <- factor(d$group, levels = c("ARDS", "Sepsis"))
  p <- st2$padj[st2$geneset == "ETosis" & st2$celltype == "Mono"]
  m1 <- mean(d$ETosis[d$group == "ARDS"])
  m2 <- mean(d$ETosis[d$group == "Sepsis"])

  ggplot(d, aes(group, ETosis, colour = group)) +
    geom_point(size = 2.0) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.45,
                 colour = pal[["ink"]], linewidth = 0.35) +
    annotate("text", x = 1.5, y = max(d$ETosis) * 1.04,
             label = sprintf("donor mean %.3f vs %.3f\ncell-level BH P = %.2f",
                             m1, m2, p),
             size = PT_BODY / .pt, family = FONT_BODY, lineheight = 0.95,
             colour = pal[["ink"]]) +
    scale_colour_manual(values = c(ARDS = pal[["brick"]], Sepsis = pal[["steel"]]),
                        name = NULL) +
    scale_y_continuous(limits = c(min(d$ETosis) * 0.98, max(d$ETosis) * 1.14),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = "ETosis program activity\n(donor mean, AUCell)") +
    theme_ajr() +
    theme(legend.position = "none")
}

# --- assemble ----------------------------------------------------------------
bot <- (panel_b() | panel_c() | panel_d()) +
  plot_layout(widths = c(1.0, 1.12, 0.85))

fig <- (panel_a() / bot) +
  plot_layout(heights = c(1, 1.06)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D"))) &
  theme(plot.tag = element_text(family = FONT_TAG, size = PT_TAG, face = "bold",
                                colour = pal[["ink"]]))

save_fig(fig, "FigureS26", width_mm = FIG_W_MM, height_mm = 158,
         dir = out_dir("supplementary"))
msg("FigureS26 written")
