# ---------------------------------------------------------------------------
# 07_fig7.R -- Figure 7, Mendelian randomisation
#   A  forward MR (sepsis liability -> ARDS), four estimators
#   B  reverse MR (ARDS liability -> sepsis), four estimators
#   C  post hoc power as a function of the causal odds ratio
#   D  SMR causal estimates for the five hub genes with an eQTL instrument
# Input  : data/frozen/04_MR_SMR/{MR_sepsis4980_ARDS_results.csv,
#          MR_reverse_ards_sepsis_results.csv, MR_power_by_OR.csv,
#          SMR_coloc_summary.txt}
# Output : Figures/main/Figure7.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

FZ <- function(...) file.path("data", "frozen", "04_MR_SMR", ...)

method_lab <- c("Inverse variance weighted" = "Inverse variance weighted",
                "MR Egger" = "MR Egger",
                "Weighted median" = "Weighted median",
                "Weighted mode" = "Weighted mode")

forest <- function(file, title, n_iv) {
  d <- read.csv(FZ(file), stringsAsFactors = FALSE)
  d <- d[d$method %in% names(method_lab), ]
  d$method <- factor(method_lab[d$method], levels = rev(unname(method_lab)))
  d$OR <- exp(d$b)
  d$lo <- exp(d$b - 1.96 * d$se)
  d$hi <- exp(d$b + 1.96 * d$se)
  d$lab <- sprintf("%.2f (%.2f-%.2f)   P %s", d$OR, d$lo, d$hi,
                   fmt_p(d$pval, 2))

  xr <- range(c(d$lo, d$hi), na.rm = TRUE)
  ggplot(d, aes(OR, method)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16,
                   linewidth = 0.4, colour = pal[["ink"]]) +
    geom_point(shape = 22, size = 1.5, fill = pal[["brick"]],
               colour = pal[["ink"]], stroke = 0.3) +
    geom_text(aes(x = Inf, label = lab), hjust = 1.04, family = FONT_BODY,
              size = 7 / .pt, colour = pal[["ink"]]) +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4),
                  limits = c(xr[1] * 0.72, xr[2] * 1.02),
                  expand = expansion(mult = c(0.01, 0.95))) +
    labs(x = "Odds ratio (log scale)", y = NULL,
         subtitle = sprintf("%s, %d instruments", title, n_iv)) +
    theme_ajr() +
    theme(axis.text.y = element_text(size = PT_BODY),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["ink"]], face = "bold"))
}

power_c <- function() {
  d <- read.csv(FZ("MR_power_by_OR.csv"), stringsAsFactors = FALSE)
  d$dir <- ifelse(d$direction == "forward_sepsis_to_ARDS",
                  "Forward: sepsis to ARDS", "Reverse: ARDS to sepsis")

  ggplot(d, aes(OR, power_simulated, colour = dir)) +
    geom_hline(yintercept = 0.8, linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    annotate("text", x = min(d$OR), y = 0.82, label = "80% power", hjust = 0,
             family = FONT_BODY, size = 7 / .pt, colour = pal[["slate"]]) +
    geom_line(linewidth = 0.45) +
    geom_point(size = 0.9) +
    scale_colour_manual(values = c(pal[["brick"]], pal[["steel"]]), name = NULL) +
    scale_y_continuous(limits = c(0, 1.06), breaks = seq(0, 1, 0.2)) +
    labs(x = "Causal odds ratio", y = "Power") +
    theme_ajr() +
    theme(legend.position = "top", legend.direction = "horizontal",
          legend.margin = margin(0, 0, -3, 0))
}

smr_d <- function() {
  txt <- readLines(FZ("SMR_coloc_summary.txt"), warn = FALSE)
  i <- grep("^基因", txt)[1]
  rows <- txt[(i + 1):(i + 5)]
  parts <- strsplit(trimws(rows), "[[:space:]]+")
  d <- do.call(rbind, lapply(parts, function(p) {
    data.frame(gene = p[1],
               b_sepsis = as.numeric(p[3]), p_sepsis = as.numeric(p[4]),
               b_ards = as.numeric(p[5]), p_ards = as.numeric(p[6]),
               stringsAsFactors = FALSE)
  }))
  se_from_p <- function(b, p) abs(b) / abs(qnorm(p / 2))
  out <- rbind(
    data.frame(gene = d$gene, est = d$b_sepsis,
               lo = d$b_sepsis - 1.96 * se_from_p(d$b_sepsis, d$p_sepsis),
               hi = d$b_sepsis + 1.96 * se_from_p(d$b_sepsis, d$p_sepsis),
               grp = "Sepsis"),
    data.frame(gene = d$gene, est = d$b_ards,
               lo = d$b_ards - 1.96 * se_from_p(d$b_ards, d$p_ards),
               hi = d$b_ards + 1.96 * se_from_p(d$b_ards, d$p_ards),
               grp = "ARDS")
  )
  out$gene <- factor(out$gene, levels = rev(d$gene))
  out$grp <- factor(out$grp, levels = c("Sepsis", "ARDS"))

  ggplot(out, aes(est, gene, colour = grp)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3,
               colour = pal[["slate"]]) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16, linewidth = 0.4,
                   position = position_dodge(width = 0.5)) +
    geom_point(size = 1.1, position = position_dodge(width = 0.5)) +
    scale_colour_manual(values = c(pal[["steel"]], pal[["brick"]]), name = NULL) +
    labs(x = "SMR causal estimate (log OR per SD of expression)", y = NULL,
         subtitle = "No gene survives HEIDI or colocalisation") +
    theme_ajr() +
    theme(axis.text.y = element_text(face = "italic"),
          legend.position = "top", legend.direction = "horizontal",
          legend.margin = margin(0, 0, -3, 0),
          plot.subtitle = element_text(size = PT_BODY, hjust = 0,
                                       colour = pal[["slate"]]))
}

lab_args <- list(label_size = PT_TAG, label_fontfamily = FONT_TAG,
                 label_fontface = "bold", label_colour = pal[["ink"]])

pA <- forest("MR_sepsis4980_ARDS_results.csv",
             "Sepsis liability to ARDS", 18)
pB <- forest("MR_reverse_ards_sepsis_results.csv",
             "ARDS liability to sepsis", 7)

row1 <- do.call(cowplot::plot_grid,
                c(list(pA, pB, ncol = 2, labels = c("A", "B"), align = "h",
                       axis = "tb"), lab_args))
row2 <- do.call(cowplot::plot_grid,
                c(list(power_c(), smr_d(), ncol = 2, labels = c("C", "D"),
                       align = "h", axis = "tb"), lab_args))

fig7 <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1.05))

save_fig(fig7, "Figure7", width_mm = FIG_W_MM, height_mm = 200)
msg("Figure7 written")
