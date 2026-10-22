#!/usr/bin/env Rscript
# step45: recompute six-hub-gene signature per-sample score with the GSVA method
# (gsvaParam, Hänzelmann 2013) in the two discovery cohorts (GSE76293, GSE69528),
# replacing the previous ssGSEA implementation to match the revised Methods.
# Also re-renders the Figure 5A panel (GSVA score, ggplot2 + ggsci NPG).
# Outputs (07_hub_immune/):
#   <ds>_hub6_GSVA.csv   per-sample GSVA score of the six-gene signature
#   hub6_GSVA_stats.csv  Wilcoxon summary (discovery cohorts)
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
suppressPackageStartupMessages({ library(GSVA); library(dplyr) })

BASE <- PROJ
IMM  <- file.path(BASE, "07_hub_immune")
OUT_MAIN <- file.path(BASE, "AJRCCM_submission/03_figures/Figures_ggplot2_NPG/main_panels")
HUB6 <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")

DS <- list(
  list(id = "GSE76293", dis = "ARDS",   tag = "ARDS discovery (PMN fractions)"),
  list(id = "GSE69528", dis = "Sepsis", tag = "Sepsis discovery (whole blood)")
)

stats <- data.frame()
ss <- do.call(rbind, lapply(DS, function(z) {
  expr <- as.matrix(readRDS(file.path(BASE, "01_DEG", paste0(z$id, "_gene_expr.rds"))))
  grp  <- readRDS(file.path(BASE, "01_DEG", paste0(z$id, "_grp.rds")))
  if (is.data.frame(grp)) grp <- grp$group
  genes <- intersect(HUB6, rownames(expr))
  sc <- gsva(gsvaParam(exprData = expr, geneSets = list(hub6 = genes),
                       kcdf = "Gaussian", minSize = 2))
  score <- as.numeric(sc[1, ])
  df <- data.frame(sample = colnames(expr), group = grp, score = score)
  write.csv(df, file.path(IMM, paste0(z$id, "_hub6_GSVA.csv")), row.names = FALSE)
  x  <- score[grp == z$dis]; y0 <- score[grp == "Control"]
  p  <- wilcox.test(x, y0, exact = FALSE, correct = FALSE)$p.value
  stats <<- rbind(stats, data.frame(dataset = z$id, disease = z$dis,
                                    genes_used = paste(genes, collapse = ","),
                                    n_disease = length(x), n_control = length(y0),
                                    median_disease = median(x), median_control = median(y0),
                                    P = p))
  cat(sprintf("%s hub6 GSVA: median %.3f vs %.3f, P = %.3g\n",
              z$id, median(x), median(y0), p))
  data.frame(dataset = paste0(z$id, " ", z$dis, " discovery"),
             group = factor(df$group, levels = c("Control", "ARDS", "Sepsis")),
             score = df$score,
             plab = sprintf("italic(P) == %.2g", p))
}))
write.csv(stats, file.path(IMM, "hub6_GSVA_stats.csv"), row.names = FALSE)

ss$dataset <- factor(ss$dataset, levels = unique(ss$dataset))
plab_df <- ss %>% distinct(dataset, plab) %>%
  mutate(y = max(ss$score) * 1.06, x = 1.5)

p5a <- ggplot(ss, aes(group, score, fill = group)) +
  geom_violin(color = "grey30", linewidth = 0.25, scale = "width", width = 0.8) +
  geom_boxplot(width = 0.14, linewidth = 0.25, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.05, size = 0.35, alpha = 0.5, color = "grey25") +
  geom_text(data = plab_df, aes(x, y, label = plab), inherit.aes = FALSE,
            parse = TRUE, size = 2.2) +
  facet_wrap(~dataset, nrow = 1, scales = "free_x") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  scale_fill_manual(values = c(Control = COL_CONTROL, ARDS = COL_DISEASE, Sepsis = COL_DISEASE)) +
  labs(title = "Six-hub-gene signature GSVA score", x = NULL,
       y = "GSVA score", fill = NULL) +
  theme_ajrccm() +
  theme(legend.position = "none", strip.text = element_text(size = 6))
save_pdf(p5a, file.path(OUT_MAIN, "Fig5A_hub6_GSVA.pdf"), 178, 50)
cat("step45 done.\n")
