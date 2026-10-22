#!/usr/bin/env Rscript
# Figure 8 panels v2 (aligned with the original composite and manuscript text):
#   A: site-specific docking affinity heatmap (6 hub genes x 10 flavonoids, co-crystal site)
#   C: glucocorticoid repositioning (CMap / GEO perturbations)
# (Panel B PyMOL binding-pose renders are kept from the original pipeline.)
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

DD <- file.path(PROJ, "09_docking_drug")

## A: site-specific docking heatmap (co-crystal site; star marks co-crystallized pair)
ds <- read.csv(file.path(DD, "docking_site_results.csv"))
site_pref <- c("cocrystal", "co-crystal", "cocrystallized")
ds$is_cc <- tolower(ds$site) %in% site_pref | grepl("crystal", tolower(ds$site))
if (any(ds$is_cc)) {
  best <- ds %>% group_by(gene, ligand) %>%
    summarise(affinity = ifelse(any(is_cc), affinity[which(is_cc)[1]],
                                affinity[which.min(abs(affinity))]),
              cc = any(is_cc), .groups = "drop")
} else {
  best <- ds %>% group_by(gene, ligand) %>%
    summarise(affinity = min(affinity), cc = FALSE, .groups = "drop")
}
lig_order <- c("EGCG", "baicalein", "berberine", "curcumin", "emodin",
               "kaempferol", "luteolin", "quercetin", "resveratrol", "tetrahydroalstonine")
best$gene <- factor(best$gene, levels = rev(HUB6))
best$ligand <- factor(best$ligand, levels = lig_order[lig_order %in% best$ligand])
best$lab <- sprintf("%.1f", best$affinity)
best$lab[best$cc] <- paste0(best$lab[best$cc], " *")
pA <- ggplot(best, aes(ligand, gene, fill = affinity)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = lab), size = 2.0) +
  scale_fill_gradient(low = "#E64B35FF", high = "#3C5488FF",
                      name = "Affinity\n(kcal/mol)") +
  scale_y_discrete(labels = function(x) parse(text = gene_italic(x))) +
  labs(title = "A  Site-specific docking (AutoDock Vina)",
       subtitle = "* co-crystal binding site", x = NULL, y = NULL) +
  theme_ajrccm() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6),
        legend.title = element_text(size = 6), legend.key.height = unit(4, "mm"),
        plot.subtitle = element_text(size = 6))
save_pdf(pA, file.path(OUT_MAIN, "Fig8A_docking_site_heatmap.pdf"), 140, 62)

## C: glucocorticoid repositioning
gc <- data.frame(
  drug = factor(c("Hydrocortisone", "Dexamethasone", "Prednisolone"),
                levels = c("Hydrocortisone", "Dexamethasone", "Prednisolone")),
  neglogp = -log10(c(4.4e-03, 3.7e-04, 4.4e-05)),
  source = c("GEO (rank 4)", "CMap up (rank 1)", "GEO (rank 1)"),
  plab = c("P = 4.4e-03", "P = 3.7e-04", "P = 4.4e-05"))
pC <- ggplot(gc, aes(drug, neglogp)) +
  geom_col(fill = "#3C5488FF", width = 0.55) +
  geom_vline(xintercept = -log10(0.05) * 0 + 0, linewidth = 0) +
  geom_text(aes(label = plab, y = neglogp - 0.08), hjust = 1, size = 2.2, color = "white") +
  geom_text(aes(label = source, y = neglogp + 0.12), hjust = 0, size = 2.0, color = "grey30") +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 6.2), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "C  Glucocorticoid repositioning (CMap / GEO perturbations)",
       x = NULL, y = expression(-log[10](P))) +
  theme_ajrccm() + theme(plot.margin = margin(2, 12, 2, 2))
save_pdf(pC, file.path(OUT_MAIN, "Fig8C_glucocorticoid_repositioning.pdf"), 95, 52)

cat("Figure 8 v2 panels done.\n")
