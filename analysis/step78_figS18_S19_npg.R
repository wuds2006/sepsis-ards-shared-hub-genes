#!/usr/bin/env Rscript
# Step 78 — Figure S18 / S19 以 ggplot2 + ggsci NPG 重绘（Helvetica，AJRCCM 发表规格）
# S18: CMap 糖皮质激素富集三背景阴性对照（A: best -log10P 分布；B: top-10 命中率）
# S19: 回对接逐模式 RMSD 验证（A: UPP1 3EUF/BAU；B: VNN1 7SLV/9TC）
# 输出: PDF 矢量 + PNG @300 dpi + TIFF @1000 dpi (ragg, LZW, RGB)

source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ragg) })

FIGDIR <- file.path(PROJ, "AJRCCM_submission/05_supplementary/Supplementary_Figures_Tables/Figures")
REAL_P <- 3.66e-4

# NPG 配色映射（与图内各组对应）
C_REAL  <- "#E64B35FF"   # 真实六基因签名 = NPG 红
C_CAND  <- "#8491B4FF"   # 非 hub 候选背景 = NPG 浅蓝
C_DEG   <- "#F39B7FFF"   # DEG-matched 背景 = NPG 浅橙
C_UNIV  <- "#4DBBD5FF"   # 全转录组背景 = NPG 青蓝
C_TOP   <- "#3C5488FF"   # 最优模式 = NPG 深蓝

save_all <- function(plot, out_base, w_mm, h_mm) {
  save_pdf(plot, paste0(out_base, ".pdf"), w_mm, h_mm)
  agg_png(paste0(out_base, ".png"), width = w_mm / 25.4, height = h_mm / 25.4,
          units = "in", res = 300, background = "white")
  print(plot); dev.off()
  agg_tiff(paste0(out_base, ".tiff"), width = w_mm / 25.4, height = h_mm / 25.4,
           units = "in", res = 1000, background = "white", compression = "lzw")
  print(plot); dev.off()
  cat("saved:", out_base, "\n")
}

## ---------------- Figure S18 ----------------
neg <- read_csv(file.path(PROJ, "09_docking_drug/cmap_negative_control_results.csv"),
                show_col_types = FALSE)
deg <- read_csv(file.path(PROJ, "09_docking_drug/cmap_deg_matched_null.csv"),
                show_col_types = FALSE)

cand <- neg %>% filter(background == "A_nonhub_candidates", library == "Old_CMAP_up")
univ <- neg %>% filter(background == "B_measured_universe", library == "Old_CMAP_up")

neglog <- function(p) suppressWarnings(-log10(as.numeric(p)))
top10 <- function(x) sum(as.logical(x), na.rm = TRUE)
cand_vals <- neglog(cand$best_gluco_p); cand_vals <- cand_vals[!is.na(cand_vals)]
deg_vals  <- neglog(deg$best_gluco_p);  deg_vals  <- deg_vals[!is.na(deg_vals)]
real_x    <- -log10(REAL_P)

cat(sprintf("S18-A: cand hits %d/%d top10 %d | deg hits %d/%d top10 %d | univ top10 %d | real %.3f\n",
            length(cand_vals), nrow(cand), top10(cand$gluco_in_top10),
            length(deg_vals), nrow(deg), top10(deg$gluco_in_top10),
            top10(univ$gluco_in_top10), real_x))

hist_df <- bind_rows(
  tibble(v = cand_vals, bg = paste0("non-hub candidates (n=", length(cand_vals), " with hits)")),
  tibble(v = deg_vals,  bg = paste0("DEG-matched (n=", length(deg_vals), " with hits)"))
)
hist_df$bg <- factor(hist_df$bg, levels = c(
  paste0("non-hub candidates (n=", length(cand_vals), " with hits)"),
  paste0("DEG-matched (n=", length(deg_vals), " with hits)")))

p18A <- ggplot(hist_df, aes(v, fill = bg)) +
  geom_histogram(binwidth = 0.1, boundary = 1.0, alpha = 0.8,
                 color = "white", linewidth = 0.15, position = "identity") +
  geom_vline(xintercept = real_x, color = C_REAL, linewidth = 0.5) +
  annotate("text", x = real_x + 0.08, y = 72, hjust = 0, size = 2.3,
           label = "six-hub-gene\nsignature", color = C_REAL, family = AJRCCM_FONT) +
  scale_fill_manual(values = setNames(c(C_CAND, C_DEG), levels(hist_df$bg)), name = NULL) +
  scale_x_continuous(limits = c(1, 6), breaks = 1:6) +
  labs(title = "A  Null distributions (Old_CMAP_up)",
       x = expression("best glucocorticoid " * -log[10] * "(P)"),
       y = "random gene sets") +
  theme_ajrccm() +
  theme(legend.position = c(0.62, 0.72),
        legend.background = element_rect(fill = "white", color = NA))

freq <- tibble(
  bg = c("random\nmeasured genes", "random non-hub\ncandidates",
         "DEG-matched\nrandom sets", "six-hub-gene\nsignature"),
  pct = c(100 * top10(univ$gluco_in_top10) / nrow(univ),
          100 * top10(cand$gluco_in_top10) / nrow(cand),
          100 * top10(deg$gluco_in_top10) / nrow(deg),
          100)
)
freq$bg <- factor(freq$bg, levels = freq$bg)

p18B <- ggplot(freq, aes(bg, pct, fill = bg)) +
  geom_col(width = 0.62, color = "black", linewidth = 0.25) +
  geom_text(aes(label = paste0(round(pct), "%")), vjust = -0.5, size = 2.3,
            family = AJRCCM_FONT) +
  scale_fill_manual(values = c(C_UNIV, C_CAND, C_DEG, C_REAL), guide = "none") +
  scale_y_continuous(limits = c(0, 112), expand = c(0, 0)) +
  labs(title = "B  Top-10 hit frequency", x = NULL,
       y = "glucocorticoid in top 10 (%)") +
  theme_ajrccm() +
  theme(axis.text.x = element_text(size = 6))

p18 <- p18A + p18B + plot_layout(widths = c(1.05, 1))
save_all(p18, file.path(FIGDIR, "FigureS18_CMap_null_distribution"), 180, 80)

## ---------------- Figure S19 ----------------
rd <- read_csv(file.path(PROJ, "09_docking_drug/redocking/redocking_per_mode.csv"),
               show_col_types = FALSE) %>%
  mutate(mode_lab = paste0("Mode ", mode),
         panel = ifelse(target == "UPP1",
                        "A  UPP1 (PDB 3EUF) - BAU",
                        "B  VNN1 (PDB 7SLV) - 9TC"),
         fill = ifelse(top_ranked, "top", "other"))
rd$panel <- factor(rd$panel, levels = unique(rd$panel))

p19 <- ggplot(rd, aes(mode_lab, rmsd_A, fill = fill)) +
  geom_col(width = 0.62, color = "black", linewidth = 0.25) +
  geom_hline(yintercept = 2.0, linetype = "dashed", color = C_REAL, linewidth = 0.45) +
  geom_text(aes(y = rmsd_A + 0.28, label = sprintf("%.2f", rmsd_A)),
            size = 2.1, family = AJRCCM_FONT) +
  geom_text(aes(y = -1.15, label = sprintf("%.2f", affinity_kcal_mol)),
            size = 1.9, color = "grey35", family = AJRCCM_FONT) +
  geom_text(data = ~ subset(.x, mode == max(mode)),
            aes(x = 5.42, y = 2.0, label = "2.0 \u00c5"), color = C_REAL,
            size = 2.1, hjust = 1, vjust = -0.35, family = AJRCCM_FONT) +
  geom_text(data = ~ subset(.x, mode == 1),
            aes(x = 0.62, y = 11.2, label = "top-ranked pose"), color = C_TOP,
            size = 1.9, hjust = 0, family = AJRCCM_FONT) +
  geom_text(data = ~ subset(.x, mode == 1),
            aes(x = 0.62, y = -1.78, label = "Vina affinity (kcal/mol) below each mode"),
            size = 1.9, hjust = 0, color = "grey35", family = AJRCCM_FONT) +
  scale_fill_manual(values = c(top = C_TOP, other = C_DEG), guide = "none") +
  scale_y_continuous(limits = c(-1.9, 11.6), breaks = seq(0, 10, 2), expand = c(0, 0)) +
  facet_wrap(~panel, nrow = 1) +
  labs(x = NULL, y = "Heavy-atom RMSD (\u00c5)") +
  coord_cartesian(clip = "off") +
  theme_ajrccm() +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0),
        plot.margin = margin(2, 12, 2, 2))

save_all(p19, file.path(FIGDIR, "FigureS19_Redocking_validation"), 180, 78)

cat("step78 done.\n")
