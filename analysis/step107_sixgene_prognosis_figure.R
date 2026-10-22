#!/usr/bin/env Rscript
# ============================================================
# step107_sixgene_prognosis_figure.R
# Figure S20：GSE65682 六基因预后综合图（2 列 × 2 行 + 跨列森林图，A-E）
#
#   A 风险评分 KM 曲线（中位分组，log-rank）
#   B 六个 hub 基因各自的 KM 曲线（中位分组）
#   C 时间依赖 ROC（7/14/28 天）
#   D 28 天校准曲线（rms::calibrate，bootstrap 400）
#   E 协变量调整 Cox 森林图（风险评分 + 年龄 + 性别）
#
# 面板字母：可编辑黑色 Helvetica Bold 9 号；其余文字黑色 Helvetica。
# 输出：13_survival_prognosis/figures + 投稿包 figures/supplementary
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(timeROC)
  library(rms)
  library(ggplot2)
  library(cowplot)
  library(ragg)
  library(ggplotify)
  library(scales)
})

set.seed(20261015)
BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
PKG  <- file.path(BASE, "AJRCCM_submission/09_review_deAI_2026-09-25")
OUT  <- file.path(BASE, "13_survival_prognosis")
for (d in c("figures", "tables", "logs")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
LOGF <- file.path(OUT, "logs", "step107_run.log")
log_line <- function(...) { msg <- paste0(...); cat(msg, "\n"); cat(msg, "\n", file = LOGF, append = TRUE) }
log_line("=== step107 Figure S20 启动 ===")

HUB <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")
T0 <- 27.99
COL_HI <- "#E64B35FF"; COL_LO <- "#3C5488FF"

# ------------------------------------------------------------
# 数据：复用 step106 保存的分析数据
# ------------------------------------------------------------
d <- readRDS(file.path(OUT, "models", "step106_data.rds"))
df <- d$df
for (g in HUB) df[[paste0(g, "_z")]] <- as.numeric(scale(df[[g]]))
log_line(sprintf("n=%d，事件 %d", nrow(df), sum(df$event)))

theme_pub <- function(base_size = 8) theme_classic(base_size = base_size) +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(colour = "black"),
        legend.title = element_blank(),
        plot.title = element_text(size = 8, hjust = 0))

# ------------------------------------------------------------
# 模型：六基因多因素 Cox + 风险评分
# ------------------------------------------------------------
fml <- as.formula(paste("Surv(tte_adj, event) ~", paste(paste0(HUB, "_z"), collapse = " + ")))
mfit <- coxph(fml, data = df)
df$riskScore <- as.numeric(predict(mfit, type = "lp"))
df$riskGroup <- factor(ifelse(df$riskScore >= median(df$riskScore), "High", "Low"),
                       levels = c("Low", "High"))
lr <- survdiff(Surv(tte_adj, event) ~ riskGroup, data = df)
p_lr <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)
log_line(sprintf("log-rank P=%.4g", p_lr))

# (A) 风险评分 KM
km <- survfit(Surv(tte_adj, event) ~ riskGroup, data = df)
km_df <- data.frame(time = km$time, surv = km$surv, lower = km$lower, upper = km$upper,
                    group = rep(c("Low", "High"), km$strata))
s28 <- summary(km, times = T0)$surv
log_line(sprintf("28 天生存率：低风险 %.3f，高风险 %.3f", s28[1], s28[2]))
pA <- ggplot(km_df, aes(time, surv, colour = group, fill = group)) +
  geom_step(linewidth = 0.6) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  scale_colour_manual(values = c("Low" = COL_LO, "High" = COL_HI)) +
  scale_fill_manual(values = c("Low" = COL_LO, "High" = COL_HI)) +
  annotate("text", x = 21, y = 0.97, size = 2.6, family = "Helvetica", hjust = 1,
           label = sprintf("log-rank P=%.3f", p_lr)) +
  annotate("text", x = 21, y = 0.90, size = 2.3, family = "Helvetica", hjust = 1,
           label = sprintf("28-day survival: low %.2f, high %.2f", s28[1], s28[2])) +
  coord_cartesian(ylim = c(0.6, 1)) +
  labs(x = "Days since ICU admission", y = "Survival probability", title = "Six-gene risk score") +
  theme_pub() + theme(legend.position = c(0.25, 0.12), legend.text = element_text(size = 6.5),
                      legend.key.size = unit(3, "pt"), legend.background = element_blank())

# (B) 单基因 KM
km_list <- lapply(HUB, function(g) {
  grp <- factor(ifelse(df[[g]] >= median(df[[g]]), "High", "Low"), levels = c("Low", "High"))
  fit <- survfit(Surv(tte_adj, event) ~ grp, data = df)
  s <- summary(fit)
  p <- pchisq(survdiff(Surv(tte_adj, event) ~ grp, data = df)$chisq, df = 1, lower.tail = FALSE)
  data.frame(time = fit$time, surv = fit$surv,
             group = rep(c("Low", "High"), fit$strata), gene = g, P = p)
})
km_all <- do.call(rbind, km_list)
km_lab <- unique(km_all[, c("gene", "P")])
km_lab$lab <- sprintf("%s (P=%.2f)", km_lab$gene, km_lab$P)
km_all$gene_lab <- km_lab$lab[match(km_all$gene, km_lab$gene)]
km_all$gene_lab <- factor(km_all$gene_lab, levels = km_lab$lab)
pB <- ggplot(km_all, aes(time, surv, colour = group)) +
  geom_step(linewidth = 0.4) + facet_wrap(~ gene_lab, nrow = 2) +
  scale_colour_manual(values = c("Low" = COL_LO, "High" = COL_HI)) +
  coord_cartesian(ylim = c(0.55, 1)) +
  labs(x = "Days since ICU admission", y = "Survival", title = "Single genes (median split)") +
  theme_pub(7) + theme(strip.text = element_text(size = 5.5),
                       legend.position = "bottom", legend.text = element_text(size = 6),
                       legend.key.size = unit(3, "pt"))

# (C) 时间依赖 ROC
TT <- c(7, 14, T0); TL <- c(7, 14, 28)
tr <- timeROC(T = df$tte_adj, delta = df$event, marker = df$riskScore, cause = 1, times = TT, iid = TRUE)
A <- as.matrix(tr$AUC)
aucs <- A[, 1]
ci_tab <- as.matrix(confint(tr, T0)$CI_AUC) / 100
log_line(sprintf("AUC 7d=%.3f, 14d=%.3f, 28d=%.3f", aucs[1], aucs[2], aucs[3]))
roc_df <- do.call(rbind, lapply(seq_along(TT), function(i) {
  data.frame(FPR = 1 - tr$FP[, i], TPR = tr$TP[, i],
             time = factor(TL[i], levels = TL))
}))
pC <- ggplot(roc_df, aes(FPR, TPR, colour = time)) + geom_line(linewidth = 0.55) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  scale_colour_manual(values = c("#4DBBD5FF", "#00A087FF", "#E64B35FF"),
                      labels = sprintf("%d d (AUC=%.3f)", TL, round(aucs, 3))) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "Time-dependent ROC") +
  theme_pub() + theme(legend.position = c(0.65, 0.18), legend.text = element_text(size = 6),
                      legend.key.size = unit(3, "pt"), legend.background = element_blank())

# (D) 28 天校准（rms）
dd <- datadist(df); options(datadist = "dd")
cph_fit <- cph(as.formula(paste("Surv(tte_adj, event) ~", paste(paste0(HUB, "_z"), collapse = " + "))),
               data = df, x = TRUE, y = TRUE, surv = TRUE, time.inc = 28)
cal <- calibrate(cph_fit, u = 28, cmethod = "hare", m = 80, B = 400)
gD <- ggplotify::as.grob(function()
  plot(cal, xlab = "Predicted 28-day survival probability", ylab = "Observed",
       subtitles = FALSE, col = COL_HI, lwd = 2, cex.axis = 0.8, cex.lab = 0.9))
pD <- ggdraw(gD) + draw_label("28-day calibration", x = 0.02, y = 0.97, hjust = 0, vjust = 1,
                              size = 8, fontfamily = "Helvetica")

# (E) 协变量调整森林图（step105 结果）
forest_df <- read.csv(file.path(OUT, "tables", "adjusted_cox_forest.csv"), stringsAsFactors = FALSE)
forest_df$var <- c("Risk score (per unit)", "Age (per year)", "Female vs male")
forest_df$lab <- sprintf("%s  HR %.2f (95%% CI %.2f-%.2f); P=%.3g",
                         forest_df$var, forest_df$HR, forest_df$lo, forest_df$hi, forest_df$P)
forest_df$var_f <- factor(forest_df$var, levels = rev(forest_df$var))
pE <- ggplot(forest_df, aes(HR, var_f)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, colour = COL_HI, linewidth = 0.5) +
  geom_point(size = 2.4, shape = 16, colour = COL_HI) +
  geom_text(aes(x = max(hi) * 1.05, label = sprintf("HR %.2f (%.2f-%.2f), P=%.3g", HR, lo, hi, P)),
            hjust = 0, size = 2.6, family = "Helvetica") +
  scale_x_log10(limits = c(0.5, 8), expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = "Hazard ratio (log scale, 95% CI)", y = NULL,
       title = "Cox model adjusted for age and sex") +
  theme_pub() + theme(plot.title = element_text(size = 8, hjust = 0))

# ------------------------------------------------------------
# 组合：A B / C D / E（跨两列）
# ------------------------------------------------------------
top <- plot_grid(pA, pB, pC, pD, ncol = 2, nrow = 2,
                 labels = c("A", "B", "C", "D"),
                 label_fontfamily = "Helvetica", label_fontface = "bold",
                 label_size = 9, label_colour = "black",
                 label_x = 0.02, label_y = 0.98, hjust = 0, vjust = 1)
comp <- plot_grid(top, pE, ncol = 1, rel_heights = c(2.6, 1),
                  labels = c("", "E"),
                  label_fontfamily = "Helvetica", label_fontface = "bold",
                  label_size = 9, label_colour = "black",
                  label_x = 0.01, label_y = 0.98, hjust = 0, vjust = 1)
W <- 7.0; H <- 8.4
pdf(file.path(OUT, "figures", "FigureS20_sixgene_prognosis_GSE65682.pdf"), width = W, height = H, family = "Helvetica")
print(comp); dev.off()
agg_png(file.path(OUT, "figures", "FigureS20_sixgene_prognosis_GSE65682.png"), width = W, height = H,
        units = "in", res = 300, background = "white")
print(comp); dev.off()
agg_tiff(file.path(OUT, "figures", "FigureS20_sixgene_prognosis_GSE65682.tiff"), width = W, height = H,
         units = "in", res = 300, compression = "lzw", background = "white")
print(comp); dev.off()

SUPP <- file.path(PKG, "figures/supplementary")
dir.create(SUPP, recursive = TRUE, showWarnings = FALSE)
for (f in c("FigureS20_sixgene_prognosis_GSE65682.pdf", "FigureS20_sixgene_prognosis_GSE65682.png",
            "FigureS20_sixgene_prognosis_GSE65682.tiff"))
  file.copy(file.path(OUT, "figures", f), file.path(SUPP, f), overwrite = TRUE)

# 结果表（供补充表使用）
write.csv(data.frame(model = "Clinical + Hub6 six-gene score",
                     AUC_7d = round(aucs[1], 3), AUC_14d = round(aucs[2], 3), AUC_28d = round(aucs[3], 3),
                     AUC28_CI = sprintf("%.3f-%.3f", ci_tab[1, 1], ci_tab[1, 2]),
                     logrank_P = signif(p_lr, 3),
                     surv28_low = round(s28[1], 3), surv28_high = round(s28[2], 3)),
          file.path(OUT, "tables", "sixgene_prognosis_summary.csv"), row.names = FALSE)
log_line("=== step107 完成 ===")
