#!/usr/bin/env Rscript
# step105b_pergene_KM_GSE65682.R
# 6 个 hub 基因各自的 KM 曲线（中位表达分组，28 天死亡）
suppressPackageStartupMessages({library(survival); library(survminer); library(ggplot2); library(patchwork)})

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "13_survival_prognosis")
df   <- read.csv(file.path(OUT, "tables", "GSE65682_risk_scores.csv"))
HUB  <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")
COL_HI <- "#E64B35FF"; COL_LO <- "#3C5488FF"

plist <- list()
ptab <- data.frame(gene = HUB, logrank_P = NA)
for (i in seq_along(HUB)) {
  g <- HUB[i]
  dd <- df
  dd$grp <- factor(ifelse(dd[[g]] >= median(dd[[g]]), "High", "Low"), levels = c("Low", "High"))
  fit <- survfit(Surv(tte, event) ~ grp, data = dd)
  lr <- survdiff(Surv(tte, event) ~ grp, data = dd)
  pv <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)
  ptab$logrank_P[i] <- pv
  # 用 survminer 核心但禁图例，子图统一精简
  p <- ggsurvplot(fit, data = dd, pval = TRUE, conf.int = FALSE,
                  palette = c(COL_LO, COL_HI), legend = "none",
                  xlab = "Days since ICU admission", ylab = "Survival probability",
                  ggtheme = theme_classic(base_size = 7) +
                            theme(plot.title = element_text(face = "bold.italic", size = 8)),
                  title = g)
  # ggsurvplot 返回对象，取 $plot
  plist[[g]] <- p$plot + theme(legend.position = "none")
}
write.csv(ptab, file.path(OUT, "tables", "pergene_KM_logrank.csv"), row.names = FALSE)

comp <- wrap_plots(plist, ncol = 3)
pdf(file.path(OUT, "figures", "KM_pergene_6hub.pdf"), width = 7.0, height = 4.7, family = "Helvetica")
print(comp); dev.off()
png(file.path(OUT, "figures", "KM_pergene_6hub.png"), width = 7.0, height = 4.7, units = "in", res = 300, family = "Helvetica")
print(comp); dev.off()
print(ptab)
cat("完成：pergene KM\n")
