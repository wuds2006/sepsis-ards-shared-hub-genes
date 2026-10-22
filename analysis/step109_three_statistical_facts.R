#!/usr/bin/env Rscript
# ============================================================
# step109_three_statistical_facts.R
# 补齐三项统计事实（供意见 3 使用）
#   ① GSE65682 Cox 模型的比例风险（PH）假设检验（cox.zph）
#   ② KM 中位切点之外是否有第 2 切点或连续项敏感性分析
#   ③ 免疫 22 亚型 BH 覆盖面（另见 step109b，python 侧）
# 输出：13_survival_prognosis/tables/step109_*.csv + logs/step109.log
# ============================================================

suppressPackageStartupMessages({
  library(survival)
})

set.seed(20261019)
BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "13_survival_prognosis")
LOGF <- file.path(OUT, "logs", "step109_run.log")
log_line <- function(...) { msg <- paste0(...); cat(msg, "\n"); cat(msg, "\n", file = LOGF, append = TRUE) }
log_line("=== step109 启动 ===")

dat <- readRDS(file.path(OUT, "models/step106_data.rds"))
df  <- dat$df
HUB <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")

# z 标化（与 step105/106 口径一致）
Z <- scale(df[, HUB])
colnames(Z) <- paste0("z_", HUB)
df2 <- cbind(df, as.data.frame(Z))

f_hub  <- as.formula(paste("Surv(tte_adj, event) ~", paste(colnames(Z), collapse = " + ")))
f_full <- as.formula(paste("Surv(tte_adj, event) ~ age + sex + endotype +", paste(colnames(Z), collapse = " + ")))

## ---------- ① PH 假设检验 ----------
log_line("")
log_line("### ① 比例风险（PH）假设检验 ###")
m_hub  <- coxph(f_hub,  data = df2, ties = "efron")
m_full <- coxph(f_full, data = df2, ties = "efron")

zp_hub  <- cox.zph(m_hub)
zp_full <- cox.zph(m_full)
cat("\n--- Hub6 Cox: cox.zph ---\n"); print(zp_hub)
cat("\n--- Hub6 + Clinical Cox: cox.zph ---\n"); print(zp_full)

write.csv(as.data.frame(zp_hub$table)  |> (\(d){ d$term <- rownames(d); d })(),
          file.path(OUT, "tables", "step109_coxzph_hub6.csv"), row.names = FALSE)
write.csv(as.data.frame(zp_full$table) |> (\(d){ d$term <- rownames(d); d })(),
          file.path(OUT, "tables", "step109_coxzph_full.csv"), row.names = FALSE)

pmin_hub  <- min(zp_hub$table[, "p"])
pmin_full <- min(zp_full$table[, "p"])
pglobal_hub  <- zp_hub$GLOBAL <- zp_hub$table[nrow(zp_hub$table), "p"]
log_line(sprintf("Hub6 模型 全局 PH P = %.4f；最小单项 P = %.4f", pglobal_hub, pmin_hub))
log_line(sprintf("Hub6+临床 模型 全局 PH P = %.4f；最小单项 P = %.4f",
                 zp_full$table[nrow(zp_full$table), "p"], pmin_full))

## ---------- ② KM 切点敏感性 ----------
log_line("")
log_line("### ② KM 切点 / 连续项敏感性分析 ###")

# risk score = Hub6 Cox 线性预测值（与 Figure S20 口径一致）
df2$risk <- as.numeric(predict(m_hub, type = "lp"))

# (a) 连续项 Cox
m_cont <- coxph(Surv(tte_adj, event) ~ risk, data = df2, ties = "efron")
s_cont <- summary(m_cont)
cat("\n--- 连续 risk score Cox ---\n"); print(s_cont$coefficients)

# (b) 不同切点方案的 log-rank
cuts <- list(
  "median (50%)"        = quantile(df2$risk, 0.50, names = FALSE),
  "tertile (33/67%)"    = quantile(df2$risk, c(1/3, 2/3), names = FALSE),
  "quartile (25/75%)"   = quantile(df2$risk, c(0.25, 0.75), names = FALSE)
)

res_rows <- list()
logrank <- function(grp) {
  s <- survdiff(Surv(tte_adj, event) ~ grp, data = df2)
  1 - pchisq(s$chisq, length(s$n) - 1)
}
# (a) 连续
res_rows[[length(res_rows) + 1]] <- data.frame(
  scheme = "continuous risk score",
  n_groups = NA,
  HR = s_cont$coefficients[1, "exp(coef)"],
  lo = s_cont$conf.int[1, "lower .95"],
  hi = s_cont$conf.int[1, "upper .95"],
  P  = s_cont$coefficients[1, "Pr(>|z|)"]
)
# (b-1) 中位二分（现稿口径）
g2 <- factor(ifelse(df2$risk > cuts[[1]], "high", "low"), levels = c("low", "high"))
m2 <- coxph(Surv(tte_adj, event) ~ g2, data = df2, ties = "efron")
s2 <- summary(m2)
res_rows[[length(res_rows) + 1]] <- data.frame(
  scheme = "binary at median", n_groups = 2,
  HR = s2$coefficients[1, "exp(coef)"], lo = s2$conf.int[1, "lower .95"],
  hi = s2$conf.int[1, "upper .95"], P = s2$coefficients[1, "Pr(>|z|)"]
)
# (b-2) 三分位（低 vs 高，剔除中位组）
g3 <- cut(df2$risk, c(-Inf, cuts[[2]][1], cuts[[2]][2], Inf), labels = c("low", "mid", "high"))
df2$g3 <- factor(ifelse(g3 == "mid", NA, as.character(g3)), levels = c("low", "high"))
sub3 <- df2[!is.na(df2$g3), ]
m3 <- coxph(Surv(tte_adj, event) ~ g3, data = sub3, ties = "efron")
s3 <- summary(m3)
res_rows[[length(res_rows) + 1]] <- data.frame(
  scheme = "tertile extremes (low vs high, mid dropped)", n_groups = 2,
  HR = s3$coefficients[1, "exp(coef)"], lo = s3$conf.int[1, "lower .95"],
  hi = s3$conf.int[1, "upper .95"], P = s3$coefficients[1, "Pr(>|z|)"]
)
# (b-3) 四分位（低 vs 高）
g4 <- cut(df2$risk, c(-Inf, cuts[[3]][1], cuts[[3]][2], Inf), labels = c("low", "mid", "high"))
df2$g4 <- factor(ifelse(g4 == "mid", NA, as.character(g4)), levels = c("low", "high"))
sub4 <- df2[!is.na(df2$g4), ]
m4 <- coxph(Surv(tte_adj, event) ~ g4, data = sub4, ties = "efron")
s4 <- summary(m4)
res_rows[[length(res_rows) + 1]] <- data.frame(
  scheme = "quartile extremes (low vs high, mid dropped)", n_groups = 2,
  HR = s4$coefficients[1, "exp(coef)"], lo = s4$conf.int[1, "lower .95"],
  hi = s4$conf.int[1, "upper .95"], P = s4$coefficients[1, "Pr(>|z|)"]
)
# (b-4) 三分位趋势检验（3 组 log-rank）
gp3 <- factor(as.character(g3), levels = c("low", "mid", "high"))
res_rows[[length(res_rows) + 1]] <- data.frame(
  scheme = "tertile 3-group log-rank (trend)", n_groups = 3,
  HR = NA, lo = NA, hi = NA,
  P = logrank(gp3)
)

res <- do.call(rbind, res_rows)
res$HR <- round(res$HR, 3); res$lo <- round(res$lo, 3)
res$hi <- round(res$hi, 3); res$P  <- signif(res$P, 3)
cat("\n--- 切点敏感性汇总 ---\n"); print(res, row.names = FALSE)
write.csv(res, file.path(OUT, "tables", "step109_cutpoint_sensitivity.csv"), row.names = FALSE)

log_line("")
log_line("=== step109 完成 ===")
