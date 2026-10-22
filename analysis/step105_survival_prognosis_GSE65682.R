#!/usr/bin/env Rscript
# ============================================================
# step105_survival_prognosis_GSE65682.R
# 生存分析与预后模型构建与评价（GSE65682, MARS 脓毒性休克队列）
#
# 设计：
#   1. 解析 GSE65682 series matrix（表达 + 临床特征）
#   2. 解析 family.soft 平台表获得探针→基因符号（GPL13667, Affy HG 2.0 ST）
#   3. 纳入 479 例有 28 天生存终点的脓毒性休克患者（剔除配对随访样本，仅保留每患者首份样本）
#   4. 6 个 hub 基因表达提取（多探针取最高平均表达）
#   5. 单因素 Cox（各基因，z 标化）→ 多因素 Cox（6 基因）→ 风险评分
#   6. KM 曲线（中位风险评分分组，log-rank）
#   7. 时间依赖 ROC（7/14/28 天，timeROC）
#   8. C 指数（bootstrap 95% CI）
#   9. 校准曲线（28 天，rms::calibrate）+ nomogram（rms）
#  10. 临床协变量调整模型（年龄 + 性别 + 风险评分）森林图
#
# 输出：13_survival_prognosis/{figures,tables,models,logs}
# 图形：PDF 矢量 + PNG@300dpi（项目惯例）
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(timeROC)
  library(rms)
  library(ggplot2)
})

BASE  <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
SM    <- file.path(BASE, "00_rawdata/validation/GSE65682/GSE65682_series_matrix.txt")
SOFT  <- file.path(BASE, "00_rawdata/validation/GSE65682/GSE65682_family.soft.gz")
OUT   <- file.path(BASE, "13_survival_prognosis")
for (d in c("figures", "tables", "models", "logs")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
LOGF  <- file.path(OUT, "logs", "step105_run.log")
log_line <- function(...) { msg <- paste0(...); cat(msg, "\n"); cat(msg, "\n", file = LOGF, append = TRUE) }

HUB <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")

# NPG 配色（与项目一致）
COL_HI <- "#E64B35FF"   # 高风险
COL_LO <- "#3C5488FF"   # 低风险
pal_npg_roc <- function() c("#E64B35FF", "#4DBBD5FF", "#00A087FF")

# ------------------------------------------------------------
# 1. 解析 series matrix：临床特征
# ------------------------------------------------------------
log_line("=== 解析 series matrix ===")
lines <- readLines(SM, warn = FALSE)

get_field <- function(tag, field) {
  ln <- lines[grepl(paste0("^", tag), lines)]
  if (nzchar(field)) hit <- ln[grepl(paste0("^", tag, ".*", field, " *:"), ln)] else hit <- ln
  if (length(hit) == 0) return(NULL)
  vals <- strsplit(hit, "\t")[[1]][-1]
  gsub('^"|"$', "", vals)
}

gsm   <- get_field("!Sample_geo_accession", "")
title <- get_field("!Sample_title", "")
event <- suppressWarnings(as.integer(gsub(".*: ", "", get_field("!Sample_characteristics_ch1", "mortality_event_28days"))))
tte   <- suppressWarnings(as.numeric(gsub(".*: ", "", get_field("!Sample_characteristics_ch1", "time_to_event_28days"))))
age   <- suppressWarnings(as.numeric(gsub(".*: ", "", get_field("!Sample_characteristics_ch1", "age"))))
sex   <- gsub(".*: ", "", get_field("!Sample_characteristics_ch1", "gender"))
endo  <- gsub(".*: ", "", get_field("!Sample_characteristics_ch1", "endotype_class"))

clin <- data.frame(gsm = gsm, title = title, event = event, tte = tte,
                   age = age, sex = sex, endotype = endo, stringsAsFactors = FALSE)
log_line(sprintf("总样本 n=%d；有 28 天终点 n=%d（死亡 %d，存活 %d）",
                 nrow(clin), sum(!is.na(clin$event)), sum(clin$event == 1, na.rm = TRUE),
                 sum(clin$event == 0, na.rm = TRUE)))

# 仅保留有生存终点的脓毒性休克患者
clin <- clin[!is.na(clin$event) & !is.na(clin$tte), ]
# 每患者仅保留首份样本（title 中方括号内为患者编号/条码，取唯一患者首个 GSM 顺序）
clin$patient <- sub(".*\\[", "", clin$title)
clin$patient <- sub("\\]$", "", clin$patient)
# 患者标识：末段数字串（如 3270）
clin$pid <- regmatches(clin$patient, regexpr("[0-9]+_[0-9]+_[0-9]+_[A-H][0-9]+_[0-9a-z]+$", clin$patient))
clin$pid <- ifelse(is.na(clin$pid), clin$patient, clin$pid)
# 提取纯编号（去板位）：B02_3270 → 3270
clin$pnum <- sub(".*[_]([0-9a-z]+)$", "\\1", clin$patient)
dup <- sum(duplicated(clin$pnum))
log_line(sprintf("疑似同患者多份样本数（按编号去重前 n=%d，重复 %d）", nrow(clin), dup))
if (dup > 0) clin <- clin[!duplicated(clin$pnum), ]
log_line(sprintf("纳入分析样本 n=%d", nrow(clin)))

# tte=0 的事件（11 例当日死亡）调整为 0.5 天避免 Cox 报错
clin$tte_adj <- ifelse(clin$tte == 0 & clin$event == 1, 0.5, clin$tte)

# ------------------------------------------------------------
# 2. 探针注释（GPL13667 = Affymetrix HG-U219；由 family.soft 前部部分解压提取）
# ------------------------------------------------------------
log_line("=== 读取平台注释（GPL13667）===")
ANN <- file.path(BASE, "00_rawdata/validation/GSE65682/GPL13667_probe_symbol.tsv")
probe_map <- read.delim(ANN, stringsAsFactors = FALSE)
probe_map <- probe_map[probe_map$symbol != "" & !is.na(probe_map$symbol), ]
log_line(sprintf("平台注释探针 n=%d（含基因符号）", nrow(probe_map)))

# ------------------------------------------------------------
# 3. 表达矩阵提取（hub 基因）
# ------------------------------------------------------------
log_line("=== 提取 hub 基因表达 ===")
expr_start <- grep("^!series_matrix_table_begin", lines)
expr_end   <- grep("^!series_matrix_table_end", lines)
hdr_line   <- lines[expr_start + 1]
expr_gsm   <- gsub('^"|"$', "", strsplit(hdr_line, "\t")[[1]][-1])
gene_probes <- probe_map[probe_map$symbol %in% HUB, ]

# 逐行扫描表达表，取目标探针（避免整表载入内存）
want <- unique(gene_probes$probe)
rows <- list(); row_ids <- character(0)
con <- file(SM, "r")
in_tab <- FALSE; done <- FALSE
while (!done && length(l <- readLines(con, n = 2000, warn = FALSE)) > 0) {
  for (i in seq_along(l)) {
    if (grepl("^!series_matrix_table_begin", l[i])) { in_tab <- TRUE; next }
    if (grepl("^!series_matrix_table_end", l[i])) { in_tab <- FALSE; done <- TRUE; break }
    if (in_tab) {
      parts <- strsplit(l[i], "\t", fixed = TRUE)[[1]]
      p <- gsub('^"|"$', "", parts[1])
      if (p %in% want) {
        rows[[length(rows) + 1]] <- as.numeric(gsub('^"|"$', "", parts[-1]))
        row_ids <- c(row_ids, p)
      }
    }
  }
  if (length(row_ids) == length(want)) done <- TRUE
}
close(con)
emat <- do.call(rbind, rows)
rownames(emat) <- row_ids
colnames(emat) <- expr_gsm
gene_probes <- probe_map[probe_map$probe %in% row_ids & probe_map$symbol %in% HUB, ]
gene_probes <- probe_map[probe_map$probe %in% rownames(emat) & probe_map$symbol %in% HUB, ]
log_line(sprintf("hub 基因相关探针 n=%d", nrow(gene_probes)))

# 多探针 → 基因（取全队列平均表达最高的探针）
samp <- clin$gsm
gexpr <- sapply(HUB, function(g) {
  pr <- gene_probes$probe[gene_probes$symbol == g]
  pr <- intersect(pr, rownames(emat))
  if (length(pr) == 0) return(NA_real_)
  if (length(pr) == 1) return(emat[pr, samp])
  means <- rowMeans(emat[pr, samp, drop = FALSE])
  emat[pr[which.max(means)], samp]
})
gexpr <- t(gexpr)
colnames(gexpr) <- samp
log_line(paste0("各基因探针数：", paste(sapply(HUB, function(g) sum(gene_probes$symbol == g)), collapse = "/")))

df <- data.frame(clin, t(gexpr), stringsAsFactors = FALSE)
df$sex_f <- factor(df$sex, levels = c("male", "female"))
stopifnot(!any(is.na(df[, HUB])))

# ------------------------------------------------------------
# 4. 单因素 Cox（z 标化表达）
# ------------------------------------------------------------
log_line("=== 单因素 Cox（28 天死亡）===")
uni <- do.call(rbind, lapply(HUB, function(g) {
  z <- scale(df[[g]])
  fit <- coxph(Surv(tte_adj, event) ~ z, data = df)
  s <- summary(fit)
  data.frame(gene = g, beta = coef(fit), HR = s$conf.int[1],
             lo = s$conf.int[3], hi = s$conf.int[4],
             z = s$coefficients[4], P = s$coefficients[5], stringsAsFactors = FALSE)
}))
write.csv(uni, file.path(OUT, "tables", "univariate_cox_hub_genes.csv"), row.names = FALSE)
log_line(capture.output(print(uni[, c("gene", "HR", "lo", "hi", "P")], digits = 3)))

# ------------------------------------------------------------
# 5. 多因素 Cox + 风险评分
# ------------------------------------------------------------
log_line("=== 多因素 Cox（6 基因）===")
for (g in HUB) df[[paste0(g, "_z")]] <- as.numeric(scale(df[[g]]))
fml <- as.formula(paste("Surv(tte_adj, event) ~", paste(paste0(HUB, "_z"), collapse = " + ")))
mfit <- coxph(fml, data = df)
ms <- summary(mfit)
print(mfit)
log_line(capture.output(print(mfit)))

coef_v <- coef(mfit)
df$riskScore <- as.numeric(as.matrix(df[, paste0(HUB, "_z")]) %*% coef_v)
df$riskGroup <- factor(ifelse(df$riskScore >= median(df$riskScore), "High", "Low"),
                       levels = c("Low", "High"))
write.csv(df[, c("gsm", "title", "patient", "event", "tte", "age", "sex", "endotype", HUB, "riskScore", "riskGroup")],
          file.path(OUT, "tables", "GSE65682_risk_scores.csv"), row.names = FALSE)
saveRDS(mfit, file.path(OUT, "models", "cox6_model.rds"))

# 基因间 Spearman 相关（共线性检查）
cmat <- cor(df[, HUB], method = "spearman")
write.csv(round(cmat, 3), file.path(OUT, "tables", "hub_gene_spearman_GSE65682.csv"))

# ------------------------------------------------------------
# 6. KM 曲线
# ------------------------------------------------------------
log_line("=== KM 曲线（中位风险评分分组）===")
km <- survfit(Surv(tte_adj, event) ~ riskGroup, data = df)
lr <- survdiff(Surv(tte_adj, event) ~ riskGroup, data = df)
p_lr <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)
log_line(sprintf("log-rank chi2=%.2f, P=%.3e", lr$chisq, p_lr))

p_km <- ggsurvplot(km, data = df, risk.table = TRUE, pval = TRUE,
                   conf.int = TRUE, palette = c(COL_LO, COL_HI),
                   legend.labs = c("Low risk", "High risk"),
                   xlab = "Days since ICU admission", ylab = "Survival probability",
                   risk.table.height = 0.28, font.family = "Helvetica",
                   ggtheme = theme_classic(base_size = 8))
pdf(file.path(OUT, "figures", "KM_riskScore_GSE65682.pdf"), width = 4.5, height = 4.5, family = "Helvetica")
print(p_km); dev.off()
png(file.path(OUT, "figures", "KM_riskScore_GSE65682.png"), width = 4.5, height = 4.5, units = "in", res = 300, family = "Helvetica")
print(p_km); dev.off()

# ------------------------------------------------------------
# 7. 时间依赖 ROC
# ------------------------------------------------------------
log_line("=== 时间依赖 ROC ===")
# 评估时点：7、14 天及 28 天（用 27.99 规避行政删失边界处 IPCW 权重退化导致的 NA）
TT <- c(7, 14, 27.99); TL <- c(7, 14, 28)
tr <- timeROC(T = df$tte_adj, delta = df$event, marker = df$riskScore,
              cause = 1, times = TT, iid = TRUE)
aucs <- sapply(seq_along(TT), function(i) tr$AUC[i])
ci_tab <- confint(tr, 27.99)$CI_AUC / 100   # 3x2，行为时点
roc_tab <- data.frame(time = TL, AUC = round(aucs, 3),
                      CI = sprintf("%.3f-%.3f", ci_tab[, 1], ci_tab[, 2]))
write.csv(roc_tab, file.path(OUT, "tables", "timeROC_riskScore.csv"), row.names = FALSE)
log_line(capture.output(print(roc_tab)))

roc_df <- do.call(rbind, lapply(seq_along(TT), function(i) {
  data.frame(FPR = 1 - tr$FP[, i], TPR = tr$TP[, i], time = TL[i])
}))
p_roc <- ggplot(roc_df, aes(FPR, TPR, color = factor(time))) +
  geom_line(linewidth = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
  scale_color_manual(values = pal_npg_roc(), name = NULL,
                     labels = sprintf("%d d (AUC=%.3f)", roc_tab$time, roc_tab$AUC)) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  coord_equal() +
  theme_classic(base_size = 8) +
  theme(legend.position = c(0.68, 0.18), legend.background = element_blank())
pdf(file.path(OUT, "figures", "timeROC_riskScore.pdf"), width = 3.2, height = 3.2, family = "Helvetica")
print(p_roc); dev.off()
png(file.path(OUT, "figures", "timeROC_riskScore.png"), width = 3.2, height = 3.2, units = "in", res = 300, family = "Helvetica")
print(p_roc); dev.off()

# ------------------------------------------------------------
# 8. C 指数（bootstrap）
# ------------------------------------------------------------
log_line("=== C 指数 ===")
set.seed(2026)
cindx <- survival::survConcordance(Surv(tte_adj, event) ~ riskScore, df)$concordance
nboot <- 500
bc <- numeric(nboot)
for (b in seq_len(nboot)) {
  idx <- sample(seq_len(nrow(df)), replace = TRUE)
  d_b <- df[idx, ]
  if (length(unique(d_b$riskGroup)) < 2 || sum(d_b$event) < 5) next
  bc[b] <- survival::survConcordance(Surv(tte_adj, event) ~ riskScore, d_b)$concordance
}
bc <- bc[bc != 0]
c_ci <- quantile(bc, c(0.025, 0.975))
log_line(sprintf("C-index = %.3f (bootstrap 95%% CI %.3f-%.3f, %d resamples)", cindx, c_ci[1], c_ci[2], length(bc)))

# ------------------------------------------------------------
# 9. 校准曲线 + nomogram（rms）
# ------------------------------------------------------------
log_line("=== 校准与 nomogram ===")
dd <- datadist(df); options(datadist = "dd")
cph_fit <- cph(Surv(tte_adj, event) ~ CDK5RAP2_z + CD177_z + GRB10_z + UPP1_z + EXOSC4_z + VNN1_z,
               data = df, x = TRUE, y = TRUE, surv = TRUE, time.inc = 28)
cal <- calibrate(cph_fit, u = 28, cmethod = "hare", m = 80, B = 400)
pdf(file.path(OUT, "figures", "calibration_28d.pdf"), width = 3.5, height = 3.5, family = "Helvetica")
plot(cal, xlab = "Predicted 28-day survival probability", ylab = "Observed",
     subtitles = FALSE, col = COL_HI, lwd = 2); dev.off()
png(file.path(OUT, "figures", "calibration_28d.png"), width = 3.5, height = 3.5, units = "in", res = 300, family = "Helvetica")
plot(cal, xlab = "Predicted 28-day survival probability", ylab = "Observed",
     subtitles = FALSE, col = COL_HI, lwd = 2); dev.off()

nom <- nomogram(cph_fit, fun = function(x) 1 - exp(-exp(x)),
                funlabel = "28-day mortality probability",
                lp = FALSE, fun.at = c(0.1, 0.2, 0.4, 0.6, 0.8))
pdf(file.path(OUT, "figures", "nomogram_28d.pdf"), width = 7, height = 4.5, family = "Helvetica")
plot(nom); dev.off()
png(file.path(OUT, "figures", "nomogram_28d.png"), width = 7, height = 4.5, units = "in", res = 300, family = "Helvetica")
plot(nom); dev.off()

# ------------------------------------------------------------
# 10. 协变量调整森林图（年龄 + 性别 + 风险评分）
# ------------------------------------------------------------
log_line("=== 协变量调整模型 ===")
adj <- coxph(Surv(tte_adj, event) ~ riskScore + age + sex_f, data = df)
as <- summary(adj)
forest_df <- data.frame(
  var = c("Risk score (per SD)", "Age (per year)", "Sex: female vs male"),
  HR = as$conf.int[, 1], lo = as$conf.int[, 3], hi = as$conf.int[, 4],
  P = as$coefficients[, 5])
write.csv(forest_df, file.path(OUT, "tables", "adjusted_cox_forest.csv"), row.names = FALSE)
log_line(capture.output(print(forest_df, digits = 3)))

forest_df$var <- factor(forest_df$var, levels = rev(forest_df$var))
p_for <- ggplot(forest_df, aes(HR, var)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
  geom_point(size = 2, shape = 16, color = COL_HI) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, color = COL_HI, linewidth = 0.5) +
  scale_x_log10() +
  labs(x = "Hazard ratio (log scale, 95% CI)", y = NULL) +
  theme_classic(base_size = 8)
pdf(file.path(OUT, "figures", "forest_adjusted_cox.pdf"), width = 3.8, height = 2.2, family = "Helvetica")
print(p_for); dev.off()
png(file.path(OUT, "figures", "forest_adjusted_cox.png"), width = 3.8, height = 2.2, units = "in", res = 300, family = "Helvetica")
print(p_for); dev.off()

# ------------------------------------------------------------
# 汇总输出
# ------------------------------------------------------------
log_line("=== 汇总 ===")
summary_tab <- data.frame(
  item = c("cohort", "n_patients", "deaths_28d", "logrank_P",
           "AUC_7d", "AUC_14d", "AUC_28d", "C_index", "C_index_CI"),
  value = c("GSE65682 septic shock (MARS)", nrow(df), sum(df$event),
            signif(p_lr, 3), round(aucs[1], 3), round(aucs[2], 3), round(aucs[3], 3),
            round(cindx, 3), sprintf("%.3f-%.3f", c_ci[1], c_ci[2])))
write.csv(summary_tab, file.path(OUT, "tables", "prognosis_summary.csv"), row.names = FALSE)
log_line(capture.output(print(summary_tab)))
log_line("完成。输出目录：", OUT)
