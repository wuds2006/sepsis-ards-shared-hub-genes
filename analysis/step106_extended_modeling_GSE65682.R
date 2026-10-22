#!/usr/bin/env Rscript
# ============================================================
# step106_extended_modeling_GSE65682.R
# GSE65682（MARS 脓毒性休克队列，n=479，28 天死亡终点）扩展预后建模
#
# 在 step105（6 hub 基因 Cox 模型）基础上的扩展：
#   1. 179 个共享候选基因表达提取（GPL13667 探针 → 基因符号）
#   2. LASSO-Cox（glmnet，10 折 CV，alpha=1）在 179 基因上建模与变量筛选
#      - 同时报告 elastic net（alpha=0.5）与 ridge（alpha=0）的 CV 偏似然偏差
#   3. 模型谱系比较：
#        M1 Clinical（年龄 + 性别 + MARS endotype）
#        M2 Hub6（6 个 hub 基因 z 标化）
#        M3 Clinical + Hub6
#        M4 LASSO-179 signature
#        M5 Clinical + LASSO-179 signature
#   4. 判别与校准评价：
#        - Harrell C（apparent）
#        - bootstrap（B=200）optimism 校正 C（LASSO 在 bootstrap 内嵌套重跑 CV）
#        - 时间依赖 AUC（timeROC，7/14/28 天）
#        - riskRegression::Score（AUC / Brier 差异检验，vs M1）
#        - NRI / IDI（28 天，Pencina 法）
#        - 决策曲线分析（28 天净获益）
#        - 分组校准（5 组，KM 观察 vs 平均预测）
#   5. 非线性对照：随机生存森林（ranger，179 基因 + 临床变量），OOB C 指数与 28 天 AUC
#
# 输出：13_survival_prognosis/{figures,tables,models,logs}
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(rms)
  library(timeROC)
  library(riskRegression)
  library(prodlim)
  library(ranger)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(ragg)
  library(scales)
})

set.seed(20261015)

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
PKG  <- file.path(BASE, "AJRCCM_submission/09_review_deAI_2026-09-25")
SM   <- file.path(BASE, "00_rawdata/validation/GSE65682/GSE65682_series_matrix.txt")
ANN  <- file.path(BASE, "00_rawdata/validation/GSE65682/GPL13667_probe_symbol.tsv")
SHARED <- file.path(BASE, "03_intersection/Shared_CGs.csv")
OUT  <- file.path(BASE, "13_survival_prognosis")
for (d in c("figures", "tables", "models", "logs")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
LOGF <- file.path(OUT, "logs", "step106_run.log")
log_line <- function(...) { msg <- paste0(...); cat(msg, "\n"); cat(msg, "\n", file = LOGF, append = TRUE) }
log_line("=== step106 扩展建模启动 ===")

HUB <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")
T0  <- 27.99          # 28 天时点（规避行政删失边界）
T0L <- 28

# NPG 配色（与项目一致）
pal <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF")

# ------------------------------------------------------------
# 1. 临床数据解析（与 step105 一致）
# ------------------------------------------------------------
log_line("=== 解析 series matrix 临床字段 ===")
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
clin <- clin[!is.na(clin$event) & !is.na(clin$tte), ]
clin$patient <- sub(".*\\[", "", clin$title); clin$patient <- sub("\\]$", "", clin$patient)
clin$pnum <- sub(".*[_]([0-9a-z]+)$", "\\1", clin$patient)
dup <- sum(duplicated(clin$pnum))
if (dup > 0) clin <- clin[!duplicated(clin$pnum), ]
clin$tte_adj <- ifelse(clin$tte == 0 & clin$event == 1, 0.5, clin$tte)
log_line(sprintf("纳入 n=%d（事件 %d，%.1f%%）", nrow(clin),
                 sum(clin$event == 1), 100 * mean(clin$event)))
log_line(paste0("endotype 分布：", paste(names(table(clin$endotype)), table(clin$endotype), sep = "=", collapse = ", ")))
log_line(sprintf("年龄：中位 %.0f（IQR %.0f-%.0f）；缺失 %d",
                 median(clin$age, na.rm = TRUE),
                 quantile(clin$age, 0.25, na.rm = TRUE), quantile(clin$age, 0.75, na.rm = TRUE),
                 sum(is.na(clin$age))))

# ------------------------------------------------------------
# 2. 179 个共享候选基因表达提取
# ------------------------------------------------------------
log_line("=== 读取平台注释与 179 基因列表 ===")
probe_map <- read.delim(ANN, stringsAsFactors = FALSE)
probe_map <- probe_map[!is.na(probe_map$symbol) & probe_map$symbol != "", ]
shared <- read.csv(SHARED, stringsAsFactors = FALSE)$gene
shared <- intersect(shared, unique(probe_map$symbol))
log_line(sprintf("共享候选基因 n=%d（在 GPL13667 上可映射 n=%d）", length(read.csv(SHARED)$gene), length(shared)))

want <- unique(probe_map$probe[probe_map$symbol %in% shared])
log_line(sprintf("目标探针 n=%d", length(want)))

log_line("=== 扫描表达矩阵（479 例样本列）===")
t0 <- Sys.time()
con <- file(SM, "r")
in_tab <- FALSE; done <- FALSE
rows <- list(); ids <- character(0)
while (!done && length(l <- readLines(con, n = 2000, warn = FALSE)) > 0) {
  for (i in seq_along(l)) {
    if (!in_tab) { if (grepl("^!series_matrix_table_begin", l[i])) in_tab <- TRUE; next }
    if (grepl("^!series_matrix_table_end", l[i])) { done <- TRUE; break }
    tp <- regexpr("\t", l[i], fixed = TRUE)
    if (tp < 0) next
    p <- substr(l[i], 2, tp - 2)
    if (p %in% want) {
      parts <- strsplit(l[i], "\t", fixed = TRUE)[[1]]
      rows[[length(rows) + 1]] <- as.numeric(gsub('^"|"$', "", parts[-1]))
      ids <- c(ids, p)
    }
  }
}
close(con)
emat <- do.call(rbind, rows); rownames(emat) <- ids
hdr <- gsub('^"|"$', "", strsplit(lines[grep("^!series_matrix_table_begin", lines) + 1], "\t")[[1]][-1])
colnames(emat) <- hdr
log_line(sprintf("提取耗时 %.1f s；矩阵 %d 探针 × %d 样本", as.numeric(difftime(Sys.time(), t0, units = "secs")), nrow(emat), ncol(emat)))

samp <- clin$gsm
emat <- emat[, samp, drop = FALSE]
pm <- probe_map[probe_map$probe %in% rownames(emat) & probe_map$symbol %in% shared, ]

# 多探针 → 基因：取全队列平均表达最高者（与 step105 一致）
gexpr <- do.call(cbind, lapply(shared, function(g) {
  pr <- intersect(pm$probe[pm$symbol == g], rownames(emat))
  if (length(pr) == 0) return(rep(NA_real_, ncol(emat)))
  if (length(pr) == 1) return(emat[pr, ])
  means <- rowMeans(emat[pr, , drop = FALSE])
  emat[pr[which.max(means)], ]
}))
colnames(gexpr) <- shared
miss <- shared[colSums(is.na(gexpr)) > 0]
if (length(miss) > 0) log_line(paste0("无法映射的基因（剔除）：", paste(miss, collapse = ", ")))
gexpr <- gexpr[, colSums(is.na(gexpr)) == 0, drop = FALSE]
log_line(sprintf("最终基因矩阵：%d 基因 × %d 样本", ncol(gexpr), nrow(gexpr)))

df <- data.frame(clin, gexpr, stringsAsFactors = FALSE)
df$sex_f <- factor(df$sex, levels = c("male", "female"))
df$endotype <- factor(df$endotype)
df$age_imp <- ifelse(is.na(df$age), median(df$age, na.rm = TRUE), df$age)  # 年龄缺失用中位数填补
log_line(sprintf("年龄缺失填补 n=%d", sum(is.na(df$age))))
saveRDS(list(df = df, gexpr = gexpr), file.path(OUT, "models", "step106_data.rds"))

Y <- Surv(df$tte_adj, df$event)
X179 <- as.matrix(gexpr)
XHUB <- as.matrix(scale(df[, HUB]))
stopifnot(!any(is.na(X179)), !any(is.na(XHUB)))

# ------------------------------------------------------------
# 3. 工具函数
# ------------------------------------------------------------
# rcorr.cens 的 C Index 定义为 x 与生存时间的一致性；此处 x 为风险评分（越高越早死亡），
# 故预后 C 指数 = 1 - C Index = (1 - Dxy)/2
harrellC <- function(lp, y = Y) as.numeric(1 - (rcorr.cens(lp, y)["Dxy"] + 1) / 2)
# glmnet 的 Cox 模型无截距项；coef() 返回 171 行（=基因数），不再含 (Intercept)
get_beta <- function(bmat) {
  nm <- rownames(bmat)
  if (length(nm) && "(Intercept)" %in% nm) return(as.numeric(bmat[-which(nm == "(Intercept)"), 1]))
  as.numeric(bmat[, 1])
}
auc_at <- function(marker, times = T0, y_time = df$tte_adj, y_event = df$event) {
  tr <- timeROC(T = y_time, delta = y_event, marker = marker, cause = 1, times = times, iid = TRUE)
  # timeROC 在单时点时会在结果首行插入 t=0（AUC=NA），故取最后一行
  A <- as.matrix(tr$AUC)
  a <- as.numeric(A[nrow(A), 1])
  ci <- tryCatch({
    m <- as.matrix(confint(tr, times)$CI_AUC)
    if (ncol(m) == 2) as.numeric(m[nrow(m), ]) / 100 else NULL
  }, error = function(e) NULL)
  list(auc = a, ci = ci)
}
pred28 <- function(fit, newdata = df) {
  sf <- survfit(fit, newdata = newdata)
  s <- summary(sf, times = T0)$surv
  as.numeric(1 - s)
}
nri_idi <- function(p_old, p_new, cutoff = NULL) {
  # 28 天二分结局：事件（28 天内死亡）vs 对照（随访至 28 天仍存活）
  case <- df$event == 1
  ctrl <- df$event == 0 & df$tte_adj >= 27
  if (is.null(cutoff)) cutoff <- sum(case) / (sum(case) + sum(ctrl))
  # IDI
  idi <- (mean(p_new[case]) - mean(p_old[case])) - (mean(p_new[ctrl]) - mean(p_old[ctrl]))
  # 分类 NRI
  up <- p_new >= cutoff & p_old < cutoff
  dn <- p_new < cutoff & p_old >= cutoff
  nri <- (mean(up[case]) - mean(dn[case])) - (mean(up[ctrl]) - mean(dn[ctrl]))
  # 各分量
  list(NRI = nri, IDI = idi,
       p_up_case = mean(up[case]), p_dn_case = mean(dn[case]),
       p_up_ctrl = mean(up[ctrl]), p_dn_ctrl = mean(dn[ctrl]),
       n_case = sum(case), n_ctrl = sum(ctrl), cutoff = cutoff)
}
# 二分 NRI/IDI 的 z 检验（Pencina）
nri_p <- function(res) {
  se <- sqrt(res$p_up_case / res$n_case + res$p_dn_case / res$n_case +
             res$p_up_ctrl / res$n_ctrl + res$p_dn_ctrl / res$n_ctrl)
  z <- res$NRI / se; 2 * (1 - pnorm(abs(z)))
}
calib_group <- function(p, ngroup = 5) {
  g <- cut(p, breaks = quantile(p, probs = seq(0, 1, length.out = ngroup + 1), names = FALSE),
           include.lowest = TRUE, labels = FALSE)
  obs <- sapply(sort(unique(g)), function(k) {
    idx <- which(g == k)
    sf <- survfit(Surv(tte_adj, event) ~ 1, data = df[idx, ])
    1 - as.numeric(summary(sf, times = T0)$surv)
  })
  data.frame(group = sort(unique(g)), n = as.integer(table(g)),
             pred = as.numeric(tapply(p, g, mean)), obs = as.numeric(obs))
}
dca_curve <- function(p, thresh = seq(0.05, 0.60, by = 0.01)) {
  case <- df$event == 1
  ctrl <- df$event == 0 & df$tte_adj >= 27
  keep <- case | ctrl
  y <- as.integer(df$event[keep] == 1); pp <- p[keep]; n <- length(y)
  sapply(thresh, function(pt) {
    tr <- pp >= pt
    sum(y == 1 & tr) / n - sum(y == 0 & tr) / n * (pt / (1 - pt))
  })
}
save_pub <- function(plot_obj, name, w, h) {
  tryCatch({
    pdf(file.path(OUT, "figures", paste0(name, ".pdf")), width = w, height = h, family = "Helvetica")
    print(plot_obj); dev.off()
    agg_png(file.path(OUT, "figures", paste0(name, ".png")), width = w, height = h, units = "in", res = 300,
            background = "white")
    print(plot_obj); dev.off()
  }, error = function(e) {
    while (dev.cur() > 1) dev.off()
    log_line(sprintf("[图形跳过] %s：%s", name, conditionMessage(e)))
  })
  invisible(NULL)
}
theme_pub <- function(base_size = 8) theme_classic(base_size = base_size) +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_text(colour = "black"),
        legend.title = element_blank())

# ------------------------------------------------------------
# 4. 模型拟合
# ------------------------------------------------------------
log_line("=== 模型拟合 ===")
m_clin  <- coxph(Surv(tte_adj, event) ~ age_imp + sex_f + endotype, data = df)
m_hub   <- coxph(Y ~ XHUB)
m_clin_hub <- coxph(Surv(tte_adj, event) ~ age_imp + sex_f + endotype + CDK5RAP2 + CD177 + GRB10 + UPP1 + EXOSC4 + VNN1,
                    data = df)  # 未标化，仅用于 LRT；主分析用 z 标化

# z 标化版临床 + hub6（与 step105 一致）
for (g in HUB) df[[paste0(g, "_z")]] <- as.numeric(scale(df[[g]]))
m_hub6z <- coxph(as.formula(paste("Surv(tte_adj, event) ~", paste(paste0(HUB, "_z"), collapse = " + "))), data = df)
m_clin_hub6z <- coxph(as.formula(paste("Surv(tte_adj, event) ~ age_imp + sex_f + endotype +",
                                       paste(paste0(HUB, "_z"), collapse = " + "))), data = df)

# LASSO-Cox（179 基因）
log_line("=== LASSO-Cox（179 基因，10 折 CV）===")
t0 <- Sys.time()
cv_lasso <- cv.glmnet(X179, Y, family = "cox", alpha = 1, nfolds = 10, standardize = TRUE,
                      control = list(maxit = 1e5), cox.ties = "efron")
log_line(sprintf("LASSO CV 耗时 %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
lam_min <- cv_lasso$lambda.min; lam_1se <- cv_lasso$lambda.1se
fit_lasso <- glmnet(X179, Y, family = "cox", alpha = 1, standardize = TRUE,
                    control = list(maxit = 1e5), cox.ties = "efron")
b_1se <- as.matrix(coef(fit_lasso, s = lam_1se))
b_min <- as.matrix(coef(fit_lasso, s = lam_min))
sel_1se <- rownames(b_1se)[b_1se[, 1] != 0]
sel_min <- rownames(b_min)[b_min[, 1] != 0]
log_line(sprintf("lambda.min=%.4f（非零 %d）；lambda.1se=%.4f（非零 %d）", lam_min, length(sel_min), lam_1se, length(sel_1se)))
log_line(paste0("lambda.1se 选中基因：", paste(sel_1se, collapse = ", ")))
log_line(paste0("lambda.min 选中基因：", paste(sel_min, collapse = ", ")))
log_line(paste0("6 个 hub 基因是否被 LASSO 选中（1se）：",
                paste(sapply(HUB, function(g) paste0(g, "=", g %in% sel_1se)), collapse = ", ")))

# elastic net / ridge 对照（CV 偏似然偏差）
cv_en    <- cv.glmnet(X179, Y, family = "cox", alpha = 0.5, nfolds = 10, standardize = TRUE,
                      control = list(maxit = 1e5), cox.ties = "efron")
cv_ridge <- cv.glmnet(X179, Y, family = "cox", alpha = 0, nfolds = 10, standardize = TRUE,
                      control = list(maxit = 1e5), cox.ties = "efron")
fit_en    <- glmnet(X179, Y, family = "cox", alpha = 0.5, standardize = TRUE,
                    control = list(maxit = 1e5), cox.ties = "efron")
fit_ridge <- glmnet(X179, Y, family = "cox", alpha = 0, standardize = TRUE,
                    control = list(maxit = 1e5), cox.ties = "efron")
penalty_tab <- data.frame(
  penalty = c("LASSO (alpha=1)", "Elastic net (alpha=0.5)", "Ridge (alpha=0)"),
  lambda_min = c(lam_min, cv_en$lambda.min, cv_ridge$lambda.min),
  lambda_1se = c(lam_1se, cv_en$lambda.1se, cv_ridge$lambda.1se),
  n_nonzero_min = c(length(sel_min),
                    sum(get_beta(as.matrix(coef(fit_en, s = cv_en$lambda.min)) != 0)),
                    sum(get_beta(as.matrix(coef(fit_ridge, s = cv_ridge$lambda.min)) != 0))),
  cv_deviance_min = c(min(cv_lasso$cvm), min(cv_en$cvm), min(cv_ridge$cvm)),
  cv_deviance_se  = c(cv_lasso$cvsd[which.min(cv_lasso$cvm)],
                      cv_en$cvsd[which.min(cv_en$cvm)],
                      cv_ridge$cvsd[which.min(cv_ridge$cvm)]))
write.csv(penalty_tab, file.path(OUT, "tables", "ext_penalty_comparison.csv"), row.names = FALSE)
log_line(capture.output(print(penalty_tab, digits = 4)))

# LASSO signature 评分
# 说明：1-SE 规则（lambda.1se）在 179 个基因上未保留任何变量，说明该 panel 的预后信号很弱；
#       主签名因此采用 CV 最小偏似然偏差对应的 lambda.min，并在正文中如实报告 1-SE 结果。
beta <- if (length(sel_1se) == 0) get_beta(b_min) else get_beta(b_1se)
lambda_used <- if (length(sel_1se) == 0) lam_min else lam_1se
rule_used <- if (length(sel_1se) == 0) "lambda.min (1-SE rule retained no gene)" else "lambda.1se"
log_line(paste0("采用规则：", rule_used))
df$lasso_lp <- as.numeric(X179 %*% beta)
m_lasso <- coxph(Surv(tte_adj, event) ~ lasso_lp, data = df)
m_clin_lasso <- coxph(Surv(tte_adj, event) ~ age_imp + sex_f + endotype + lasso_lp, data = df)

b_used <- if (length(sel_1se) == 0) b_min else b_1se
lasso_coef_tab <- data.frame(gene = rownames(b_used),
                             beta_standardised = get_beta(b_used),
                             HR = exp(get_beta(b_used)),
                             stringsAsFactors = FALSE)
lasso_coef_tab <- lasso_coef_tab[abs(lasso_coef_tab$beta_standardised) > 0, ]
lasso_coef_tab <- lasso_coef_tab[order(-abs(lasso_coef_tab$beta_standardised)), ]
lasso_coef_tab$selected_hub <- lasso_coef_tab$gene %in% HUB
write.csv(lasso_coef_tab, file.path(OUT, "tables", "ext_lasso_coefficients.csv"), row.names = FALSE)

# 随机生存森林（ranger）
log_line("=== 随机生存森林（ranger）===")
rsf_dat <- data.frame(tte_adj = df$tte_adj, event = df$event, age_imp = df$age_imp,
                      sex_f = df$sex_f, endotype = df$endotype, gexpr)
rsf <- ranger(Surv(tte_adj, event) ~ ., data = rsf_dat, num.trees = 1000,
              mtry = floor(sqrt(ncol(rsf_dat) - 2)), importance = "impurity",
              min.node.size = 15, seed = 20261015, num.threads = 4)
rsf_C <- 1 - rsf$prediction.error
# 个体水平的 RSF 风险采用 10 折交叉验证的折外预测（in-bag 预测在此样本量下严重乐观）
set.seed(20261015)
fold <- sample(rep(1:10, length.out = nrow(rsf_dat)))
rsf_cv <- rep(NA_real_, nrow(rsf_dat))
for (k in 1:10) {
  tr_i <- which(fold != k); te_i <- which(fold == k)
  rf_k <- ranger(Surv(tte_adj, event) ~ ., data = rsf_dat[tr_i, ], num.trees = 500,
                 mtry = floor(sqrt(ncol(rsf_dat) - 2)), min.node.size = 15,
                 seed = 20261015 + k, num.threads = 4)
  pk <- predict(rf_k, data = rsf_dat[te_i, ])
  ck <- max(which(pk$unique.death.times <= T0))
  rsf_cv[te_i] <- 1 - as.numeric(pk$survival[, ck])
}
df$rsf_risk <- rsf_cv
log_line(sprintf("RSF 折外预测 n=%d（缺失 %d）；折外 C 指数=%.3f",
                 sum(!is.na(df$rsf_risk)), sum(is.na(df$rsf_risk)), harrellC(df$rsf_risk)))
log_line(sprintf("RSF OOB C 指数 = %.3f", rsf_C))
imp <- sort(rsf$variable.importance, decreasing = TRUE)
write.csv(data.frame(variable = names(imp), importance = as.numeric(imp))[1:20, ],
          file.path(OUT, "tables", "ext_rsf_importance.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 5. 判别性能比较
# ------------------------------------------------------------
log_line("=== 判别性能（apparent C / timeROC）===")
models <- list(`Clinical` = m_clin, `Hub6` = m_hub6z,
               `Clinical + Hub6` = m_clin_hub6z,
               `LASSO-179` = m_lasso, `Clinical + LASSO-179` = m_clin_lasso)
lp_list <- list(`Clinical` = predict(m_clin, type = "lp"),
                `Hub6` = predict(m_hub6z, type = "lp"),
                `Clinical + Hub6` = predict(m_clin_hub6z, type = "lp"),
                `LASSO-179` = predict(m_lasso, type = "lp"),
                `Clinical + LASSO-179` = predict(m_clin_lasso, type = "lp"),
                `RSF (OOB)` = df$rsf_risk)

auc1 <- function(marker, t) as.numeric(auc_at(marker, t)$auc)[1]
ci1  <- function(marker, t) {
  v <- auc_at(marker, t)$ci
  if (is.null(v) || length(v) < 2) NA_character_ else paste0(sprintf("%.3f", v[1]), "-", sprintf("%.3f", v[2]))
}
perf <- do.call(rbind, lapply(names(lp_list), function(k) {
  lp <- lp_list[[k]]
  data.frame(model = k, C_apparent = round(as.numeric(harrellC(lp))[1], 3),
             AUC_7d = round(auc1(lp, 7), 3), AUC_14d = round(auc1(lp, 14), 3), AUC_28d = round(auc1(lp, T0), 3),
             AUC28_CI = ci1(lp, T0), stringsAsFactors = FALSE)
}))
aic_tab <- do.call(rbind, lapply(names(models), function(k) {
  f <- models[[k]]
  data.frame(model = k, AIC = round(AIC(f), 1), logLik = round(as.numeric(logLik(f)), 1),
             n_par = length(coef(f)), stringsAsFactors = FALSE)
}))
perf <- merge(perf, aic_tab, by = "model", all.x = TRUE)
write.csv(perf, file.path(OUT, "tables", "ext_model_performance.csv"), row.names = FALSE)
log_line(capture.output(print(perf)))

# 似然比检验：加入基因特征是否改善拟合
lrt_tab <- do.call(rbind, list(
  data.frame(comparison = "Clinical + Hub6 vs Clinical",
             chisq = round(as.numeric(2 * (as.numeric(logLik(m_clin_hub6z)) - as.numeric(logLik(m_clin)))), 2),
             df = 6,
             P = signif(pchisq(2 * (as.numeric(logLik(m_clin_hub6z)) - as.numeric(logLik(m_clin))), 6, lower.tail = FALSE), 3)),
  data.frame(comparison = "Clinical + LASSO-179 vs Clinical",
             chisq = round(as.numeric(2 * (as.numeric(logLik(m_clin_lasso)) - as.numeric(logLik(m_clin)))), 2),
             df = 1,
             P = signif(pchisq(2 * (as.numeric(logLik(m_clin_lasso)) - as.numeric(logLik(m_clin))), 1, lower.tail = FALSE), 3))
))
write.csv(lrt_tab, file.path(OUT, "tables", "ext_lrt.csv"), row.names = FALSE)
log_line(capture.output(print(lrt_tab)))

# ------------------------------------------------------------
# 6. Bootstrap 内部验证（optimism 校正 C 指数）
# ------------------------------------------------------------
B <- 300
log_line(sprintf("=== Bootstrap optimism 校正（B=%d）===", B))
log_line("LASSO 惩罚系数固定在全队列 CV 选定值（bootstrap 内不重跑 λ 选择；这会轻度低估 optimism，见 Methods）")
n <- nrow(df)
boot_res <- list(
  `Clinical` = numeric(B), `Hub6` = numeric(B),
  `Clinical + Hub6` = numeric(B), `LASSO-179` = numeric(B),
  `Clinical + LASSO-179` = numeric(B))
opt <- boot_res; app <- numeric(5)
names(app) <- names(boot_res)
app["Clinical"] <- harrellC(lp_list$`Clinical`)
app["Hub6"] <- harrellC(lp_list$`Hub6`)
app["Clinical + Hub6"] <- harrellC(lp_list$`Clinical + Hub6`)
app["LASSO-179"] <- harrellC(lp_list$`LASSO-179`)
app["Clinical + LASSO-179"] <- harrellC(lp_list$`Clinical + LASSO-179`)

t0 <- Sys.time()
for (b in seq_len(B)) {
  idx <- sample(n, replace = TRUE)
  db <- df[idx, ]
  yb <- Surv(db$tte_adj, db$event)
  # 临床 / hub6 固定模型
  f1 <- try(coxph(Surv(tte_adj, event) ~ age_imp + sex_f + endotype, data = db), silent = TRUE)
  f2 <- try(coxph(as.formula(paste("Surv(tte_adj, event) ~", paste(paste0(HUB, "_z"), collapse = " + "))), data = db), silent = TRUE)
  f3 <- try(coxph(as.formula(paste("Surv(tte_adj, event) ~ age_imp + sex_f + endotype +",
                                   paste(paste0(HUB, "_z"), collapse = " + "))), data = db), silent = TRUE)
  if (!inherits(f1, "try-error")) {
    boot_res$`Clinical`[b] <- harrellC(predict(f1, newdata = db, type = "lp"), yb)
    opt$`Clinical`[b] <- boot_res$`Clinical`[b] - harrellC(predict(f1, newdata = df, type = "lp"), Y)
  }
  if (!inherits(f2, "try-error")) {
    boot_res$`Hub6`[b] <- harrellC(predict(f2, newdata = db, type = "lp"), yb)
    opt$`Hub6`[b] <- boot_res$`Hub6`[b] - harrellC(predict(f2, newdata = df, type = "lp"), Y)
  }
  if (!inherits(f3, "try-error")) {
    boot_res$`Clinical + Hub6`[b] <- harrellC(predict(f3, newdata = db, type = "lp"), yb)
    opt$`Clinical + Hub6`[b] <- boot_res$`Clinical + Hub6`[b] - harrellC(predict(f3, newdata = df, type = "lp"), Y)
  }
  # LASSO：惩罚系数固定，重拟合系数
  fb <- try(glmnet(X179[idx, ], yb, family = "cox", alpha = 1, standardize = TRUE,
                   lambda = lambda_used, cox.ties = "efron"), silent = TRUE)
  if (!inherits(fb, "try-error")) {
    bb <- get_beta(as.matrix(coef(fb)))
    lpb_boot <- as.numeric(X179[idx, ] %*% bb)
    lpb_orig <- as.numeric(X179 %*% bb)
    c_boot <- harrellC(lpb_boot, yb); c_orig <- harrellC(lpb_orig, Y)
    boot_res$`LASSO-179`[b] <- c_boot; opt$`LASSO-179`[b] <- c_boot - c_orig
    db$lpb <- lpb_boot
    f5 <- try(coxph(Surv(tte_adj, event) ~ age_imp + sex_f + endotype + lpb, data = db), silent = TRUE)
    if (!inherits(f5, "try-error")) {
      df_lp <- df; df_lp$lpb <- lpb_orig
      boot_res$`Clinical + LASSO-179`[b] <- harrellC(predict(f5, newdata = db, type = "lp"), yb)
      opt$`Clinical + LASSO-179`[b] <- boot_res$`Clinical + LASSO-179`[b] -
        harrellC(predict(f5, newdata = df_lp, type = "lp"), Y)
    }
  }
  if (b %% 25 == 0) log_line(sprintf("  bootstrap %d/%d（%.0f s）", b, B, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
boot_tab <- do.call(rbind, lapply(names(app), function(k) {
  o <- opt[[k]]
  data.frame(model = k, C_apparent = round(app[k], 3),
             optimism = round(mean(o, na.rm = TRUE), 3),
             optimism_lo = round(quantile(o, 0.025, na.rm = TRUE), 3),
             optimism_hi = round(quantile(o, 0.975, na.rm = TRUE), 3),
             C_corrected = round(app[k] - mean(o, na.rm = TRUE), 3),
             stringsAsFactors = FALSE)
}))
write.csv(boot_tab, file.path(OUT, "tables", "ext_bootstrap_optimism.csv"), row.names = FALSE)
log_line(capture.output(print(boot_tab)))

# ------------------------------------------------------------
# 7. riskRegression::Score（AUC / Brier 差异检验，vs Clinical）
# ------------------------------------------------------------
log_line("=== riskRegression Score（AUC / Brier，vs Clinical）===")
sc <- try(Score(list(`Clinical` = m_clin, `Hub6` = m_hub6z, `Clinical + Hub6` = m_clin_hub6z,
                     `LASSO-179` = m_lasso, `Clinical + LASSO-179` = m_clin_lasso),
                formula = Hist(tte_adj, event) ~ 1, data = df,
                times = c(7, 14, T0), metrics = c("auc", "brier"),
                contrasts = TRUE, B = 200, seed = 1), silent = TRUE)
score_rows <- NULL
if (!inherits(sc, "try-error")) {
  s <- summary(sc, times = T0)
  auc_df <- as.data.frame(s$AUC$score)
  brs_df <- as.data.frame(s$Brier$score)
  score_rows <- data.frame(model = rownames(auc_df),
                           AUC_28d = round(auc_df[, 1], 3),
                           Brier_28d = round(brs_df[, 1], 4))
  write.csv(score_rows, file.path(OUT, "tables", "ext_score_auc_brier.csv"), row.names = FALSE)
  log_line(capture.output(print(score_rows)))
  # 差异检验
  dif <- try(as.data.frame(s$AUC$contrasts), silent = TRUE)
  if (!is.null(dif) && nrow(dif) > 0) {
    write.csv(dif, file.path(OUT, "tables", "ext_score_auc_contrasts.csv"))
    log_line("AUC 差异检验（vs Clinical）：")
    log_line(capture.output(print(dif, digits = 3)))
  }
}

# ------------------------------------------------------------
# 8. NRI / IDI（28 天）
# ------------------------------------------------------------
log_line("=== NRI / IDI（28 天）===")
p_clin <- pred28(m_clin)
p_hub  <- pred28(m_hub6z)
p_las  <- pred28(m_lasso)
p_ch   <- pred28(m_clin_hub6z)
p_cl   <- pred28(m_clin_lasso)
p_rsf  <- df$rsf_risk
cut <- sum(df$event == 1) / sum(df$event == 1 | (df$event == 0 & df$tte_adj >= 27))
res_cn <- nri_idi(p_clin, p_ch); res_cl <- nri_idi(p_clin, p_cl)
nri_tab <- data.frame(
  comparison = c("Clinical + Hub6 vs Clinical", "Clinical + LASSO-179 vs Clinical"),
  NRI = round(c(res_cn$NRI, res_cl$NRI), 3),
  NRI_P = signif(c(nri_p(res_cn), nri_p(res_cl)), 3),
  IDI = round(c(res_cn$IDI, res_cl$IDI), 4),
  cutoff = round(cut, 3),
  n_case = c(res_cn$n_case, res_cl$n_case), n_ctrl = c(res_cn$n_ctrl, res_cl$n_ctrl))
write.csv(nri_tab, file.path(OUT, "tables", "ext_nri_idi.csv"), row.names = FALSE)
log_line(capture.output(print(nri_tab)))

# ------------------------------------------------------------
# 9. 校准与 DCA
# ------------------------------------------------------------
log_line("=== 分组校准与决策曲线 ===")
cal_list <- list(`Clinical` = p_clin, `Hub6` = p_hub, `LASSO-179` = p_las,
                 `Clinical + Hub6` = p_ch, `Clinical + LASSO-179` = p_cl, `RSF` = p_rsf)
cal_df <- do.call(rbind, lapply(names(cal_list), function(k) {
  d <- calib_group(cal_list[[k]]); cbind(model = k, d)
}))
write.csv(cal_df, file.path(OUT, "tables", "ext_calibration.csv"), row.names = FALSE)
log_line(capture.output(print(cal_df, digits = 3)))

thr <- seq(0.05, 0.60, by = 0.01)
dca_df <- do.call(rbind, lapply(c(names(cal_list), "Treat all", "Treat none"), function(k) {
  if (k == "Treat all") nb <- sapply(thr, function(pt) {
      case <- df$event == 1; ctrl <- df$event == 0 & df$tte_adj >= 27
      keep <- case | ctrl; y <- as.integer(df$event[keep] == 1); n <- length(y)
      mean(y) - (1 - mean(y)) * (pt / (1 - pt))
    })
  else if (k == "Treat none") nb <- rep(0, length(thr))
  else nb <- dca_curve(cal_list[[k]], thr)
  data.frame(model = k, threshold = thr, net_benefit = as.numeric(nb))
}))
write.csv(dca_df, file.path(OUT, "tables", "ext_dca.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 10. 图形输出
# ------------------------------------------------------------
log_line("=== 生成图形 ===")

# (A) LASSO CV 曲线
dfcv <- data.frame(log_lambda = log(cv_lasso$lambda), cvm = cv_lasso$cvm, cvsd = cv_lasso$cvsd,
                   nz = cv_lasso$nzero)
pA <- ggplot(dfcv, aes(log_lambda, cvm)) +
  geom_errorbar(aes(ymin = cvm - cvsd, ymax = cvm + cvsd), width = 0.02, colour = "grey70") +
  geom_point(size = 1.2, colour = pal[1]) +
  geom_vline(xintercept = log(c(lam_min, lam_1se)), linetype = c(2, 3), colour = c(pal[2], pal[3])) +
  annotate("text", x = log(lam_min), y = max(dfcv$cvm) * 0.97, label = "lambda.min",
           vjust = 1.1, hjust = 1.05, size = 2.3, family = "Helvetica") +
  annotate("text", x = log(lam_1se), y = max(dfcv$cvm) * 0.80, label = "lambda.1se",
           vjust = 1.1, hjust = 1.05, size = 2.3, family = "Helvetica") +
  labs(x = expression(log(lambda)), y = "Partial likelihood deviance (10-fold CV)") +
  theme_pub() + theme(axis.title.x = element_text(size = 8))
save_pub(pA, "ext_lasso_cv", 3.4, 2.6)

# (B) LASSO 选择基因的效应（标准化 HR）
if (nrow(lasso_coef_tab) > 0) {
  tb <- lasso_coef_tab[order(lasso_coef_tab$beta_standardised), ]
  tb$gene <- factor(tb$gene, levels = tb$gene)
  tb$col <- ifelse(tb$gene %in% HUB, "Hub gene", "Shared candidate gene")
  pB <- ggplot(tb, aes(gene, beta_standardised, fill = col)) +
    geom_col(colour = "black", linewidth = 0.15, width = 0.75) +
    coord_flip() + scale_fill_manual(values = c("Hub gene" = pal[1], "Shared candidate gene" = pal[2])) +
    labs(x = NULL, y = "LASSO coefficient (log2 scale)") +
    theme_pub() + theme(legend.position = "bottom", legend.key.size = unit(3, "pt"),
                        legend.text = element_text(size = 6))
} else {
  pB <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No gene selected at lambda.1se", size = 3) +
    theme_void()
}
save_pub(pB, "ext_lasso_signature", 3.4, 2.6)

# (C) 模型 28 天时间依赖 ROC 对比
roc_list <- list(`Clinical` = lp_list$`Clinical`, `Hub6` = lp_list$`Hub6`,
                 `Clinical + Hub6` = lp_list$`Clinical + Hub6`,
                 `LASSO-179` = lp_list$`LASSO-179`, `Clinical + LASSO-179` = lp_list$`Clinical + LASSO-179`,
                 `RSF (OOB)` = df$rsf_risk)
roc_df <- do.call(rbind, lapply(names(roc_list), function(k) {
  tr <- timeROC(T = df$tte_adj, delta = df$event, marker = roc_list[[k]], cause = 1, times = T0, iid = TRUE)
  j <- ncol(as.matrix(tr$FP))
  data.frame(model = k, FPR = 1 - as.numeric(tr$FP[, j]), TPR = as.numeric(tr$TP[, j]),
             AUC = round(as.numeric(as.matrix(tr$AUC)[nrow(as.matrix(tr$AUC)), 1]), 3))
}))
roc_df <- roc_df[is.finite(roc_df$FPR) & is.finite(roc_df$TPR), ]
roc_df$model <- factor(roc_df$model, levels = names(roc_list))
pC <- ggplot(roc_df, aes(FPR, TPR, colour = model)) + geom_line(linewidth = 0.55, na.rm = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  scale_colour_manual(values = setNames(c("#999999FF", pal[1], pal[4], pal[3], pal[2], pal[5]), names(roc_list))) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_pub() + theme(legend.position = "bottom", legend.text = element_text(size = 6),
                      legend.key.size = unit(3, "pt"))
save_pub(pC, "ext_model_roc_28d", 3.4, 3.0)

# (D) 校准（5 组）
cal_df$model <- factor(cal_df$model, levels = names(cal_list))
cal_df <- cal_df[is.finite(cal_df$pred) & is.finite(cal_df$obs), ]
pD <- ggplot(cal_df, aes(pred, obs, colour = model, shape = model)) +
  geom_point(size = 1.6) + geom_line(linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  scale_colour_manual(values = setNames(c("#999999FF", pal[1], pal[3], pal[4], pal[2], pal[5]), names(cal_list))) +
  coord_fixed(ratio = 1, xlim = c(0, 0.6), ylim = c(0, 0.6), expand = FALSE) +
  labs(x = "Predicted 28-day mortality", y = "Observed (Kaplan-Meier)") +
  theme_pub() + theme(legend.position = "bottom", legend.text = element_text(size = 6),
                      legend.key.size = unit(3, "pt"))
save_pub(pD, "ext_calibration_models", 3.4, 3.0)

# (E) 决策曲线
dca_df$model <- factor(dca_df$model, levels = c(names(cal_list), "Treat all", "Treat none"))
dca_df <- dca_df[is.finite(dca_df$net_benefit), ]
pE <- ggplot(dca_df, aes(threshold, net_benefit, colour = model)) + geom_line(linewidth = 0.55) +
  coord_cartesian(ylim = c(-0.02, 0.25)) +
  scale_colour_manual(values = setNames(c("#999999FF", pal[1], pal[3], pal[4], pal[2], pal[5], "grey40", "grey40"),
                                        c(names(cal_list), "Treat all", "Treat none"))) +
  labs(x = "Threshold probability", y = "Net benefit") +
  theme_pub() + theme(legend.position = "bottom", legend.text = element_text(size = 5.5),
                      legend.key.size = unit(3, "pt"))
save_pub(pE, "ext_dca_28d", 3.4, 3.0)

# (F1) 增量价值（ΔC 校正后 / NRI / IDI）
inc <- data.frame(
  metric = c("Delta C-index", "Delta C-index", "NRI", "NRI", "IDI (x100)", "IDI (x100)"),
  model = c("Clinical + Hub6", "Clinical + LASSO-179", "Clinical + Hub6", "Clinical + LASSO-179",
            "Clinical + Hub6", "Clinical + LASSO-179"),
  value = c(boot_tab$C_corrected[boot_tab$model == "Clinical + Hub6"] - boot_tab$C_corrected[boot_tab$model == "Clinical"],
            boot_tab$C_corrected[boot_tab$model == "Clinical + LASSO-179"] - boot_tab$C_corrected[boot_tab$model == "Clinical"],
            res_cn$NRI, res_cl$NRI, res_cn$IDI * 100, res_cl$IDI * 100))
pF1 <- ggplot(inc, aes(model, value, fill = metric)) + geom_col(colour = "black", linewidth = 0.15, width = 0.6) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = setNames(pal[c(1, 2, 3)], c("Delta C-index", "NRI", "IDI (x100)"))) +
  labs(x = NULL, y = "Incremental value") +
  theme_pub(7) + theme(legend.position = "none", strip.text = element_text(size = 6),
                       axis.text.x = element_text(angle = 30, hjust = 1, size = 6))
save_pub(pF1, "ext_incremental_value", 3.4, 2.4)

# (F2) RSF vs Cox：C 指数对比 + Top 变量重要性
cmp <- data.frame(model = c("Clinical", "Hub6", "LASSO-179", "RSF"),
                  C = c(boot_tab$C_corrected[boot_tab$model == "Clinical"],
                        boot_tab$C_corrected[boot_tab$model == "Hub6"],
                        boot_tab$C_corrected[boot_tab$model == "LASSO-179"], rsf_C))
cmp$model <- factor(cmp$model, levels = cmp$model)
pF2a <- ggplot(cmp, aes(model, C)) + geom_col(fill = pal[4], colour = "black", width = 0.6) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
  coord_cartesian(ylim = c(0.45, 0.75)) +
  labs(x = NULL, y = "C-index (RSF: OOB)") +
  theme_pub(7) + theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6))
imp_df <- data.frame(variable = names(imp)[1:12], importance = as.numeric(imp)[1:12])
imp_df$variable <- factor(imp_df$variable, levels = rev(imp_df$variable))
imp_df$hub <- imp_df$variable %in% HUB
pF2b <- ggplot(imp_df, aes(importance, variable, fill = hub)) + geom_col(colour = "black", width = 0.7) +
  scale_fill_manual(values = c("TRUE" = pal[1], "FALSE" = pal[2]), guide = "none") +
  labs(x = "Variable importance (impurity)", y = NULL) +
  theme_pub(7) + theme(axis.text.y = element_text(size = 6))
save_pub(pF2a, "ext_rsf_cindex", 3.4, 2.4)
save_pub(pF2b, "ext_rsf_importance", 3.4, 2.4)

# ------------------------------------------------------------
# 11. 组合图 Figure S20（2 列 × 4 行，A-H）
# ------------------------------------------------------------
log_line("=== 组装 Figure S21 ===")
comp <- plot_grid(pA, pB, pC, pD, pE, pF1, pF2a, pF2b, ncol = 2, nrow = 4,
                  labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
                  label_fontfamily = "Helvetica", label_fontface = "bold",
                  label_size = 9, label_colour = "black",
                  label_x = 0.02, label_y = 0.98, hjust = 0, vjust = 1)
W <- 7.0; H <- 11.0
pdf(file.path(OUT, "figures", "FigureS21_extended_prognostic_modeling.pdf"), width = W, height = H, family = "Helvetica")
print(comp); dev.off()
agg_png(file.path(OUT, "figures", "FigureS21_extended_prognostic_modeling.png"), width = W, height = H,
        units = "in", res = 300, background = "white")
print(comp); dev.off()
agg_tiff(file.path(OUT, "figures", "FigureS21_extended_prognostic_modeling.tiff"), width = W, height = H,
         units = "in", res = 300, compression = "lzw", background = "white")
print(comp); dev.off()

# 复制到投稿包
SUPP <- file.path(PKG, "figures/supplementary")
dir.create(SUPP, recursive = TRUE, showWarnings = FALSE)
for (f in c("FigureS21_extended_prognostic_modeling.pdf", "FigureS21_extended_prognostic_modeling.png",
            "FigureS21_extended_prognostic_modeling.tiff")) {
  file.copy(file.path(OUT, "figures", f), file.path(SUPP, f), overwrite = TRUE)
}

# ------------------------------------------------------------
# 12. 汇总
# ------------------------------------------------------------
summary_lines <- c(
  sprintf("LASSO 保留基因数：%d（其中 hub 基因 %d 个）",
          nrow(lasso_coef_tab), sum(lasso_coef_tab$gene %in% HUB)),
  sprintf("LASSO 选中基因：%s", paste(lasso_coef_tab$gene, collapse = ", ")),
  sprintf("C 指数（校正后）：%s",
          paste(sprintf("%s=%.3f", boot_tab$model, boot_tab$C_corrected), collapse = "; ")),
  sprintf("28 天 AUC：%s", paste(sprintf("%s=%.3f", perf$model, perf$AUC_28d), collapse = "; ")),
  sprintf("NRI（+LASSO vs Clinical）=%.3f（P=%.3g）；IDI=%.4f", res_cl$NRI, nri_p(res_cl), res_cl$IDI),
  sprintf("RSF OOB C 指数=%.3f", rsf_C))
log_line("=== 汇总 ===")
for (s in summary_lines) log_line(s)
writeLines(summary_lines, file.path(OUT, "tables", "ext_summary.txt"))
saveRDS(list(m_clin = m_clin, m_hub6z = m_hub6z, m_clin_hub6z = m_clin_hub6z,
             m_lasso = m_lasso, m_clin_lasso = m_clin_lasso, cv_lasso = cv_lasso, rsf = rsf),
        file.path(OUT, "models", "step106_models.rds"))
log_line("=== step106 完成 ===")
