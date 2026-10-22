# ============================================================
# Step 51 —— MR 统计功效分析（审稿 W3 / 问题2）
# 1) 正向 sepsis->ARDS 与反向 ARDS->sepsis 的 IV F 统计量、R2
# 2) 解析功效（Brion 2013 / mRnd 等价）：给定 OR 的把握度
# 3) 模拟功效验证（5000 次）
# 4) 80% 把握度下最小可检 OR
# ============================================================
suppressPackageStartupMessages({ library(data.table) })
BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "04_MR_SMR")

# ---------- 重建 harmonised IV（镜像 step1d_MR_4980.R / step48b_reverse_MR.R） ----------
build_forward <- function() {
  expo <- fread(file.path(OUT, "sepsis_sig_4980.tsv"))
  outc <- fread(file.path(OUT, "ards_match_4980.tsv"))
  m <- merge(expo, outc[, .(rsids, beta_o = beta, sebeta_o = sebeta, alt_o = alt)],
             by.x = "rsid", by.y = "rsids")
  m[, flip := (ea != alt_o)]
  m[, beta_o := ifelse(flip, -beta_o, beta_o)]
  setorder(m, chrom, pos)
  clumped <- m[0]; tmp <- m
  while (nrow(tmp) > 0) {
    top <- tmp[1]; clumped <- rbind(clumped, top)
    tmp <- tmp[!(chrom == top$chrom & abs(pos - top$pos) < 1e6)]
  }
  data.frame(snp = clumped$rsid, bx = clumped$beta, sex = clumped$se,
             eaf = clumped$eaf, by = clumped$beta_o, sey = clumped$sebeta_o,
             stringsAsFactors = FALSE)
}
build_reverse <- function() {
  expo <- fread(file.path(OUT, "ards_sig.tsv"))
  outc <- fread(file.path(OUT, "sepsis_match_ards.tsv"))
  m <- merge(expo, outc[, .(rs_id, beta_o = beta, sebeta_o = standard_error, ea_o = effect_allele)],
             by.x = "rsids", by.y = "rs_id")
  m[, flip := (ea_o != alt)]
  m[, beta_o := ifelse(flip, -beta_o, beta_o)]
  m[, chrom := as.numeric(`#chrom`)]
  setorder(m, chrom, pos)
  clumped <- m[0]; tmp <- m
  while (nrow(tmp) > 0) {
    top <- tmp[1]; clumped <- rbind(clumped, top)
    tmp <- tmp[!(chrom == top$chrom & abs(pos - top$pos) < 1e6)]
  }
  data.frame(snp = clumped$rsids, bx = clumped$beta, sex = clumped$sebeta,
             eaf = clumped$af_alt, by = clumped$beta_o, sey = clumped$sebeta_o,
             stringsAsFactors = FALSE)
}
clean <- function(d) { d <- d[complete.cases(d) & d$sex > 0 & d$sey > 0, ]; d }
# 精确镜像 TwoSampleMR harmonise_data 的移除清单（等位基因不兼容 SNP）
drop_incompatible <- list(
  forward_sepsis_to_ARDS = "rs7103228",
  reverse_ARDS_to_sepsis = c("rs11392863", "rs138100165"))

ivw_se <- function(d) sqrt(1 / sum(d$bx^2 / d$sey^2))
# TwoSampleMR mr_ivw 默认乘法随机效应 SE：fixed SE × sqrt(max(1, Q/df))
ivw_se_mult <- function(d) {
  w  <- d$bx^2 / d$sey^2
  b  <- sum(d$bx * d$by / d$sey^2) / sum(w)
  Q  <- sum((d$by - b * d$bx)^2 / d$sey^2)
  se_f <- sqrt(1 / sum(w))
  se_f * sqrt(max(1, Q / (nrow(d) - 1)))
}
power_analytic <- function(beta_true, se_ivw, alpha = 0.05) {
  stats::pnorm(abs(beta_true) / se_ivw - stats::qnorm(1 - alpha/2))
}
power_sim <- function(d, beta_true, nsim = 5000, alpha = 0.05) {
  set.seed(11)
  zcrit <- stats::qnorm(1 - alpha/2); hit <- 0
  w  <- d$bx / d$sey^2                       # IVW 权重 bx/sigma_Y^2
  sw <- sum(d$bx * w)                        # 分母 Σ bx^2/sey^2
  for (i in 1:nsim) {
    by_sim <- stats::rnorm(nrow(d), mean = beta_true * d$bx, sd = d$sey)
    b_hat  <- sum(w * by_sim) / sw
    se_hat <- sqrt(1 / sw)
    if (abs(b_hat) / se_hat > zcrit) hit <- hit + 1
  }
  hit / nsim
}

summarize <- function(tag, d) {
  d <- clean(d)
  d <- d[!(d$snp %in% drop_incompatible[[tag]]), ]   # 镜像 harmonise 移除
  F <- (d$bx / d$sex)^2
  r2 <- 2 * d$eaf * (1 - d$eaf) * d$bx^2
  se_ivw   <- ivw_se(d)        # 固定效应（mRnd/Brion 解析式用）
  se_mult  <- ivw_se_mult(d)   # 乘法随机效应（与正文报告口径一致）
  pow_mult <- function(ors, se) sapply(log(ors), function(b) power_analytic(b, se))
  beta_min_f  <- (stats::qnorm(0.975) + stats::qnorm(0.80)) * se_ivw
  beta_min_m  <- (stats::qnorm(0.975) + stats::qnorm(0.80)) * se_mult
  ors <- c(1.1, 1.2, 1.5, 2.0)
  pow_ana <- pow_mult(ors, se_ivw)
  pow_sim <- sapply(log(ors), function(b) power_sim(d, b, nsim = 5000))
  pow_rpt <- pow_mult(ors, se_mult)   # 按正文报告口径（乘法 SE）
  res <- data.frame(
    direction = tag, n_IV = nrow(d),
    F_mean = round(mean(F), 1), F_median = round(median(F), 1),
    F_min = round(min(F), 1), F_max = round(max(F), 1),
    total_R2 = round(sum(r2), 6),
    IVW_SE_logOR_fixed = round(se_ivw, 4),
    IVW_SE_logOR_multiplicative = round(se_mult, 4),
    min_detectable_OR_80power_fixed = round(exp(beta_min_f), 2),
    min_detectable_OR_80power_as_reported = round(exp(beta_min_m), 2))
  cat(sprintf("\n===== %s: %d IVs, F mean %.1f, SE(logOR) fixed=%.4f / mult=%.4f =====\n",
              tag, nrow(d), mean(F), se_ivw, se_mult))
  cat(sprintf("  80%%功效最小可检 OR: 固定效应 %.2f | 按报告口径(乘法) %.2f\n",
              exp(beta_min_f), exp(beta_min_m)))
  tab <- data.frame(OR = ors, power_fixed = round(pow_ana, 3),
                    power_simulated = round(pow_sim, 3),
                    power_as_reported_SE = round(pow_rpt, 3))
  print(tab)
  write.csv(d, file.path(OUT, sprintf("power_%s_IVs.csv", tag)), row.names = FALSE)
  list(summary = res, power_tab = tab)
}

fw <- summarize("forward_sepsis_to_ARDS", build_forward())
rv <- summarize("reverse_ARDS_to_sepsis", build_reverse())

sm <- rbind(fw$summary, rv$summary)
pt <- rbind(data.frame(direction = "forward_sepsis_to_ARDS", fw$power_tab),
            data.frame(direction = "reverse_ARDS_to_sepsis", rv$power_tab))
write.csv(sm, file.path(OUT, "MR_power_summary.csv"), row.names = FALSE)
write.csv(pt, file.path(OUT, "MR_power_by_OR.csv"), row.names = FALSE)
cat("\n=== Step 51 完成 ===\n")
