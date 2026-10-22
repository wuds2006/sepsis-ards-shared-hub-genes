# ============================================================
# Step 48b: 反向 MR —— ARDS (FinnGen J10_ARDS) -> sepsis (ieu-b-4980)
# 工具变量: p<5e-6 (与正向 MR 同一放宽阈值), rsID 匹配, 1Mb 距离 clumping
# 流程镜像 step1d_MR_4980.R
# ============================================================
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "04_MR_SMR")

expo <- fread(file.path(OUT, "ards_sig.tsv"))            # FinnGen: beta 对应 alt
outc <- fread(file.path(OUT, "sepsis_match_ards.tsv"))   # Hamilton: beta 对应 effect_allele
cat("暴露工具变量:", nrow(expo), "| 结局匹配:", nrow(outc), "\n")

# 按 rsid 对齐
m <- merge(expo, outc[, .(rs_id, pval_o = p_value, beta_o = beta, sebeta_o = standard_error,
                          eaf_o = effect_allele_frequency, ea_o = effect_allele, oa_o = other_allele)],
           by.x = "rsids", by.y = "rs_id")
cat("对齐后 SNP 数:", nrow(m), "\n")

# 等位基因方向: sepsis effect_allele vs ARDS alt
m[, flip := (ea_o != alt)]
m[, beta_o := ifelse(flip, -beta_o, beta_o)]
m[, eaf_o := ifelse(flip, 1 - eaf_o, eaf_o)]

# 距离 clumping (1Mb, 用暴露 ARDS 坐标)
m[, chrom := as.numeric(`#chrom`)]
setorder(m, chrom, pos)
clumped <- m[0]; tmp <- m
while (nrow(tmp) > 0) {
  top <- tmp[1]; clumped <- rbind(clumped, top)
  tmp <- tmp[!(chrom == top$chrom & abs(pos - top$pos) < 1e6)]
}
cat("clumping 后工具变量数:", nrow(clumped), "\n")

dat <- data.frame(SNP = clumped$rsids,
                  effect_allele.exposure = clumped$alt, other_allele.exposure = clumped$ref,
                  beta.exposure = clumped$beta, se.exposure = clumped$sebeta,
                  eaf.exposure = clumped$af_alt, pval.exposure = clumped$pval,
                  effect_allele.outcome = clumped$ea_o, other_allele.outcome = clumped$oa_o,
                  beta.outcome = clumped$beta_o, se.outcome = clumped$sebeta_o,
                  eaf.outcome = clumped$eaf_o, pval.outcome = clumped$pval_o,
                  exposure = "ARDS", id.exposure = "J10_ARDS",
                  mr_keep.exposure = TRUE, pval_origin.exposure = "reported",
                  samplesize.exposure = 493732, ncase.exposure = 431,
                  outcome = "Sepsis_UKB", id.outcome = "ieu-b-4980",
                  mr_keep.outcome = TRUE, pval_origin.outcome = "reported",
                  samplesize.outcome = 462927, ncase.outcome = 11643,
                  stringsAsFactors = FALSE)
dat <- dat[complete.cases(dat), ]
dat <- dat[is.finite(dat$beta.exposure) & is.finite(dat$beta.outcome), ]
cat("可用 SNP 数:", nrow(dat), "\n")

if (nrow(dat) >= 2) {
  exposure_dat <- dat[, c("SNP","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure",
                          "eaf.exposure","pval.exposure","exposure","id.exposure","mr_keep.exposure",
                          "pval_origin.exposure","samplesize.exposure","ncase.exposure")]
  outcome_dat <- dat[, c("SNP","effect_allele.outcome","other_allele.outcome","beta.outcome","se.outcome",
                         "eaf.outcome","pval.outcome","outcome","id.outcome","mr_keep.outcome",
                         "pval_origin.outcome","samplesize.outcome","ncase.outcome")]
  harm <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
  cat("harmonise 后 SNP 数:", nrow(harm), " (keep:", sum(harm$mr_keep), ")\n")

  mr_res <- mr(harm, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
  write.csv(mr_res, file.path(OUT, "MR_reverse_ards_sepsis_results.csv"), row.names = FALSE)
  cat("\n===== 反向 MR 结果 (ARDS -> sepsis) =====\n")
  print(mr_res[, c("method","nsnp","b","se","pval")])

  het <- mr_heterogeneity(harm)
  ple <- mr_pleiotropy_test(harm)
  write.csv(het, file.path(OUT, "MR_reverse_heterogeneity.csv"), row.names = FALSE)
  write.csv(ple, file.path(OUT, "MR_reverse_pleiotropy.csv"), row.names = FALSE)
  cat("\n异质性(IVW Cochran Q):\n"); print(het[, c("method","Q","Q_pval")])
  cat("\n多效性(Egger intercept):\n"); print(ple[, c("egger_intercept","pval")])

  steiger <- directionality_test(harm)
  write.csv(steiger, file.path(OUT, "MR_reverse_steiger.csv"), row.names = FALSE)
  cat("\n===== Steiger =====\n")
  print(steiger[, c("snp_r2.exposure","snp_r2.outcome","correct_causal_direction","steiger_pval")])

  F_stat <- (harm$beta.exposure / harm$se.exposure)^2
  cat("\n工具变量 F 统计量: 均值", round(mean(F_stat), 1), " 最小", round(min(F_stat), 1), "\n")
} else {
  cat("\n可用工具变量不足 2 个, 反向 MR 不可行\n")
}
cat("\n===== Step 48b 反向 MR 完成 =====\n")
