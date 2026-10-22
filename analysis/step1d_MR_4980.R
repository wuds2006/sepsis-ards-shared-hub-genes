# ============================================================
# Step 1d-MR —— sepsis (ieu-b-4980 = GCST90281174, Hamilton UKB) -> ARDS (FinnGen)
# 工具变量: p<5e-6 (脓毒症无全基因组显著位点), rsID匹配, 距离clumping
# ============================================================
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(ggplot2)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "04_MR_SMR")

expo <- fread(file.path(OUT, "sepsis_sig_4980.tsv"))
outc <- fread(file.path(OUT, "ards_match_4980.tsv"))
cat("暴露工具变量:", nrow(expo), "| 结局匹配:", nrow(outc), "\n")

# 按 rsid 对齐
m <- merge(expo, outc[, .(rsids, pval_o = pval, beta_o = beta, sebeta_o = sebeta,
                          af_alt_o = af_alt, alt_o = alt, ref_o = ref)],
           by.x = "rsid", by.y = "rsids")
cat("对齐后 SNP 数:", nrow(m), "\n")

# 等位基因方向处理: sepsis ea vs ARDS alt
m[, flip := (ea != alt_o)]
m[, beta_o := ifelse(flip, -beta_o, beta_o)]
m[, af_alt_o := ifelse(flip, 1 - af_alt_o, af_alt_o)]

# 距离 clumping(1Mb, 用暴露 GRCh37 坐标)
setorder(m, chrom, pos)
clumped <- m[0]; tmp <- m
while (nrow(tmp) > 0) {
  top <- tmp[1]; clumped <- rbind(clumped, top)
  tmp <- tmp[!(chrom == top$chrom & abs(pos - top$pos) < 1e6)]
}
cat("clumping 后工具变量数:", nrow(clumped), "\n")

dat <- data.frame(SNP = clumped$rsid,
                  effect_allele.exposure = clumped$ea, other_allele.exposure = clumped$oa,
                  beta.exposure = clumped$beta, se.exposure = clumped$se,
                  eaf.exposure = clumped$eaf, pval.exposure = clumped$pval,
                  effect_allele.outcome = clumped$alt_o, other_allele.outcome = clumped$ref_o,
                  beta.outcome = clumped$beta_o, se.outcome = clumped$sebeta_o,
                  eaf.outcome = clumped$af_alt_o, pval.outcome = clumped$pval_o,
                  exposure = "Sepsis_UKB", id.exposure = "ieu-b-4980",
                  mr_keep.exposure = TRUE, pval_origin.exposure = "reported",
                  samplesize.exposure = 462927, ncase.exposure = 11643,
                  outcome = "ARDS", id.outcome = "J10_ARDS",
                  mr_keep.outcome = TRUE, pval_origin.outcome = "reported",
                  samplesize.outcome = 493732, ncase.outcome = 431,
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
  cat("harmonise 后 SNP 数:", nrow(harm), "\n")

  mr_res <- mr(harm, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
  mr_res$exposure <- "Sepsis_UKB(ieu-b-4980)"; mr_res$outcome <- "ARDS"
  write.csv(mr_res, file.path(OUT, "MR_sepsis4980_ARDS_results.csv"), row.names = FALSE)
  cat("\n===== MR 结果 =====\n")
  print(mr_res[, c("method","nsnp","b","se","pval")])

  # 异质性 + 多效性
  het <- mr_heterogeneity(harm)
  ple <- mr_pleiotropy_test(harm)
  write.csv(het, file.path(OUT, "MR_sepsis4980_heterogeneity.csv"), row.names = FALSE)
  write.csv(ple, file.path(OUT, "MR_sepsis4980_pleiotropy.csv"), row.names = FALSE)
  cat("\n异质性(IVW Cochran Q):\n"); print(het[, c("method","Q","Q_pval")])
  cat("\n多效性(Egger intercept):\n"); print(ple[, c("egger_intercept","pval")])

  # Steiger
  steiger <- directionality_test(harm)
  write.csv(steiger, file.path(OUT, "MR_sepsis4980_steiger.csv"), row.names = FALSE)
  cat("\n===== Steiger =====\n")
  print(steiger[, c("snp_r2.exposure","snp_r2.outcome","correct_causal_direction","steiger_pval")])

  # F 统计量(工具强度)
  F_stat <- (harm$beta.exposure / harm$se.exposure)^2
  cat("\n工具变量 F 统计量: 均值", round(mean(F_stat), 1), " 最小", round(min(F_stat), 1), "\n")

  tryCatch({ p <- mr_scatter_plot(mr_res, harm)[[1]]; save_fig(p, "MR_sepsis4980_scatter", w=7, h=6, dir=OUT) }, error=function(e) cat("散点图跳过\n"))
  tryCatch({ p2 <- mr_forest_plot(mr_singlesnp(harm))[[1]]; save_fig(p2, "MR_sepsis4980_forest", w=8, h=6, dir=OUT) }, error=function(e) cat("森林图跳过\n"))
}
cat("\n===== Step1d-MR (ieu-b-4980) 完成 =====\n")
