# ============================================================
# Step 1d-MR —— 读取提取文件, 按 chrom:pos 对齐, MR + Steiger
# ============================================================
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(ggplot2)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "04_MR_SMR")

expo <- fread(file.path(OUT, "expo_sig.tsv"), header = TRUE)
outc <- fread(file.path(OUT, "outc_match.tsv"), header = TRUE)
cat("暴露显著 SNP:", nrow(expo), "| 结局匹配 SNP:", nrow(outc), "\n")

# 距离法 clumping (1Mb)
expo[, chrom := as.integer(chrom)]
setorder(expo, chrom, pos)
clumped <- expo[0]; tmp <- expo
while (nrow(tmp) > 0) {
  top <- tmp[1]; clumped <- rbind(clumped, top)
  tmp <- tmp[!(chrom == top$chrom & abs(pos - top$pos) < 1e6)]
}
cat("clumping 后工具变量数:", nrow(clumped), "\n")

# 按 chrom:pos 对齐
expo_use <- clumped
outc[, chrom := as.integer(chrom)]
m <- merge(expo_use, outc, by = c("chrom","pos"), suffixes = c("", "_o"))
cat("对齐后 SNP 数:", nrow(m), "\n")

if (nrow(m) >= 1) {
  # 处理链方向(ref/alt 一致则直接用, 不一致翻转 beta)
  m[, flip := (ref != ref_o) | (alt != alt_o)]
  m[, beta_o := ifelse(flip, -beta_o, beta_o)]
  # 交换 af
  m[, af_alt_o := ifelse(flip, 1 - af_alt_o, af_alt_o)]

  dat <- data.frame(SNP = paste0("chr", m$chrom, "_", m$pos),
                    effect_allele.exposure = m$alt, other_allele.exposure = m$ref,
                    beta.exposure = m$beta, se.exposure = m$sebeta,
                    eaf.exposure = m$af_alt, pval.exposure = m$pval,
                    effect_allele.outcome = m$alt, other_allele.outcome = m$ref,
                    beta.outcome = m$beta_o, se.outcome = m$sebeta_o,
                    eaf.outcome = m$af_alt_o, pval.outcome = m$pval_o,
                    exposure = "Bacterial_infection", id.exposure = "AB1_OTHER_BACTERIAL",
                    mr_keep.exposure = TRUE, pval_origin.exposure = "reported",
                    samplesize.exposure = 360000, ncase.exposure = 20000,
                    outcome = "ARDS", id.outcome = "J10_ARDS",
                    mr_keep.outcome = TRUE, pval_origin.outcome = "reported",
                    samplesize.outcome = 493732, ncase.outcome = 431,
                    stringsAsFactors = FALSE)
  dat <- dat[complete.cases(dat), ]
  dat <- dat[is.finite(dat$beta.exposure) & is.finite(dat$beta.outcome), ]
  cat("可用 SNP 数:", nrow(dat), "\n")

  if (nrow(dat) >= 2) {
    exposure_dat <- dat[, c("SNP","effect_allele.exposure","other_allele.exposure",
                            "beta.exposure","se.exposure","eaf.exposure","pval.exposure",
                            "exposure","id.exposure","mr_keep.exposure","pval_origin.exposure")]
    outcome_dat <- dat[, c("SNP","effect_allele.outcome","other_allele.outcome",
                           "beta.outcome","se.outcome","eaf.outcome","pval.outcome",
                           "outcome","id.outcome","mr_keep.outcome","pval_origin.outcome")]
    harm <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
    cat("harmonise 后 SNP 数:", nrow(harm), "\n")
    mr_res <- mr(harm, method_list = c("mr_ivw","mr_egger_regression",
                                       "mr_weighted_median","mr_weighted_mode"))
    mr_res$exposure <- "Bacterial_infection(sepsis_proxy)"; mr_res$outcome <- "ARDS"
    write.csv(mr_res, file.path(OUT, "MR_sepsis_ARDS_results.csv"), row.names = FALSE)
    cat("\n===== MR 结果 =====\n")
    print(mr_res[, c("method","nsnp","b","se","pval")])

    steiger <- directionality_test(harm)
    write.csv(steiger, file.path(OUT, "MR_steiger.csv"), row.names = FALSE)
    cat("\n===== Steiger =====\n")
    print(steiger[, c("snp_r2.exposure","snp_r2.outcome","correct_causal_direction","steiger_pval")])

    tryCatch({
      p <- mr_scatter_plot(mr_res, harm)[[1]]; save_fig(p, "MR_scatter", w=7, h=6, dir=OUT)
    }, error = function(e) cat("散点图跳过\n"))
    tryCatch({
      p2 <- mr_forest_plot(mr_singlesnp(harm))[[1]]; save_fig(p2, "MR_forest", w=8, h=6, dir=OUT)
    }, error = function(e) cat("森林图跳过\n"))
  } else cat("可用 SNP 过少(<2)\n")
}
cat("\n===== Step1d-MR 完成 =====\n")
