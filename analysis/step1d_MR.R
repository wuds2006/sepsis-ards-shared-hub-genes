# ============================================================
# Step 1d-MR —— 两样本 MR: 细菌感染(脓毒症代理) -> ARDS + Steiger
# 暴露: FinnGen R11 AB1_OTHER_BACTERIAL; 结局: FinnGen R11 J10_ARDS
# 内存优化: awk 流式预过滤, 只读显著 SNP
# 注: 暴露与结局同源(FinnGen), 存在样本重叠
# ============================================================
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(ggplot2)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "04_MR_SMR")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

EXPO <- file.path(BASE, "00_rawdata/gwas/finngen_R11_AB1_OTHER_BACTERIAL.gz")
OUTC <- file.path(BASE, "00_rawdata/gwas/finngen_R11_J10_ARDS.gz")
# 列: #chrom(1) pos(2) ref(3) alt(4) rsids(5) nearest_genes(6) pval(7) mlogp(8) beta(9) sebeta(10) af_alt(11)

# 1) 暴露显著 SNP (p<5e-8) 流式提取
cat("提取暴露显著 SNP...\n")
expo <- fread(cmd = sprintf("gunzip -c %s | awk -F'\\t' 'NR==1 || $7+0 < 5e-8'", shQuote(EXPO)),
              header = TRUE, showProgress = FALSE)
cat("暴露显著 SNP 数:", nrow(expo), "\n")
if (nrow(expo) < 2) {
  expo <- fread(cmd = sprintf("gunzip -c %s | awk -F'\\t' 'NR==1 || $7+0 < 5e-6'", shQuote(EXPO)),
                header = TRUE, showProgress = FALSE)
  cat("放宽到 5e-6, SNP 数:", nrow(expo), "\n")
}

# 2) 距离法 clumping (1Mb)
names(expo) <- sub("^#", "", names(expo))
setorder(expo, chrom, pos)
clumped <- expo[0]
tmp <- expo
while (nrow(tmp) > 0) {
  top <- tmp[1]
  clumped <- rbind(clumped, top)
  tmp <- tmp[!(chrom == top$chrom & abs(pos - top$pos) < 1e6)]
}
cat("clumping 后工具变量数:", nrow(clumped), "\n")

# 3) 结局: 用 rsids 匹配(单次 awk 流式)
rs <- clumped$rsids
rs <- rs[rs != ""]
grep_pat <- paste(rs, collapse = "|")
outc <- fread(cmd = sprintf("gunzip -c %s | awk -F'\\t' 'NR==1 || $5 ~ /%s/'", shQuote(OUTC), grep_pat),
              header = TRUE, showProgress = FALSE)
names(outc) <- sub("^#", "", names(outc))
cat("结局匹配 SNP 数:", nrow(outc), "\n")

# 4) 对齐
expo_use <- clumped[rsids %in% outc$rsids]
outc[, c("pval_o","beta_o","sebeta_o","af_alt_o") := .(pval, beta, sebeta, af_alt)]
m <- merge(expo_use, outc[, .(rsids, beta_o, sebeta_o, af_alt_o, pval_o)], by = "rsids")
cat("对齐后 SNP 数:", nrow(m), "\n")

if (nrow(m) >= 1) {
  # 5) harmonise + MR
  dat <- data.frame(SNP = m$rsids,
                    effect_allele.exposure = m$alt, other_allele.exposure = m$ref,
                    beta.exposure = m$beta, se.exposure = m$sebeta,
                    eaf.exposure = m$af_alt, pval.exposure = m$pval,
                    effect_allele.outcome = m$alt, other_allele.outcome = m$ref,
                    beta.outcome = m$beta_o, se.outcome = m$sebeta_o,
                    eaf.outcome = m$af_alt_o, pval.outcome = m$pval_o)
  dat <- harmonise_data(exposure_dat = dat, outcome_dat = dat, action = 2)

  mr_res <- mr(dat, method_list = c("mr_ivw","mr_egger_regression",
                                    "mr_weighted_median","mr_weighted_mode"))
  mr_res$exposure <- "Bacterial_infection(sepsis_proxy)"; mr_res$outcome <- "ARDS"
  write.csv(mr_res, file.path(OUT, "MR_sepsis_ARDS_results.csv"), row.names = FALSE)
  print(mr_res[, c("method","nsnp","b","se","pval")])

  steiger <- directionality_test(dat)
  write.csv(steiger, file.path(OUT, "MR_steiger.csv"), row.names = FALSE)
  cat("\nSteiger:\n"); print(steiger[, c("snp_r2.exposure","snp_r2.outcome","correct_causal_direction","steiger_pval")])

  if (nrow(dat) > 1) {
    p <- mr_scatter_plot(mr_res, dat)[[1]]
    save_fig(p, "MR_scatter", w = 7, h = 6, dir = OUT)
  }
}
cat("\n===== Step1d-MR 完成 =====\n")
