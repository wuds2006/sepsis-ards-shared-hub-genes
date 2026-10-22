# Step1d-SMR-HEIDI: HEIDI 异质性检验 (R 近似实现, 用已有 eQTL+GWAS 数据)
# 原理(Zhu et al. 2016): 检验 top cis-eQTL 是否通过单一因果变异解释 GWAS 信号
#   d_i = b_GWAS_i - b_SMR * b_eQTL_i, HEIDI统计量 = Σ d_i^2/var(d_i) ~ χ²(m)
# 局限: 无 LD 参考, 假设 cis-eQTL 近似独立(var(d_i) 忽略 LD 协方差项, 保守方向)
suppressPackageStartupMessages(library(data.table))

setwd("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis")

genes <- c("CDK5RAP2","CD177","GRB10","UPP1","VNN1")
outcomes <- c("sepsis","ards")
EQTL_P_THRESH <- 1.57e-3  # SMR 默认 HEIDI eQTL p 阈值

results <- data.frame()
for (g in genes) {
  for (o in outcomes) {
    f <- sprintf("04_MR_SMR/coloc_input/%s_%s_coloc.txt", g, o)
    if (!file.exists(f)) next
    d <- fread(f)
    d <- d[is.finite(e_beta) & is.finite(e_se) & is.finite(g_beta) & is.finite(g_se)]
    d <- d[!duplicated(SNP)]
    if (nrow(d) < 10) next

    # top cis-eQTL
    top <- d[which.min(e_p)]
    b_SMR <- top$g_beta / top$e_beta

    # 其他 cis-eQTL (eQTL p < 阈值, 排除 top)
    others <- d[e_p < EQTL_P_THRESH & SNP != top$SNP]
    if (nrow(others) < 1) {
      results <- rbind(results, data.frame(gene=g, outcome=o, topSNP=top$SNP,
        b_SMR=b_SMR, n_HEIDI=0, T=NA, p_HEIDI=NA, df=NA))
      next
    }

    # HEIDI 检验 (近似: 忽略 LD 协方差)
    d_i <- others$g_beta - b_SMR * others$e_beta
    var_d <- others$g_se^2 + b_SMR^2 * others$e_se^2
    T <- sum(d_i^2 / var_d)
    df <- nrow(others)
    p_HEIDI <- pchisq(T, df = df, lower.tail = FALSE)

    results <- rbind(results, data.frame(gene=g, outcome=o, topSNP=top$SNP,
      b_SMR=b_SMR, n_HEIDI=df, T=T, p_HEIDI=p_HEIDI, df=df))
    cat(sprintf("%-10s %-7s topSNP=%s b_SMR=%.3f n=%d T=%.1f p_HEIDI=%.4f\n",
        g, o, top$SNP, b_SMR, df, T, p_HEIDI))
  }
}

cat("\n===== HEIDI 检验结果 (p_HEIDI<0.05 表示存在多效性, SMR估计不可靠) =====\n")
print(results, row.names = FALSE)
write.csv(results, "04_MR_SMR/HEIDI_results.csv", row.names = FALSE)
cat("\n已保存 HEIDI_results.csv\n")
