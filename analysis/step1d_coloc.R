# Step1d-coloc：eQTL × GWAS 贝叶斯共定位
suppressPackageStartupMessages({library(coloc); library(data.table)})

setwd("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/04_MR_SMR")

genes <- c("CDK5RAP2","CD177","GRB10","UPP1","VNN1")
outcomes <- c("sepsis","ards")

results <- data.frame()
for (g in genes) {
  for (o in outcomes) {
    f <- sprintf("coloc_input/%s_%s_coloc.txt", g, o)
    if (!file.exists(f)) next
    d <- fread(f)
    d <- d[is.finite(e_beta) & is.finite(e_se) & is.finite(g_beta) & is.finite(g_se)]
    d <- d[e_maf > 0.001 & e_maf < 0.999 & g_maf > 0.001 & g_maf < 0.999]
    d <- d[!duplicated(SNP)]
    if (nrow(d) < 100) { cat(g, o, "SNP太少, 跳过\n"); next }

    D1 <- list(beta = d$e_beta, varbeta = d$e_se^2, MAF = d$e_maf,
               N = 30596, type = "quant", snp = d$SNP)
    if (o == "sepsis") {
      N <- 486484; s <- 11643/486484
    } else {
      N <- 493732; s <- 431/493732
    }
    D2 <- list(beta = d$g_beta, varbeta = d$g_se^2, MAF = d$g_maf,
               N = N, type = "cc", s = s, snp = d$SNP)

    res <- tryCatch(coloc.abf(D1, D2), error = function(e) NULL)
    if (is.null(res)) { cat(g, o, "coloc 失败\n"); next }
    s <- res$summary
    PP.H0 <- as.numeric(s["PP.H0.abf"]); PP.H1 <- as.numeric(s["PP.H1.abf"])
    PP.H2 <- as.numeric(s["PP.H2.abf"]); PP.H3 <- as.numeric(s["PP.H3.abf"])
    PP.H4 <- as.numeric(s["PP.H4.abf"])
    results <- rbind(results, data.frame(gene = g, outcome = o, nsnp = nrow(d),
        PP.H0 = PP.H0, PP.H1 = PP.H1, PP.H2 = PP.H2, PP.H3 = PP.H3, PP.H4 = PP.H4))
    cat(sprintf("%-10s %-7s nsnp=%d  PP.H4=%.3f  PP.H3=%.3f\n", g, o, nrow(d),
        PP.H4, PP.H3))
  }
}

cat("\n===== coloc 共定位结果 =====\n")
print(results, row.names = FALSE)
write.csv(results, "coloc_results.csv", row.names = FALSE)
cat("\n已保存 coloc_results.csv\n")
