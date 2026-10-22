# 补装缺失 R 包（不污染既有环境，仅追加缺失项）
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 600)

need_cran <- c("fastshap","survminer","VennDiagram","UpSetR","glmnetUtils","naivebayes","gbm","patchwork")
for (p in need_cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("安装 CRAN:", p, "\n")
    tryCatch(install.packages(p, quiet = TRUE), error = function(e) cat("  失败:", p, conditionMessage(e), "\n"))
  } else cat("已有:", p, "\n")
}

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", quiet = TRUE)
need_bioc <- c("DESeq2","sva","edgeR","MCPcounter")
for (p in need_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("安装 Bioc:", p, "\n")
    tryCatch(BiocManager::install(p, ask = FALSE, update = FALSE), error = function(e) cat("  失败:", p, conditionMessage(e), "\n"))
  } else cat("已有:", p, "\n")
}

cat("\n=== 最终状态 ===\n")
for (p in c(need_cran, need_bioc)) cat(sprintf("%-15s %s\n", p, requireNamespace(p, quietly = TRUE)))
cat("DONE\n")
