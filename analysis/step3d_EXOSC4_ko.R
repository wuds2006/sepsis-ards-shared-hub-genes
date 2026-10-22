# Step3d 补充: scTenifoldKnk 虚拟敲除 EXOSC4
# EXOSC4 是除 UPP1 外唯一在 PBMC 单核细胞有可检测表达的 Hub 基因(排名4492, 668细胞)
# 策略: top1500 高表达基因 + 强制纳入 EXOSC4, 保持内存可控
suppressPackageStartupMessages({library(scTenifoldKnk); library(Matrix)})

setwd("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis")
OUT <- "08_singlecell/knockout"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

counts <- readMM("08_singlecell/cellchat_input/counts.mtx")
genes <- readLines("08_singlecell/cellchat_input/genes.tsv")
rownames(counts) <- genes
cat("原始 counts:", dim(counts), "\n")

gene_mean <- rowMeans(counts)
top_genes <- names(sort(gene_mean, decreasing = TRUE))[1:1500]
# 强制纳入 EXOSC4
extra <- setdiff(c("UPP1","EXOSC4"), top_genes)
if (length(extra) > 0) {
  top_genes <- c(top_genes, extra)
  cat("强制纳入基因:", paste(extra, collapse=", "), "\n")
}
counts <- counts[top_genes, ]
set.seed(42)
cell_keep <- sample(seq_len(ncol(counts)), min(3000, ncol(counts)))
counts <- counts[, cell_keep]
counts <- counts[rowSums(counts) > 0, ]
cat("过滤后 counts:", dim(counts), "\n")

g <- "EXOSC4"
if (!(g %in% rownames(counts))) stop("EXOSC4 不在过滤后矩阵中")
cat(sprintf("\n===== 敲除 %s =====\n", g))
t0 <- Sys.time()
res <- tryCatch(
  scTenifoldKnk(countMatrix = counts, gKO = g, nc_nNet = 2, nc_nCells = 200, nCores = 1,
                qc_minPCT = 0.02),
  error = function(e) { cat("敲除失败:", conditionMessage(e), "\n"); NULL })
if (!is.null(res)) {
  cat("耗时:", round(as.numeric(difftime(Sys.time(), t0, units="mins")), 2), "分钟\n")
  dr <- res$diffRegulation
  dr <- dr[order(dr$distance, decreasing = TRUE), ]
  write.csv(dr, file.path(OUT, paste0("KO_", g, "_diffRegulation.csv")), row.names = FALSE)
  cat(sprintf("受影响基因数: %d, Top10: %s\n", nrow(dr), paste(head(dr$gene, 10), collapse=", ")))
  cat("===== EXOSC4 敲除完成 =====\n")
} else {
  cat("===== EXOSC4 敲除失败 =====\n")
}
