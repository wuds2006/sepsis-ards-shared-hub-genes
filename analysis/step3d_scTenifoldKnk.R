# Step3d: scTenifoldKnk 虚拟敲除 (6个Hub基因)
suppressPackageStartupMessages({library(scTenifoldKnk); library(Matrix)})

setwd("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis")
OUT <- "08_singlecell/knockout"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# 读 raw counts
counts <- readMM("08_singlecell/cellchat_input/counts.mtx")  # 基因 x 细胞
genes <- readLines("08_singlecell/cellchat_input/genes.tsv")
rownames(counts) <- genes
cat("原始 counts:", dim(counts), "\n")

# 预过滤: top 1500 高表达基因 + 下采样 3000 细胞 (tensor分解内存O(genes^2*nNet))
gene_mean <- rowMeans(counts)
top_genes <- names(sort(gene_mean, decreasing = TRUE))[1:1500]
counts <- counts[top_genes, ]
set.seed(42)
cell_keep <- sample(seq_len(ncol(counts)), min(3000, ncol(counts)))
counts <- counts[, cell_keep]
counts <- counts[rowSums(counts) > 0, ]
cat("过滤后 counts:", dim(counts), "\n")

hub_genes <- c("UPP1")  # UPP1 是唯一在PBMC单核细胞中高表达的Hub基因
hub_genes <- hub_genes[hub_genes %in% rownames(counts)]
cat("可敲除 Hub 基因:", paste(hub_genes, collapse=", "), "\n")

# 对每个 Hub 基因做虚拟敲除(小规模网络, 避免OOM)
all_results <- list()
for (g in hub_genes) {
  cat(sprintf("\n===== 敲除 %s =====\n", g))
  t0 <- Sys.time()
  res <- tryCatch(
    scTenifoldKnk(countMatrix = counts, gKO = g, nc_nNet = 2, nc_nCells = 200, nCores = 1),
    error = function(e) { cat("敲除失败:", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) next
  cat("耗时:", round(as.numeric(difftime(Sys.time(), t0, units="mins")), 2), "分钟\n")

  # 受影响基因(dRegulation 距离)
  dr <- res$diffRegulation
  dr <- dr[order(dr$distance, decreasing = TRUE), ]
  write.csv(dr, file.path(OUT, paste0("KO_", g, "_diffRegulation.csv")), row.names = FALSE)
  all_results[[g]] <- dr
  cat(sprintf("受影响基因数: %d, Top5: %s\n", nrow(dr),
              paste(head(dr$gene, 5), collapse=", ")))
}

# 汇总: 每个基因敲除后的 Top 受影响基因
cat("\n===== 汇总: 各Hub基因敲除的Top10受影响基因 =====\n")
for (g in names(all_results)) {
  dr <- all_results[[g]]
  cat(sprintf("%-10s: %s\n", g, paste(head(dr$gene, 10), collapse=", ")))
}

saveRDS(all_results, file.path(OUT, "scTenifoldKnk_results.rds"))
cat("\n===== Step3d scTenifoldKnk 完成 =====\n")
