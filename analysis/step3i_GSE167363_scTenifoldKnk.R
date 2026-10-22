# ============================================================
# Step 3i —— GSE167363 scTenifoldKnk 虚拟敲除（hub 基因）
# 敲除候选: UPP1(42.7%), CDK5RAP2(5.8%), EXOSC4(3.0%), VNN1(1.1%, 尝试)
# GRB10(0.6%)/CD177(0.4%) 表达过低, 不可敲除 -> 记录
# ============================================================
suppressPackageStartupMessages({library(scTenifoldKnk); library(Matrix); library(ggplot2); library(ragg)})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "08_singlecell/GSE167363/knockout")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

save_fig3 <- function(p, name, w = 7, h = 6, dir = ".", line_art = FALSE) {
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf"))
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300))
  tdpi <- if (line_art) 1000 else 500
  agg_tiff(file.path(dir, paste0(name, ".tiff")), width = w, height = h,
           units = "in", res = tdpi, compression = "lzw")
  print(p); invisible(dev.off())
  invisible(p)
}

# 读 raw counts（step3f 导出）
cnt <- readMM(file.path(BASE, "08_singlecell/GSE167363/input/counts.mtx"))
genes <- readLines(file.path(BASE, "08_singlecell/GSE167363/input/genes.tsv"))
rownames(cnt) <- genes
cat("原始 counts:", dim(cnt), "\n")

# 预过滤: top 1500 高表达基因 + 强制保留 hub 候选
HUB <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")
gene_mean <- rowMeans(cnt)
top_genes <- names(sort(gene_mean, decreasing = TRUE))[1:1500]
counts <- cnt[union(top_genes, HUB), ]

# 表达统计
pct <- sapply(HUB, function(g) if (g %in% rownames(cnt)) mean(cnt[g, ] > 0) else 0)
cat("\nhub 基因表达细胞占比:\n"); print(round(pct, 3))

# 敲除策略: 占比 >= 1% 尝试敲除
knockable <- HUB[pct >= 0.01]
not_knockable <- HUB[pct < 0.01]
cat("可敲除:", paste(knockable, collapse = ", "), "\n")
cat("不可敲除(表达<1%):", paste(not_knockable, collapse = ", "), "\n")

# 下采样 3000 细胞
set.seed(42)
cell_keep <- sample(seq_len(ncol(counts)), min(3000, ncol(counts)))
counts <- counts[, cell_keep]
counts <- counts[rowSums(counts) > 0, ]
cat("过滤后 counts:", dim(counts), "\n")

all_results <- list()
for (g in knockable) {
  if (!(g %in% rownames(counts))) { cat(g, "不在过滤后矩阵, 跳过\n"); next }
  csv_out <- file.path(OUT, paste0("KO_", g, "_diffRegulation.csv"))
  if (file.exists(csv_out)) {
    cat(sprintf("跳过 %s (已存在结果)\n", g))
    all_results[[g]] <- read.csv(csv_out, stringsAsFactors = FALSE)
    next
  }
  cat(sprintf("\n===== 敲除 %s =====\n", g))
  t0 <- Sys.time()
  res <- tryCatch(
    scTenifoldKnk(countMatrix = counts, gKO = g, nc_nNet = 2, nc_nCells = 200, nCores = 1),
    error = function(e) { cat("敲除失败:", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) next
  cat("耗时:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "分钟\n")
  dr <- res$diffRegulation
  dr <- dr[order(dr$distance, decreasing = TRUE), ]
  write.csv(dr, csv_out, row.names = FALSE)
  all_results[[g]] <- dr
  sig <- subset(dr, p.adj < 0.05)
  cat(sprintf("受影响基因数: %d | 显著(FDR<0.05): %d | Top5: %s\n",
              nrow(dr), nrow(sig), paste(head(dr$gene, 5), collapse = ", ")))

  # Top10 柱状图
  top10 <- head(dr, 10)
  top10$gene <- factor(top10$gene, levels = rev(top10$gene))
  top10$dir <- ifelse(top10$distance > 0, "up", "down")
  pb <- ggplot(top10, aes(x = distance, y = gene, fill = dir)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(up = "#D99AAE", down = "#8FA2D6")) +
    ggtitle(paste0("scTenifoldKnk virtual KO: ", g, " (Top10 perturbed genes)")) +
    xlab("Regulatory distance (dRegulation)") + ylab(NULL) +
    theme_pub() + theme(legend.position = "none")
  save_fig3(pb, paste0("KO_", g, "_top10"), w = 6.5, h = 4.5, dir = OUT, line_art = TRUE)
}

# 为跳过的基因补画 Top10 图
for (g in names(all_results)) {
  if (file.exists(file.path(OUT, paste0("KO_", g, "_top10.png")))) next
  dr <- all_results[[g]]
  top10 <- head(dr, 10)
  top10$gene <- factor(top10$gene, levels = rev(top10$gene))
  top10$dir <- ifelse(top10$distance > 0, "up", "down")
  pb <- ggplot(top10, aes(x = distance, y = gene, fill = dir)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(up = "#D99AAE", down = "#8FA2D6")) +
    ggtitle(paste0("scTenifoldKnk virtual KO: ", g, " (Top10 perturbed genes)")) +
    xlab("Regulatory distance (dRegulation)") + ylab(NULL) +
    theme_pub() + theme(legend.position = "none")
  save_fig3(pb, paste0("KO_", g, "_top10"), w = 6.5, h = 4.5, dir = OUT, line_art = TRUE)
}

# 汇总
cat("\n===== 汇总: 各 hub 基因敲除 Top10 受影响基因 =====\n")
for (g in names(all_results)) {
  dr <- all_results[[g]]
  cat(sprintf("%-10s: %s\n", g, paste(head(dr$gene, 10), collapse = ", ")))
}

# 敲除结果交叉分析: 共同受影响基因
if (length(all_results) >= 2) {
  sig_lists <- lapply(all_results, function(dr) subset(dr, p.adj < 0.05)$gene)
  common_genes <- Reduce(intersect, sig_lists)
  cat("\n各敲除的显著受影响基因交集:", length(common_genes), "个\n")
  write.csv(data.frame(gene = common_genes),
            file.path(OUT, "KO_common_affected_genes.csv"), row.names = FALSE)

  # 各基因敲除后对每个 hub 的影响矩阵（检查 hub 间相互影响）
  hub_eff <- matrix(NA, nrow = length(all_results), ncol = length(HUB),
                    dimnames = list(names(all_results), HUB))
  for (g in names(all_results)) {
    dr <- all_results[[g]]
    idx <- match(HUB, dr$gene)
    hub_eff[g, ] <- ifelse(is.na(idx), NA, dr$distance[idx])
  }
  write.csv(hub_eff, file.path(OUT, "KO_hub_cross_effects.csv"))
  cat("\nhub 基因间敲除交叉影响(dRegulation 距离):\n"); print(round(hub_eff, 3))
}

# 敲除可行性记录
feas <- data.frame(gene = HUB, pct_expressing = round(pct, 4),
                   knockable = HUB %in% names(all_results))
write.csv(feas, file.path(OUT, "KO_feasibility.csv"), row.names = FALSE)

saveRDS(all_results, file.path(OUT, "scTenifoldKnk_results.rds"))
cat("\n===== Step3i scTenifoldKnk (GSE167363) 完成 =====\n")
