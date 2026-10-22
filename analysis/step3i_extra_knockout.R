# ============================================================
# Step 3i-extra —— EXOSC4 / VNN1 补敲（下调 qc_minPCT）
# ============================================================
suppressPackageStartupMessages({library(scTenifoldKnk); library(Matrix); library(ggplot2); library(ragg)})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "08_singlecell/GSE167363/knockout")

save_fig3 <- function(p, name, w = 7, h = 6, dir = ".", line_art = FALSE) {
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf"))
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300))
  tdpi <- if (line_art) 1000 else 500
  agg_tiff(file.path(dir, paste0(name, ".tiff")), width = w, height = h,
           units = "in", res = tdpi, compression = "lzw")
  print(p); invisible(dev.off())
  invisible(p)
}

cnt <- readMM(file.path(BASE, "08_singlecell/GSE167363/input/counts.mtx"))
genes <- readLines(file.path(BASE, "08_singlecell/GSE167363/input/genes.tsv"))
rownames(cnt) <- genes

HUB <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")
gene_mean <- rowMeans(cnt)
top_genes <- names(sort(gene_mean, decreasing = TRUE))[1:1500]
counts <- cnt[union(top_genes, HUB), ]
set.seed(42)
cell_keep <- sample(seq_len(ncol(counts)), min(3000, ncol(counts)))
counts <- counts[, cell_keep]
counts <- counts[rowSums(counts) > 0, ]
cat("counts:", dim(counts), "\n")

targets <- c(EXOSC4 = 0.01, VNN1 = 0.005)
for (g in names(targets)) {
  csv_out <- file.path(OUT, paste0("KO_", g, "_diffRegulation.csv"))
  if (file.exists(csv_out)) { cat(g, "已存在, 跳过\n"); next }
  cat(sprintf("\n===== 敲除 %s (minPCT=%.3f) =====\n", g, targets[g]))
  t0 <- Sys.time()
  res <- tryCatch(
    scTenifoldKnk(countMatrix = counts, gKO = g, qc_minPCT = targets[g],
                  nc_nNet = 2, nc_nCells = 200, nCores = 1),
    error = function(e) { cat("敲除失败:", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) next
  cat("耗时:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "分钟\n")
  dr <- res$diffRegulation
  dr <- dr[order(dr$distance, decreasing = TRUE), ]
  write.csv(dr, csv_out, row.names = FALSE)
  sig <- subset(dr, p.adj < 0.05)
  cat(sprintf("受影响基因数: %d | 显著(FDR<0.05): %d | Top5: %s\n",
              nrow(dr), nrow(sig), paste(head(dr$gene, 5), collapse = ", ")))
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

cat("\n===== Step3i-extra 完成 =====\n")
