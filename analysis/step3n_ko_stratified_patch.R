# ============================================================
# Step 3n-3b —— 补跑低表达 hub 基因的分组 scTenifoldKnk（nc_nCells=500）
# 原因: nc_nCells=200 时 EXOSC4/VNN1/CDK5RAP2 在组内子集中被 QC 剔除
# 策略: 这 3 个基因的所有分组统一用 nc_nNet=2, nc_nCells=500 重跑；
#       UPP1（高表达）保留 nc_nCells=200 的既有结果（已成功）
# ============================================================
suppressPackageStartupMessages({library(scTenifoldKnk); library(Matrix)})
set.seed(42)
BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT <- file.path(BASE, "08_singlecell/knockout_stratified")
HUB <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")

load_counts <- function(which) {
  if (which == "GSE151263") cd <- file.path(BASE, "08_singlecell/cellchat_input") else
    cd <- file.path(BASE, "08_singlecell/GSE167363/input")
  cnt <- readMM(file.path(cd, "counts.mtx"))
  genes <- sub("\t.*$", "", readLines(file.path(cd, "genes.tsv")))
  rownames(cnt) <- genes
  meta <- read.csv(file.path(cd, "meta.csv"), stringsAsFactors = FALSE)
  cell_col <- if ("cell" %in% names(meta)) "cell" else names(meta)[1]
  rownames(meta) <- meta[[cell_col]]
  stopifnot(ncol(cnt) == nrow(meta))
  colnames(cnt) <- rownames(meta)
  list(cnt = cnt, meta = meta)
}

RUNS <- list(
  c("GSE151263", "ARDS",  "EXOSC4"),
  c("GSE151263", "Sepsis", "EXOSC4"),
  c("GSE167363", "Sepsis", "EXOSC4"),
  c("GSE167363", "HC",     "EXOSC4"),
  c("GSE167363", "Sepsis", "VNN1"),
  c("GSE167363", "HC",     "VNN1"),
  c("GSE167363", "Sepsis", "CDK5RAP2"),
  c("GSE167363", "HC",     "CDK5RAP2")
)

for (run in RUNS) {
  which <- run[1]; grp <- run[2]; g <- run[3]
  csv_out <- file.path(OUT, sprintf("KO_%s_%s_%s_diffRegulation.csv", which, grp, g))
  if (file.exists(csv_out)) {
    # PBMC Sepsis EXOSC4 已有 200-cells 版本 -> 删除统一用 500 重跑
    if (!(which == "GSE167363" && grp == "HC" && g == "CDK5RAP2" && FALSE)) {
      cat(sprintf("删除旧结果以统一参数: %s\n", basename(csv_out)))
      file.remove(csv_out)
    }
  }
}

cache <- new.env(parent = emptyenv())
for (run in RUNS) {
  which <- run[1]; grp <- run[2]; g <- run[3]
  csv_out <- file.path(OUT, sprintf("KO_%s_%s_%s_diffRegulation.csv", which, grp, g))
  if (file.exists(csv_out)) { cat("跳过(已存在):", basename(csv_out), "\n"); next }

  if (is.null(cache[[which]])) cache[[which]] <- load_counts(which)
  dat <- cache[[which]]
  cells_keep <- rownames(dat$meta)[dat$meta$group == grp]
  cg <- dat$cnt[, cells_keep]
  gene_mean <- rowMeans(cg)
  top_genes <- names(sort(gene_mean, decreasing = TRUE))[1:1500]
  cg <- cg[union(top_genes, HUB), ]
  pct <- if (g %in% rownames(cg)) mean(cg[g, ] > 0) else 0
  cat(sprintf("\n===== %s %s 组敲除 %s（组内检出率 %.1f%%, nc_nCells=500）=====\n",
              which, grp, g, pct * 100))
  if (pct < 0.005) {
    cat("  组内表达过低(<0.5%), 如实记录为不可敲\n")
    write.csv(data.frame(gene = g, note = "not_expressed_in_group"),
              file.path(OUT, sprintf("KO_%s_%s_%s_NOTKNOCKABLE.csv", which, grp, g)),
              row.names = FALSE)
    next
  }
  set.seed(42)
  keep <- sample(seq_len(ncol(cg)), min(3500, ncol(cg)))
  cgm <- cg[, keep]
  cgm <- cgm[rowSums(cgm) > 0, ]
  if (!(g %in% rownames(cgm))) {
    cat("  过滤后不在矩阵, 记录不可敲\n")
    write.csv(data.frame(gene = g, note = "removed_by_QC"),
              file.path(OUT, sprintf("KO_%s_%s_%s_NOTKNOCKABLE.csv", which, grp, g)),
              row.names = FALSE)
    next
  }
  t0 <- Sys.time()
  res <- tryCatch(
    scTenifoldKnk(countMatrix = cgm, gKO = g, qc_minPCT = 0.005,
                  nc_nNet = 2, nc_nCells = 500, nCores = 1),
    error = function(e) { cat("  敲除失败:", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) next
  cat("耗时:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "分钟\n")
  dr <- res$diffRegulation
  dr <- dr[order(dr$distance, decreasing = TRUE), ]
  write.csv(dr, csv_out, row.names = FALSE)
  cat(sprintf("受影响基因: %d | Top5: %s\n", nrow(dr), paste(head(dr$gene, 5), collapse = ", ")))
}
cat("\n===== 补跑完成 =====\n")
