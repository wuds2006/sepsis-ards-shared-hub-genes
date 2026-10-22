# Step3c: slingshot 拟时序 (Mono -> M2巨噬细胞分化轨迹)
suppressPackageStartupMessages({
  library(slingshot)
  library(ggplot2)
})

setwd("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis")
OUT <- "08_singlecell/slingshot"

pca <- read.csv(file.path(OUT, "mono_pca.csv"), row.names = 1)
umap <- read.csv(file.path(OUT, "mono_umap.csv"), row.names = 1)
meta <- read.csv(file.path(OUT, "mono_meta.csv"), row.names = 1)
expr <- read.csv(file.path(OUT, "mono_gene_expr.csv"), row.names = 1)

cat("细胞数:", nrow(pca), "| 亚群:", paste(sort(unique(meta$cluster)), collapse=","), "\n")

# slingshot 用 PCA 坐标, 起始簇=0(classical monocyte, M1高)
rd <- as.matrix(pca[, 1:20])
cl <- meta$cluster
sds <- slingshot(rd, cl, start.clus = "0", end.clus = "7")

# 伪时间
pt <- slingPseudotime(sds)
colnames(pt) <- paste0("Lineage", seq_len(ncol(pt)))
pseudo <- apply(pt, 1, min, na.rm = TRUE)  # 取最小伪时间作为主轨迹

cat("\n=== 轨迹信息 ===\n")
print(slingLineages(sds))

# 保存
saveRDS(sds, file.path(OUT, "slingshot_sds.rds"))
write.csv(pt, file.path(OUT, "slingshot_pseudotime.csv"))
write.csv(data.frame(cell = rownames(pca), pseudotime = pseudo),
          file.path(OUT, "slingshot_main_pseudotime.csv"), row.names = FALSE)

# 可视化轨迹(UMAP)
png(file.path(OUT, "slingshot_trajectory_UMAP.png"), width = 2200, height = 1800, res = 300)
par(mar = c(4,4,2,1))
plot(umap$UMAP1, umap$UMAP2, col = as.integer(as.factor(cl)) + 1, pch = 16, cex = 0.5,
     xlab = "UMAP1", ylab = "UMAP2", main = "Mono differentiation trajectory (slingshot)")
lines(SlingshotDataSet(sds), lwd = 2, col = "black")
dev.off()
cat("轨迹图已保存\n")

# Hub 基因沿伪时间的动态 (相关性 + 分段均值)
hub_genes <- intersect(c("CDK5RAP2","CD177","GRB10","UPP1","EXOSC4","VNN1"), colnames(expr))
cor_res <- data.frame()
for (g in c(hub_genes, "M1_score", "M2_score")) {
  if (!(g %in% colnames(expr))) next
  x <- expr[[g]]
  ct <- cor.test(pseudo, x, method = "spearman", use = "complete.obs")
  cor_res <- rbind(cor_res, data.frame(gene = g, rho = ct$estimate, p = ct$p.value))
}
cor_res$fdr <- p.adjust(cor_res$p, method = "BH")
cat("\n=== Hub基因/M1/M2 与伪时间相关性 ===\n")
print(cor_res)
write.csv(cor_res, file.path(OUT, "slingshot_gene_trajectory_cor.csv"), row.names = FALSE)

# 伪时间分箱的 Hub 基因表达趋势
pseudo_bin <- cut(pseudo, breaks = 5, labels = c("T1","T2","T3","T4","T5"))
trend <- data.frame()
for (g in c(hub_genes, "M1_score", "M2_score")) {
  if (!(g %in% colnames(expr))) next
  m <- tapply(expr[[g]], pseudo_bin, mean, na.rm = TRUE)
  trend <- rbind(trend, data.frame(gene = g, t(m)))
}
write.csv(trend, file.path(OUT, "slingshot_gene_trend.csv"), row.names = FALSE)
cat("\n=== Hub基因沿伪时间表达趋势(5分箱均值) ===\n")
print(trend, row.names = FALSE)

cat("\n===== Step3c slingshot 完成 =====\n")
