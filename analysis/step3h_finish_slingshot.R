# ============================================================
# Step 3h-finish —— 补全 slingshot 剩余图（从已保存 sds 对象）
# ============================================================
suppressPackageStartupMessages({
  library(slingshot)
  library(ggplot2)
  library(ragg)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "08_singlecell/GSE167363/slingshot")

save_fig3 <- function(p, name, w = 7, h = 6, dir = ".", line_art = FALSE) {
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf"))
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300))
  tdpi <- if (line_art) 1000 else 500
  agg_tiff(file.path(dir, paste0(name, ".tiff")), width = w, height = h,
           units = "in", res = tdpi, compression = "lzw")
  print(p); invisible(dev.off())
  invisible(p)
}

HUB <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")

sds  <- readRDS(file.path(OUT, "slingshot_sds.rds"))
pca  <- read.csv(file.path(OUT, "mono_pca.csv"), row.names = 1)
meta <- read.csv(file.path(OUT, "mono_meta.csv"), row.names = 1)
expr <- read.csv(file.path(OUT, "mono_gene_expr.csv"), row.names = 1)
pt   <- read.csv(file.path(OUT, "slingshot_pseudotime.csv"), row.names = 1)
pseudo <- apply(pt, 1, min, na.rm = TRUE)

CT_COLS <- c(
  "ClassicalMono"      = "#6A5C95",
  "IntermediateMono"   = "#D99AAE",
  "NonclassicalMono"   = "#E7A889",
  "InflammatoryMono"   = "#5B608C",
  "ISGMono"            = "#8FA2D6",
  "MonoDC"             = "#A9CBB8"
)
CT_COLS <- CT_COLS[intersect(names(CT_COLS), unique(meta$celltype))]

# ---------- PCA 空间轨迹 ----------
lin <- slingCurves(sds)
pca_df <- data.frame(PC1 = pca[, 1], PC2 = pca[, 2],
                     celltype = meta$celltype, pseudotime = pseudo)
base <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = celltype), size = 0.3, stroke = 0) +
  scale_color_manual(values = CT_COLS) + theme_pub()
for (i in seq_along(lin)) {
  lc <- as.data.frame(lin[[i]]$s[lin[[i]]$ord, 1:2, drop = FALSE])
  colnames(lc) <- c("PC1", "PC2")
  base <- base + geom_path(data = lc, aes(x = PC1, y = PC2),
                           color = "#5B608C", linewidth = 1.1, inherit.aes = FALSE)
}
base <- base + ggtitle(sprintf("Slingshot trajectory in PCA space (%d lineages)", length(lin)))
save_fig3(base, "GSE167363_slingshot_trajectory_PCA", w = 7, h = 6, dir = OUT)

# ---------- hub 基因沿伪时间 ----------
hub_avail <- intersect(HUB, colnames(expr))
cor_res <- data.frame()
for (g in hub_avail) {
  ct <- suppressWarnings(cor.test(pseudo, expr[, g], method = "spearman"))
  cor_res <- rbind(cor_res, data.frame(gene = g, rho = ct$estimate, p = ct$p.value))
}
cor_res$fdr <- p.adjust(cor_res$p, method = "BH")
cat("\n=== hub 基因与伪时间 Spearman 相关 ===\n"); print(cor_res)
write.csv(cor_res, file.path(OUT, "slingshot_gene_trajectory_cor.csv"), row.names = FALSE)

pseudo_bin <- cut(pseudo, breaks = 5, labels = c("T1","T2","T3","T4","T5"))
trend <- data.frame()
for (g in hub_avail) {
  m <- tapply(expr[, g], pseudo_bin, mean, na.rm = TRUE)
  trend <- rbind(trend, data.frame(gene = g, t(m)))
}
write.csv(trend, file.path(OUT, "slingshot_gene_trend.csv"), row.names = FALSE)
cat("\n=== hub 基因沿伪时间表达趋势(5 分箱均值) ===\n"); print(trend, row.names = FALSE)

gdf <- data.frame(
  pseudotime = rep(pseudo, length(hub_avail)),
  expr = as.vector(as.matrix(expr[, hub_avail, drop = FALSE])),
  gene = rep(hub_avail, each = nrow(expr))
)
gene_cols <- c("#6A5C95", "#D99AAE", "#E7A889", "#A9CBB8", "#8FA2D6", "#5B608C")
names(gene_cols) <- hub_avail
pg <- ggplot(gdf, aes(x = pseudotime, y = expr)) +
  geom_point(size = 0.2, color = "#CFE3E6", alpha = 0.4) +
  geom_smooth(aes(color = gene), method = "loess", se = FALSE, linewidth = 0.9,
              method.args = list(span = 0.3)) +
  scale_color_manual(values = gene_cols) +
  facet_wrap(~gene, ncol = 3, scales = "free_y") +
  ggtitle("Hub Gene Expression along Pseudotime (GSE167363)") + theme_pub()
save_fig3(pg, "GSE167363_slingshot_gene_trend", w = 12, h = 7, dir = OUT, line_art = TRUE)

cat("\n===== Step3h-finish 完成 =====\n")
