# HPA 免疫细胞 RNA 表达: 核心基因的免疫细胞特异性验证
# 数据源: HPA Immune Cell (健康供者原代免疫细胞 RNA-seq, nTPM), 长格式
suppressPackageStartupMessages({library(data.table); library(pheatmap); library(ggplot2)})

setwd("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis")
OUT <- "10_HPA_validation"

# 1) 读长格式数据: Gene, Gene name, Immune cell, TPM, pTPM, nTPM
ic <- fread(file.path(OUT, "rna_immune_cell.tsv"))
cat("长格式行数:", nrow(ic), "| 列:", paste(colnames(ic), collapse=", "), "\n")
cat("免疫细胞类型:", length(unique(ic$`Immune cell`)), "种\n")

# 2) 核心基因(含 VNN2 旁系同源)
target <- c("CDK5RAP2","CD177","GRB10","UPP1","EXOSC4","VNN1","VNN2")
sub <- ic[`Gene name` %in% target]
cat("核心基因命中:", paste(unique(sub$`Gene name`), collapse=", "), "\n")

# 3) 转宽格式(基因 x 细胞类型, 值=nTPM)
mat <- dcast(sub, `Gene name` ~ `Immune cell`, value.var = "nTPM", fun.aggregate = mean)
gene_names <- mat$`Gene name`
mat <- as.matrix(mat[, -1, with=FALSE])
rownames(mat) <- gene_names
storage.mode(mat) <- "double"
mat[is.na(mat)] <- 0
cat("\n核心基因表达矩阵:", nrow(mat), "基因 x", ncol(mat), "细胞类型\n")

# 4) 每个基因的免疫细胞特异性
specificity <- data.frame()
for (g in rownames(mat)) {
  v <- mat[g, ]
  top_idx <- which.max(v)
  specificity <- rbind(specificity, data.frame(
    gene = g,
    top_cell = colnames(mat)[top_idx],
    top_nTPM = round(v[top_idx], 1),
    specificity_ratio = round(v[top_idx] / sum(v), 3)
  ))
}
cat("\n===== 核心基因免疫细胞特异性(最高表达细胞) =====\n")
print(specificity, row.names = FALSE)
write.csv(specificity, file.path(OUT, "HPA_immune_cell_specificity.csv"), row.names = FALSE)

# 5) 热图(z-score 行标准化)
mat_log <- log10(mat + 1)
mat_z <- t(scale(t(mat_log)))
mat_z[is.na(mat_z)] <- 0

png(file.path(OUT, "HPA_immune_cell_heatmap.png"), width = 2200, height = 1100, res = 300)
pheatmap(mat_z, cluster_cols = TRUE, cluster_rows = TRUE,
         color = colorRampPalette(c("#7FB3A4", "white", "#6A5C95"))(100),
         main = "HPA: Hub genes across immune cells (log10 nTPM z-score)",
         fontsize = 11, angle_col = 45)
dev.off()
cat("热图已保存\n")

# 6) 核心基因在各免疫细胞的柱状图(nTPM)
plot_df <- melt(sub, id.vars = "Gene name", measure.vars = "nTPM",
                variable.name = "tmp", value.name = "nTPM")
# 重新构建长表
plot_df <- sub[, .(gene = `Gene name`, cell = `Immune cell`, nTPM = nTPM)]
# 只保留有表达
plot_df <- plot_df[nTPM > 0]
# 排序: 按 CD177 表达
cell_order <- plot_df[gene == "CD177"][order(-nTPM)]$cell
plot_df$cell <- factor(plot_df$cell, levels = cell_order)

p <- ggplot(plot_df, aes(x = cell, y = nTPM, fill = gene)) +
  geom_col(position = position_dodge(0.85), width = 0.75) +
  scale_fill_manual(values = c("CD177"="#6A5C95","VNN1"="#8FA2D6","VNN2"="#D99AAE",
                               "UPP1"="#E7A889","CDK5RAP2"="#A9CBB8","GRB10"="#CFE3E6",
                               "EXOSC4"="#5B608C")) +
  labs(title = "HPA: Hub genes RNA expression in immune cells (healthy donors)",
       x = "Immune cell type", y = "nTPM", fill = "Gene") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUT, "HPA_immune_cell_expression.png"), p, width = 12, height = 6.5, dpi = 300)
cat("柱状图已保存\n")

# 7) 关键结论输出
cat("\n===== 关键生物学发现 =====\n")
for (g in rownames(mat)) {
  v <- mat[g, ]
  top3 <- names(sort(v, decreasing = TRUE))[1:3]
  cat(sprintf("%-10s: 最高表达于 %s, 其次 %s\n", g, top3[1], paste(top3[-1], collapse=", ")))
}

cat("\n===== HPA 免疫细胞表达验证完成 =====\n")
