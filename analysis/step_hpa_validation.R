# HPA 功能验证: CD177 + VNN2 在感染性疾病 vs 健康对照的血浆蛋白差异表达
# 数据源: Human Protein Atlas blood_pea_disease_de (Olink PEA, 泛疾病队列)
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(data.table)})

setwd("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis")
OUT <- "10_HPA_validation"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# 配色(与论文一致)
PAL_UP <- "#6A5C95"   # 上调紫
PAL_DOWN <- "#7FB3A4" # 下调绿

# 1) 读数据
de <- fread(file.path(OUT, "blood_pea_disease_de.tsv"))
cat("总行数:", nrow(de), "\n")

# 2) 提取核心基因 + 旁系同源, 感染性疾病 vs 健康对照
target <- c("CD177", "VNN2")  # VNN2 = VNN1 旁系同源(唯一在Olink面板的VNN家族)
sub <- de[Gene %in% target & Class == "Infection" & Control == "Healthy"]
cat("感染 vs 健康 数据行数:", nrow(sub), "\n")

# 3) 数据整理
sub[, neglogp := -log10(as.numeric(`p-value adjusted`))]
sub[, sig := ifelse(`p-value adjusted` < 0.001, "***",
             ifelse(`p-value adjusted` < 0.01, "**",
             ifelse(`p-value adjusted` < 0.05, "*", "ns")))]
sub[, Disease := gsub("Staphylococcus aureus bacteremia", "S. aureus bacteremia", Disease)]
sub[, Disease := gsub("Streptococcal soft tissue infection", "Streptococcal soft tissue inf.", Disease)]
sub[, Disease := gsub("E.coli pyelonephritis", "E. coli pyelonephritis", Disease)]

# 按 CD177 logFC 排序疾病
disease_order <- sub[Gene == "CD177"][order(-logFC)]$Disease
sub[, Disease := factor(Disease, levels = disease_order)]

write.csv(sub[, .(Gene, Disease, logFC, `p-value adjusted`, neglogp, sig)],
          file.path(OUT, "HPA_CD177_VNN2_infection_vs_healthy.csv"), row.names = FALSE)

# 4) 统计学分析
cat("\n===== CD177 感染性疾病 vs 健康对照 (logFC 统计) =====\n")
cd177 <- sub[Gene == "CD177"]
cat(sprintf("CD177: %d 种感染性疾病, logFC 范围 %.3f ~ %.3f, 中位数 %.3f\n",
    nrow(cd177), min(cd177$logFC), max(cd177$logFC), median(cd177$logFC)))
cat(sprintf("CD177 显著上调(padj<0.05): %d/%d 种疾病\n",
    sum(cd177$`p-value adjusted` < 0.05 & cd177$logFC > 0), nrow(cd177)))

cat("\n===== VNN2 (VNN1旁系同源) =====\n")
vnn2 <- sub[Gene == "VNN2"]
cat(sprintf("VNN2: %d 种感染性疾病, logFC 范围 %.3f ~ %.3f, 中位数 %.3f\n",
    nrow(vnn2), min(vnn2$logFC), max(vnn2$logFC), median(vnn2$logFC)))
cat(sprintf("VNN2 显著上调(padj<0.05): %d/%d 种疾病\n",
    sum(vnn2$`p-value adjusted` < 0.05 & vnn2$logFC > 0), nrow(vnn2)))

# CD177 vs VNN2 相关性
wide <- merge(cd177[, .(Disease, cd177_logFC = logFC)],
              vnn2[, .(Disease, vnn2_logFC = logFC)], by = "Disease")
ct <- cor.test(wide$cd177_logFC, wide$vnn2_logFC, method = "pearson")
cat(sprintf("\nCD177 vs VNN2 logFC Pearson相关: r=%.3f, p=%.3e\n", ct$estimate, ct$p.value))

# 5) 图表 1: 分组柱状图 (logFC + 显著性)
p1 <- ggplot(sub, aes(x = Disease, y = logFC, fill = Gene)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = sig, y = logFC + sign(logFC)*0.15), position = position_dodge(0.8),
            size = 3.2, vjust = ifelse(sub$logFC > 0, 0, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("CD177" = PAL_UP, "VNN2" = PAL_DOWN)) +
  labs(title = "HPA: CD177 & VNN2 plasma protein levels in infections vs healthy",
       subtitle = "Blood PEA (Olink) disease differential expression, control = Healthy",
       x = "Infectious disease", y = "log2 fold change vs healthy",
       fill = "Protein") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
        plot.title = element_text(face = "bold")) +
  coord_flip()
ggsave(file.path(OUT, "HPA_CD177_VNN2_infection_logFC.png"), p1, width = 11, height = 7, dpi = 300)
ggsave(file.path(OUT, "HPA_CD177_VNN2_infection_logFC.pdf"), p1, width = 11, height = 7)
cat("\n柱状图已保存\n")

# 6) 图表 2: CD177 vs VNN2 散点图 + 相关性
p2 <- ggplot(wide, aes(x = cd177_logFC, y = vnn2_logFC)) +
  geom_point(size = 3, color = PAL_UP) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", fill = "grey80", alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  annotate("text", x = max(wide$cd177_logFC)*0.6, y = min(wide$vnn2_logFC)*0.9,
           label = sprintf("r = %.3f\np = %.2e", ct$estimate, ct$p.value), size = 4) +
  labs(title = "CD177 vs VNN2 logFC across infections",
       x = "CD177 log2 fold change (vs healthy)", y = "VNN2 log2 fold change (vs healthy)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(OUT, "HPA_CD177_vs_VNN2_correlation.png"), p2, width = 7, height = 6, dpi = 300)
ggsave(file.path(OUT, "HPA_CD177_vs_VNN2_correlation.pdf"), p2, width = 7, height = 6)
cat("散点图已保存\n")

# 7) 统计总结输出
cat("\n===== HPA 功能验证统计总结 =====\n")
cat("CD177 在脓毒症相关感染中显著上调:\n")
for (d in c("Streptococcal soft tissue inf.", "S. aureus bacteremia", "E. coli pyelonephritis", "Pneumococcal pneumonia")) {
  r <- cd177[Disease == d]
  if (nrow(r) > 0) cat(sprintf("  %-30s logFC=%.2f (%.1f倍) padj=%.1e\n", d, r$logFC, 2^r$logFC, r$`p-value adjusted`))
}
cat("\n===== HPA 验证完成 =====\n")
