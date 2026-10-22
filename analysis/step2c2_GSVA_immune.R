# Step2c2: GSVA 免疫通路富集评分 + pheatmap 热图
# 数据: MSigDB c5.all.v7.5.1 GO 基因集 -> 筛选免疫相关 GO BP
suppressPackageStartupMessages({
  library(GSVA)
  library(GSEABase)
  library(pheatmap)
  library(limma)
})

setwd("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis")
OUT <- "07_hub_immune/gsva"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ============ 1. 读 C5 GMT，筛选免疫相关 GO BP ============
cat("[1] 读取 C5 基因集 ...\n")
gmt_file <- "00_rawdata/msigdb/c5.all.v7.5.1.symbols.gmt"
c5 <- getGmt(gmt_file, geneIdType = SymbolIdentifier())
c5_names <- names(c5)
cat("C5 总基因集:", length(c5_names), "\n")

# 免疫关键词（精确，避免 CYTOKINE 误匹配 CYTOKINESIS）
kw <- c("IMMUNE", "IMMUNITY", "IMMUNOLOGICAL", "INFLAMMAT", "INFLAMMASOME",
        "T_CELL", "T_LYMPHOCYTE", "B_CELL", "B_LYMPHOCYTE", "LYMPHOCYTE",
        "MACROPHAGE", "MONOCYTE", "NEUTROPHIL", "MYELOID", "PHAGOCYT",
        "CYTOKINE_PRODUCTION", "CYTOKINE_MEDIATED", "CYTOKINE_SIGNALING",
        "CYTOKINE_SECRETION", "CYTOKINE_RECEPTOR", "INTERFERON", "INTERLEUKIN",
        "CHEMOKINE", "ANTIGEN", "LEUKOCYTE", "COMPLEMENT", "TOLL_LIKE",
        "NF_KAPPAB", "ADAPTIVE_IMMUNE", "INNATE_IMMUNE", "IMMUNOGLOBULIN",
        "MAST_CELL", "DENDRITIC_CELL", "NATURAL_KILLER", "NK_T_CELL",
        "REGULATORY_T_CELL", "IMMUNE_RESPONSE", "IMMUNE_SYSTEM", "IMMUNE_EFFECTOR",
        "IMMUNE_REGULAT", "HUMORAL_IMMUNE", "CELL_MEDIATED_IMMUNITY",
        "DEFENSE_RESPONSE", "T_CELL_ACTIVATION", "T_CELL_DIFFERENTIATION",
        "T_CELL_PROLIFERATION", "T_CELL_RECEPTOR", "B_CELL_ACTIVATION",
        "B_CELL_RECEPTOR", "ANTIBODY", "T_CELL_MEDIATED", "LEUKOCYTE_ACTIVATION",
        "LEUKOCYTE_MIGRATION", "LEUKOCYTE_DIFFERENTIATION")
# 误匹配排除词
exclude <- c("CYTOKINESIS", "CILIUM", "FLAGELLUM", "NEURAL_CREST", "CELL_MOTILITY",
             "CELL_MIGRATION", "CELL_CYCLE", "CELL_DIVISION")

is_immune <- function(name) {
  if (!grepl("^GOBP_", name)) return(FALSE)
  any(sapply(kw, function(k) grepl(k, name))) && !any(sapply(exclude, function(e) grepl(e, name)))
}
immune_names <- c5_names[sapply(c5_names, is_immune)]
cat("免疫相关 GO 基因集:", length(immune_names), "\n")

# 基因数筛选（10-500，排除过窄/过宽）
immune_sizes <- sapply(immune_names, function(n) length(geneIds(c5[[n]])))
immune_names <- immune_names[immune_sizes >= 10 & immune_sizes <= 500]
cat("基因数筛选(10-500)后:", length(immune_names), "\n")

# ============ 2. GSVA 评分（两数据集） ============
run_gsva <- function(ds, disease) {
  cat(sprintf("\n[2] GSVA 评分 %s (%s) ...\n", ds, disease))
  expr <- readRDS(paste0("01_DEG/", ds, "_gene_expr.rds"))
  grp <- readRDS(paste0("01_DEG/", ds, "_grp.rds"))
  # 基因集 -> 基因列表
  gs <- lapply(immune_names, function(n) geneIds(c5[[n]]))
  names(gs) <- immune_names
  # GSVA（log2 表达矩阵，method="gsva"）
  gsva_par <- gsvaParam(as.matrix(expr), gs, kcdf = "Gaussian", maxDiff = TRUE)
  gsva_res <- gsva(gsva_par)
  # 组间差异（limma）
  group <- factor(grp$group, levels = c("Control", disease))
  design <- model.matrix(~ group)
  fit <- lmFit(gsva_res, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  res$pathway <- rownames(res)
  res$direction <- ifelse(res$logFC > 0, "Up_in_disease", "Down_in_disease")
  list(gsva = gsva_res, group = group, limma = res)
}

r1 <- run_gsva("GSE76293", "ARDS")
r2 <- run_gsva("GSE69528", "Sepsis")

# ============ 3. 保存 GSVA 得分 + limma 结果 ============
write.csv(r1$gsva, file.path(OUT, "GSE76293_GSVA_scores.csv"))
write.csv(r2$gsva, file.path(OUT, "GSE69528_GSVA_scores.csv"))
write.csv(r1$limma, file.path(OUT, "GSE76293_GSVA_limma.csv"), row.names = FALSE)
write.csv(r2$limma, file.path(OUT, "GSE69528_GSVA_limma.csv"), row.names = FALSE)

# ============ 4. pheatmap 热图 ============
plot_heatmap <- function(gsva_res, group, limma_res, ds, disease, n_top = 40) {
  # 选 top 显著通路（|logFC| 排序）
  top_pathways <- head(limma_res$pathway[order(-abs(limma_res$logFC))], n_top)
  mat <- gsva_res[top_pathways, , drop = FALSE]
  # 行名简化（去 GOBP_ 前缀）
  rownames(mat) <- sub("^GOBP_", "", rownames(mat))
  # 样本按组排序
  ord <- order(group)
  mat <- mat[, ord]
  ann <- data.frame(Group = group[ord], row.names = colnames(mat))
  ann_colors <- list(Group = c(Control = "#4DBBD5", ARDS = "#E64B35", Sepsis = "#E64B35"))
  # 配色（红=上调，蓝=下调）
  palette <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  p <- pheatmap(mat, annotation_col = ann, annotation_colors = ann_colors,
                color = palette, scale = "row",
                cluster_cols = FALSE, show_colnames = FALSE,
                fontsize_row = 6, fontsize_col = 8,
                main = sprintf("%s: Immune GO-BP GSVA (top %d)", disease, n_top),
                silent = TRUE)
  ggsave_pheatmap <- function(p, path, w, h) {
    png(path, width = w, height = h, units = "in", res = 300)
    grid::grid.newpage(); grid::grid.draw(p$gtable)
    dev.off()
  }
  ggsave_pheatmap(p, file.path(OUT, sprintf("%s_GSVA_heatmap.png", ds)), 12, 10)
  cat(sprintf("  %s 热图已保存（%d 通路）\n", ds, n_top))
  invisible(p)
}

cat("\n[4] 生成 pheatmap 热图 ...\n")
plot_heatmap(r1$gsva, r1$group, r1$limma, "GSE76293", "ARDS", 40)
plot_heatmap(r2$gsva, r2$group, r2$limma, "GSE69528", "Sepsis", 40)

# ============ 5. 两数据集共享的显著免疫通路 ============
cat("\n[5] 共享显著免疫通路 ...\n")
sig1 <- r1$limma$pathway[r1$limma$adj.P.Val < 0.05]
sig2 <- r2$limma$pathway[r2$limma$adj.P.Val < 0.05]
shared <- intersect(sig1, sig2)
cat("ARDS 显著通路:", length(sig1), "| Sepsis 显著通路:", length(sig2), "| 共享:", length(shared), "\n")
if (length(shared) > 0) {
  writeLines(shared, file.path(OUT, "shared_significant_immune_pathways.txt"))
  cat("共享通路:", paste(sub("^GOBP_", "", head(shared, 20)), collapse = "; "), "\n")
}

cat("\n===== GSVA 免疫通路评分完成 =====\n")
