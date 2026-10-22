# ============================================================
# Step2c3: GSVA 免疫通路富集分析（修改版）
#   (1) 严格仅保留「疾病 vs 对照」样本：
#       GSE76293 = 12 ARDS + 12 Control（剔除 70 个体外培养 PMN）
#       GSE69528 = 83 Sepsis + 28 Control（剔除 27 个 2 型糖尿病未感染）
#   (2) 免疫相关 GO-BP 通路取 top20（校正 P 值升序，最显著）
#   (3) pheatmap 热图 + 汇总
# ============================================================
suppressPackageStartupMessages({
  library(GSVA)
  library(GSEABase)
  library(pheatmap)
  library(limma)
})

setwd("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis")
OUT <- "07_hub_immune/gsva_top20"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ============ 1. 读 C5 GMT，筛选免疫相关 GO BP ============
cat("[1] 读取 C5 基因集 ...\n")
gmt_file <- "00_rawdata/msigdb/c5.all.v7.5.1.symbols.gmt"
c5 <- getGmt(gmt_file, geneIdType = SymbolIdentifier())
c5_names <- names(c5)
cat("C5 总基因集:", length(c5_names), "\n")

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
exclude <- c("CYTOKINESIS", "CILIUM", "FLAGELLUM", "NEURAL_CREST", "CELL_MOTILITY",
             "CELL_MIGRATION", "CELL_CYCLE", "CELL_DIVISION")

is_immune <- function(name) {
  if (!grepl("^GOBP_", name)) return(FALSE)
  any(sapply(kw, function(k) grepl(k, name))) && !any(sapply(exclude, function(e) grepl(e, name)))
}
immune_names <- c5_names[sapply(c5_names, is_immune)]
immune_sizes <- sapply(immune_names, function(n) length(geneIds(c5[[n]])))
immune_names <- immune_names[immune_sizes >= 10 & immune_sizes <= 500]
cat("免疫相关 GO 基因集(10-500基因):", length(immune_names), "\n")

# ============ 2. 严格样本过滤：仅疾病 vs 对照 ============
# 从原始 groups.csv 重新过滤，剔除 EXCLUDE_cultured / EXCLUDE_t2dm 等混杂样本
load_clean <- function(ds, disease, control) {
  grp_raw <- read.csv(file.path("01_DEG", paste0(ds, "_groups.csv")),
                      stringsAsFactors = FALSE, check.names = FALSE)
  # 记录被剔除的样本
  dropped <- grp_raw$group[!grp_raw$group %in% c(control, disease)]
  # 仅保留疾病 + 对照
  keep <- grp_raw$group %in% c(control, disease)
  grp <- grp_raw[keep, , drop = FALSE]
  # 读取基因级表达矩阵（已探针去重），按样本对齐
  expr <- readRDS(file.path("01_DEG", paste0(ds, "_gene_expr.rds")))
  common <- intersect(grp$sample, colnames(expr))
  expr <- expr[, common, drop = FALSE]
  grp <- grp[match(common, grp$sample), , drop = FALSE]
  grp$group <- factor(grp$group, levels = c(control, disease))
  cat(sprintf("[2] %s: 剔除混杂样本 %d 个 (%s)\n", ds, length(dropped),
              paste(sort(unique(dropped)), collapse = ",")))
  cat(sprintf("    保留: %s=%d, %s=%d, 总样本=%d\n",
              disease, sum(grp$group == disease),
              control, sum(grp$group == control), ncol(expr)))
  list(expr = expr, grp = grp)
}

d1 <- load_clean("GSE76293", "ARDS", "Control")
d2 <- load_clean("GSE69528", "Sepsis", "Control")

# ============ 3. GSVA 评分 + limma 差异 ============
run_gsva <- function(ds, disease, dat) {
  cat(sprintf("\n[3] GSVA 评分 %s (%s) ...\n", ds, disease))
  gs <- lapply(immune_names, function(n) geneIds(c5[[n]]))
  names(gs) <- immune_names
  gsva_par <- gsvaParam(as.matrix(dat$expr), gs, kcdf = "Gaussian", maxDiff = TRUE)
  gsva_res <- gsva(gsva_par)
  design <- model.matrix(~ dat$grp$group)
  fit <- lmFit(gsva_res, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  res$pathway <- rownames(res)
  res$direction <- ifelse(res$logFC > 0, "Up_in_disease", "Down_in_disease")
  list(gsva = gsva_res, group = dat$grp$group, limma = res)
}

r1 <- run_gsva("GSE76293", "ARDS", d1)
r2 <- run_gsva("GSE69528", "Sepsis", d2)

# ============ 4. top20 免疫通路（校正 P 值升序） ============
TOP_N <- 20
select_top <- function(limma_res) {
  ord <- order(limma_res$adj.P.Val, na.last = TRUE)
  head(limma_res[ord, ], TOP_N)
}
top1 <- select_top(r1$limma)
top2 <- select_top(r2$limma)

# 保存完整 limma + top20
write.csv(r1$gsva, file.path(OUT, "GSE76293_GSVA_scores.csv"))
write.csv(r2$gsva, file.path(OUT, "GSE69528_GSVA_scores.csv"))
write.csv(r1$limma, file.path(OUT, "GSE76293_GSVA_limma_all.csv"), row.names = FALSE)
write.csv(r2$limma, file.path(OUT, "GSE69528_GSVA_limma_all.csv"), row.names = FALSE)
write.csv(top1, file.path(OUT, "GSE76293_GSVA_top20.csv"), row.names = FALSE)
write.csv(top2, file.path(OUT, "GSE69528_GSVA_top20.csv"), row.names = FALSE)

# ============ 5. pheatmap 热图（top20） ============
plot_heatmap <- function(gsva_res, group, top_res, ds, disease) {
  top_pathways <- top_res$pathway
  mat <- gsva_res[top_pathways, , drop = FALSE]
  rownames(mat) <- sub("^GOBP_", "", rownames(mat))
  ord <- order(group)
  mat <- mat[, ord]
  grp_ord <- droplevels(group[ord])
  ann <- data.frame(Group = grp_ord, row.names = colnames(mat))
  # 图例仅保留当前数据集实际存在的两个组（ARDS图: Control+ARDS; Sepsis图: Control+Sepsis）
  lv <- levels(grp_ord)
  ann_colors <- list(Group = c(Control = "#8FA2D6", ARDS = "#D99AAE", Sepsis = "#D99AAE")[lv])
  palette <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  p <- pheatmap(mat, annotation_col = ann, annotation_colors = ann_colors,
                color = palette, scale = "row",
                cluster_cols = FALSE, show_colnames = FALSE,
                fontsize_row = 7, fontsize_col = 8,
                main = sprintf("%s: Immune GO-BP GSVA (top %d)", disease, TOP_N),
                silent = TRUE)
  png(file.path(OUT, sprintf("%s_GSVA_heatmap_top20.png", ds)),
      width = 11, height = 9, units = "in", res = 300)
  grid::grid.newpage(); grid::grid.draw(p$gtable)
  dev.off()
  cat(sprintf("  %s 热图已保存（top %d 通路）\n", ds, TOP_N))
  invisible(p)
}

cat("\n[4] 生成 pheatmap 热图 (top20) ...\n")
plot_heatmap(r1$gsva, r1$group, top1, "GSE76293", "ARDS")
plot_heatmap(r2$gsva, r2$group, top2, "GSE69528", "Sepsis")

# ============ 6. 两数据集 top20 的共享通路 ============
cat("\n[5] top20 共享免疫通路 ...\n")
shared <- intersect(top1$pathway, top2$pathway)
cat("ARDS top20 / Sepsis top20 共享:", length(shared), "\n")
if (length(shared) > 0) {
  writeLines(shared, file.path(OUT, "shared_top20_immune_pathways.txt"))
  cat("共享通路:", paste(sub("^GOBP_", "", shared), collapse = "; "), "\n")
}

cat("\n===== GSVA 免疫通路 top20 分析完成 =====\n")
