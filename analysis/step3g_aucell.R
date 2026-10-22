# ============================================================
# Step 3g —— AUCell 分析：6 个 hub 基因的细胞内调控活性
# 方法: 每个 hub 基因取 top50 正共表达基因构成基因集,
#       AUCell 对每个细胞打分 = 该 hub 基因共表达程序的活性;
#       另加 6 基因 HubSignature 评分。
# 适用: GSE167363 (Sepsis CD14+ Mono) 与 GSE151263 (ARDS PBMC)
# 输出: PDF + PNG(300) + TIFF(500)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(AUCell)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(reshape2)
  library(ragg)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"

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
TOPN <- 50   # 共表达基因集大小

run_aucell <- function(rds_path, out_dir, dataset_name,
                       input_dir = NULL) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  obj <- tryCatch(readRDS(rds_path), error = function(e) NULL)
  umap_df_ext <- NULL
  if (is.null(obj)) {
    # rds 损坏/缺失 -> 从 counts.mtx + meta.csv 重建
    cat(sprintf("[%s] rds 无法读取, 从 %s 重建对象\n", dataset_name, input_dir))
    cnt <- readMM(file.path(input_dir, "counts.mtx"))
    genes <- readLines(file.path(input_dir, "genes.tsv"))
    meta0 <- read.csv(file.path(input_dir, "meta.csv"), row.names = 1)
    rownames(cnt) <- genes
    colnames(cnt) <- rownames(meta0)
    obj <- CreateSeuratObject(counts = cnt, meta.data = meta0)
    obj <- NormalizeData(obj)
    if (file.exists(file.path(input_dir, "umap.csv"))) {
      umap_df_ext <- read.csv(file.path(input_dir, "umap.csv"), row.names = 1)
      colnames(umap_df_ext) <- c("UMAP_1", "UMAP_2")
    }
  }
  if (length(Layers(obj[["RNA"]])) > 1) obj <- JoinLayers(obj)
  set.seed(42)

  norm <- as.matrix(GetAssayData(obj, layer = "data"))
  meta <- obj[[]]
  hub_avail <- intersect(HUB, rownames(norm))
  cat(sprintf("\n\n########## %s ##########\n", dataset_name))
  cat("细胞数:", ncol(obj), "| 可分析 hub 基因:", paste(hub_avail, collapse = ", "), "\n")

  # ---------- 构建基因集 ----------
  # 过滤低表达基因再算相关（加速 + 降噪）
  keep_genes <- rownames(norm)[Matrix::rowMeans(norm > 0) > 0.03]
  keep_genes <- setdiff(keep_genes, c(paste0("RPL", 1:60), paste0("RPS", 1:60), "MT-"))
  keep_genes <- union(keep_genes, hub_avail)
  sub <- norm[intersect(keep_genes, rownames(norm)), ]

  gene_sets <- list()
  for (g in hub_avail) {
    if (!(g %in% rownames(sub))) next
    # 下采样细胞加速相关计算
    set.seed(42)
    cells_cor <- sample(colnames(sub), min(4000, ncol(sub)))
    cors <- apply(sub[-which(rownames(sub) == g), cells_cor, drop = FALSE], 1,
                  function(x) suppressWarnings(cor(x, sub[g, cells_cor], method = "spearman")))
    top <- names(sort(cors, decreasing = TRUE))[seq_len(min(TOPN, sum(!is.na(cors))))]
    gene_sets[[paste0(g, "_program")]] <- union(top, g)
  }
  gene_sets[["HubSignature"]] <- hub_avail
  cat("基因集大小:", paste(sapply(gene_sets, length), collapse = ", "), "\n")

  # ---------- AUCell 打分 ----------
  rankings <- AUCell_buildRankings(norm, plotStats = FALSE, nCores = 1)
  auc <- AUCell_calcAUC(gene_sets, rankings, nCores = 1)
  auc_mat <- as.matrix(getAUC(auc))  # 基因集 x 细胞

  # 写出评分矩阵（细胞 x 基因集）
  scores <- t(auc_mat)
  write.csv(data.frame(cell = rownames(scores), meta[rownames(scores), c("celltype", "group")], scores),
            file.path(out_dir, paste0(dataset_name, "_AUCell_scores.csv")), row.names = FALSE)

  # 保存基因集
  gs_df <- data.frame(gene_set = rep(names(gene_sets), sapply(gene_sets, length)),
                      gene = unlist(gene_sets))
  write.csv(gs_df, file.path(out_dir, paste0(dataset_name, "_AUCell_genesets.csv")), row.names = FALSE)

  # ---------- 可视化 ----------
  sc <- data.frame(cell = rownames(scores),
                   celltype = meta[rownames(scores), "celltype"],
                   group = meta[rownames(scores), "group"],
                   scores, check.names = FALSE)
  melt_sc <- melt(sc, id.vars = c("cell", "celltype", "group"),
                  variable.name = "GeneSet", value.name = "AUC")
  melt_sc$GeneSet <- gsub("_program$", "", as.character(melt_sc$GeneSet))

  # 1) 小提琴图: 各基因集活性 x 细胞类型（分组拆分）
  p_v <- ggplot(melt_sc, aes(x = celltype, y = AUC, fill = group)) +
    geom_violin(scale = "width", position = position_dodge(0.9), linewidth = 0.3) +
    facet_wrap(~GeneSet, ncol = 4, scales = "free_y") +
    scale_fill_manual(values = c(HC = "#8FA2D6", Sepsis = "#D99AAE",
                                 ARDS = "#D99AAE", Control = "#8FA2D6")) +
    ggtitle(paste(dataset_name, "— Hub Gene Program Activity (AUCell)")) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig3(p_v, paste0(dataset_name, "_AUCell_violin"), w = 14, h = 8, dir = out_dir)

  # 2) 热图: 基因集 x 细胞类型（分组）平均 AUC
  agg <- aggregate(melt_sc[, "AUC", drop = FALSE],
                   by = list(GeneSet = melt_sc$GeneSet,
                             celltype = melt_sc$celltype,
                             group = melt_sc$group), mean)
  agg$lab <- paste(agg$celltype, agg$group, sep = "\n")
  hm <- ggplot(agg, aes(x = lab, y = GeneSet, fill = AUC)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradientn(colors = c("#CFE3E6", "#8FA2D6", "#6A5C95", "#5B608C")) +
    ggtitle(paste(dataset_name, "— Mean AUCell Activity")) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = NULL)
  save_fig3(hm, paste0(dataset_name, "_AUCell_heatmap"), w = 9, h = 5, dir = out_dir)

  # 3) 统计: 每基因集 x 细胞类型, HC/Control vs 疾病 Wilcoxon
  stat_res <- do.call(rbind, lapply(split(melt_sc, paste(melt_sc$GeneSet, melt_sc$celltype)), function(d) {
    groups <- unique(d$group)
    if (length(groups) != 2) return(NULL)
    g1 <- d$AUC[d$group == groups[1]]; g2 <- d$AUC[d$group == groups[2]]
    if (length(g1) < 10 || length(g2) < 10) return(NULL)
    wt <- suppressWarnings(wilcox.test(g1, g2))
    data.frame(GeneSet = d$GeneSet[1], celltype = d$celltype[1],
               ref = groups[1], case = groups[2],
               mean_ref = mean(g1), mean_case = mean(g2),
               diff = mean(g2) - mean(g1), W = wt$statistic, p = wt$p.value)
  }))
  if (!is.null(stat_res)) {
    stat_res$padj <- p.adjust(stat_res$p, method = "BH")
    write.csv(stat_res, file.path(out_dir, paste0(dataset_name, "_AUCell_stats.csv")), row.names = FALSE)
    cat("\nAUCell 组间差异（BH<0.05 的条目）:\n")
    print(subset(stat_res, padj < 0.05))
  }

  # 4) UMAP 着色
  tryCatch({
    if (!is.null(umap_df_ext)) {
      umap_df <- umap_df_ext
    } else {
      umap_df <- data.frame(Embeddings(obj, "umap"))
      colnames(umap_df) <- c("UMAP_1", "UMAP_2")
    }
    umap_df$cell <- rownames(umap_df)
    plot_df <- merge(umap_df, sc[, c("cell", colnames(scores))], by = "cell")
    plots <- lapply(colnames(scores), function(gs) {
      ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[gs]])) +
        geom_point(size = 0.3, stroke = 0) +
        scale_color_gradientn(colors = c("#CFE3E6", "#8FA2D6", "#6A5C95", "#5B608C")) +
        ggtitle(gsub("_program$", "", gs)) +
        theme_pub() + NoLegend() +
        theme(plot.title = element_text(size = 9))
    })
    save_fig3(wrap_plots(plots, ncol = 4), paste0(dataset_name, "_AUCell_umap"),
              w = 13, h = 8, dir = out_dir)
  }, error = function(e) cat("UMAP 图失败:", conditionMessage(e), "\n"))

  cat(sprintf("\n===== %s AUCell 完成 =====\n", dataset_name))
}

# ---------- GSE167363 ----------
run_aucell(
  rds_path = file.path(BASE, "08_singlecell/GSE167363/GSE167363_seurat_annotated.rds"),
  out_dir  = file.path(BASE, "08_singlecell/GSE167363/aucell"),
  dataset_name = "GSE167363"
)

# ---------- GSE151263（补缺口; rds 损坏 -> cellchat_input 重建）----------
run_aucell(
  rds_path = file.path(BASE, "08_singlecell/GSE151263_seurat_annotated.rds"),
  out_dir  = file.path(BASE, "08_singlecell/aucell"),
  dataset_name = "GSE151263",
  input_dir = file.path(BASE, "08_singlecell/cellchat_input")
)
