# ============================================================
# Step 3o —— Monocle3 拟时序分析（GSE151263 + GSE167363）
# 每数据集: preprocess(50 PCs) -> UMAP -> leiden 聚类 -> learn_graph
#           -> order_cells(root=经典单核) -> 拟时序/基因动态/graph_test
# 输出: monocle3/ 三格式图 + pseudotime/Moran's I CSV
# ============================================================
suppressPackageStartupMessages({
  library(monocle3); library(ggplot2); library(Matrix); library(patchwork); library(ragg)
})

# ---- 热修复: monocle3 1.3.1 my.moran.test 在 R 4.6 下全基因 FAIL ----
# 根因: NAOK <- deparse(substitute(na.action)) == "na.pass"，
#       R 4.6 的 deparse 把 stats::na.fail 拆成 c("stats","na.fail") 长度 2，
#       NAOK 变为长度 2 的 logical，spdep::moran 内部 if(NAOK) 报错
#       ('length = 2' in coercion to 'logical(1)')，graph_test 全部基因 FAIL。
# 修复: na.action 固定为 na.fail，NAOK 恒为 FALSE，直接写死标量。
my.moran.test.fixed <- function(x, listw, wc, alternative = "greater", randomisation = TRUE) {
  zero.policy <- TRUE; adjust.n <- TRUE; na.action <- stats::na.fail; drop.EI2 <- FALSE
  NAOK <- FALSE
  x <- na.action(x)
  na.act <- attr(x, "na.action")
  if (!is.null(na.act)) {
    if (length(listw$neighbours) < 1) warning("bad loop: length(listw$neighbours) < 1")
    subset <- !(seq_along(listw$neighbours) %in% na.act)
    listw <- subset(listw, subset, zero.policy = zero.policy)
  }
  n <- length(listw$neighbours)
  if (n != length(x)) stop("objects of different length")
  S02 <- wc$S0 * wc$S0
  res <- spdep::moran(x, listw, wc$n, wc$S0, zero.policy = zero.policy, NAOK = NAOK)
  I <- res$I; K <- res$K; EI <- (-1)/wc$n1
  if (randomisation) {
    VI <- wc$n * (wc$S1 * (wc$nn - 3 * wc$n + 3) - wc$n * wc$S2 + 3 * S02)
    tmp <- K * (wc$S1 * (wc$nn - wc$n) - 2 * wc$n * wc$S2 + 6 * S02)
    if (tmp > VI) warning("Kurtosis overflow")
    VI <- (VI - tmp)/(wc$n1 * wc$n2 * wc$n3 * S02)
    if (!drop.EI2) VI <- (VI - EI^2)
    if (VI < 0) warning("Negative variance")
  } else {
    VI <- (wc$nn * wc$S1 - wc$n * wc$S2 + 3 * S02)/(S02 * (wc$nn - 1))
    if (!drop.EI2) VI <- (VI - EI^2)
    if (VI < 0) warning("Negative variance")
  }
  ZI <- (I - EI)/sqrt(VI)
  statistic <- ZI
  names(statistic) <- "Moran I statistic standard deviate"
  if (alternative == "two.sided") PrI <- 2 * stats::pnorm(abs(ZI), lower.tail = FALSE)
  else if (alternative == "greater") PrI <- stats::pnorm(ZI, lower.tail = FALSE)
  else PrI <- stats::pnorm(ZI)
  if (!is.finite(PrI) || PrI < 0 || PrI > 1) warning("Out-of-range p-value")
  vec <- c(I, EI, VI)
  names(vec) <- c("Moran I statistic", "Expectation", "Variance")
  res <- list(statistic = statistic, p.value = PrI, estimate = vec)
  if (!is.null(na.act)) attr(res, "na.action") <- na.act
  class(res) <- "htest"
  res
}
environment(my.moran.test.fixed) <- asNamespace("monocle3")
assignInNamespace("my.moran.test", my.moran.test.fixed, ns = "monocle3")
cat("my.moran.test 已热修复 (R 4.6 deparse 兼容)\\n")

source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/common.R")
set.seed(42)
BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
SC <- file.path(BASE, "08_singlecell")
OUT <- file.path(SC, "monocle3")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
HUB <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")
PAL7 <- c("#6A5C95", "#D99AAE", "#E7A889", "#A9CBB8", "#CFE3E6", "#8FA2D6", "#5B608C")

save_fig3 <- function(p, name, w = 7, h = 6, line_art = FALSE) {
  suppressWarnings(ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf"))
  suppressWarnings(ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300))
  agg_tiff(file.path(OUT, paste0(name, ".tiff")), width = w, height = h,
           units = "in", res = if (line_art) 1000 else 500, compression = "lzw")
  print(p); invisible(dev.off())
  invisible(p)
}

load_dataset <- function(which) {
  if (which == "GSE151263") cd <- file.path(SC, "cellchat_input") else
    cd <- file.path(SC, "GSE167363/input")
  cnt <- readMM(file.path(cd, "counts.mtx"))
  genes <- sub("\\t.*$", "", readLines(file.path(cd, "genes.tsv")))
  rownames(cnt) <- genes
  meta <- read.csv(file.path(cd, "meta.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  rownames(meta) <- meta[[1]]
  colnames(cnt) <- rownames(meta)
  list(cnt = as(cnt, "dgCMatrix"), meta = meta)
}

get_earliest_principal_node <- function(cds, celltype) {
  cell_ids <- which(colData(cds)$celltype == celltype)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))]
  root_pr_nodes
}

# 每个 partition 各选一个根节点：根细胞类型细胞数 >=5 时用其最近主图节点众数，
# 否则用该 partition 全部细胞（保证所有连通分量都有根，拟时序全覆盖）
get_root_nodes_per_partition <- function(cds, celltype) {
  part <- partitions(cds)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  vnames <- igraph::V(principal_graph(cds)[["UMAP"]])$name
  roots <- character(0)
  for (p in unique(as.character(part))) {
    cells_p <- names(part)[as.character(part) == p]
    ct_cells <- cells_p[colData(cds)[cells_p, ]$celltype == celltype]
    pool <- if (length(ct_cells) >= 5) ct_cells else cells_p
    tab <- table(closest_vertex[pool, ])
    roots <- c(roots, vnames[as.numeric(names(which.max(tab)))])
  }
  cat("  partitions:", length(unique(as.character(part))), "| roots:", length(roots), "\n")
  roots
}

run_monocle3 <- function(which, root_type, group_cols) {
  cat("\n=====", which, "=====\n")
  dat <- load_dataset(which)
  cds <- new_cell_data_set(dat$cnt,
                           cell_metadata = DataFrame(dat$meta),
                           gene_metadata = DataFrame(gene_short_name = rownames(dat$cnt),
                                                     row.names = rownames(dat$cnt)))
  cat("cells:", ncol(cds), "genes:", nrow(cds), "\n")
  t0 <- Sys.time()
  cds <- preprocess_cds(cds, num_dim = 50)
  cds <- reduce_dimension(cds, preprocess_method = "PCA", reduction_method = "UMAP",
                          umap.min_dist = 0.3, umap.n_neighbors = 30L)
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds, use_partition = TRUE)
  cds <- order_cells(cds, root_pr_nodes = get_root_nodes_per_partition(cds, root_type))
  cat("graph+order 耗时:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
  # R 4.6 序列化深嵌套 S4 对象会触发 C stack 溢出（沙箱栈限 8MB 不可调），保存失败不阻断流程
  tryCatch({ saveRDS(cds, file.path(OUT, paste0("monocle3_cds_", which, ".rds"))) },
           error = function(e) cat("cds rds 保存失败(不影响图件输出):", conditionMessage(e), "\n"))

  grp_col <- group_cols[["group"]]; ct_col <- group_cols[["celltype"]]
  ct_uniq <- sort(unique(colData(cds)[[ct_col]]))
  ct_cols <- setNames(PAL7[(seq_along(ct_uniq) - 1) %% 7 + 1], ct_uniq)
  g_cols <- group_cols[["gcols"]]

  # ---- 图 1: 轨迹 by celltype ----
  p1 <- plot_cells(cds, color_cells_by = ct_col, label_groups_by_cluster = FALSE,
                   label_leaves = FALSE, label_branch_points = TRUE,
                   graph_label_size = 3, cell_size = 0.6) +
    scale_color_manual(values = ct_cols, name = NULL) +
    ggtitle(paste(which, ": Monocle3 trajectory (", root_type, " as root)", sep = "")) +
    theme_pub(10) + theme(legend.position = "right")
  save_fig3(p1, paste0(which, "_monocle3_trajectory_celltype"), w = 8.2, h = 6.2)

  # ---- 图 2: 轨迹 by group ----
  p2 <- plot_cells(cds, color_cells_by = grp_col, label_groups_by_cluster = FALSE,
                   label_leaves = FALSE, label_branch_points = TRUE,
                   graph_label_size = 3, cell_size = 0.6) +
    scale_color_manual(values = g_cols, name = NULL) +
    ggtitle(paste(which, ": trajectory colored by group")) +
    theme_pub(10) + theme(legend.position = "right")
  save_fig3(p2, paste0(which, "_monocle3_trajectory_group"), w = 8.2, h = 6.2)

  # ---- 图 3: pseudotime ----
  p3 <- plot_cells(cds, color_cells_by = "pseudotime", label_branch_points = TRUE,
                   graph_label_size = 3, cell_size = 0.6) +
    ggtitle(paste(which, ": Monocle3 pseudotime")) + theme_pub(10)
  save_fig3(p3, paste0(which, "_monocle3_pseudotime"), w = 8, h = 6.2)

  # ---- 图 4: hub 基因沿轨迹表达 ----
  hg <- intersect(HUB, rownames(cds))
  p4 <- plot_cells(cds, genes = hg, label_cell_groups = FALSE, show_trajectory_graph = TRUE,
                   cell_size = 0.6, graph_label_size = 2.5) +
    ggtitle(paste(which, ": hub gene expression along trajectory")) + theme_pub(10)
  save_fig3(p4, paste0(which, "_monocle3_hub_expression"), w = 12, h = 7.5)

  # ---- 图 5: hub 基因随拟时序动态 ----
  pt <- pseudotime(cds)
  keep <- is.finite(pt) & pt > 0
  small <- cds[, keep]
  p5 <- plot_genes_in_pseudotime(small[hg, ], min_expr = 0.5, color_cells_by = grp_col) +
    scale_color_manual(values = g_cols, name = NULL) +
    ggtitle(paste(which, ": hub gene dynamics along pseudotime")) + theme_pub(10)
  save_fig3(p5, paste0(which, "_monocle3_hub_pseudotime_dynamics"), w = 12, h = 7.5)

  # ---- 图 6: 分组拟时序分布 ----
  df <- data.frame(pseudotime = as.numeric(pt),
                   group = as.character(colData(cds)[[grp_col]]),
                   celltype = as.character(colData(cds)[[ct_col]]))
  df <- df[is.finite(df$pseudotime), ]
  p6 <- ggplot(df, aes(x = group, y = pseudotime, fill = group)) +
    geom_violin(alpha = 0.7, color = "grey25", linewidth = 0.3, trim = FALSE) +
    geom_boxplot(width = 0.14, outlier.size = 0.4, fill = "white",
                 color = "grey25", linewidth = 0.3) +
    scale_fill_manual(values = g_cols, name = NULL) +
    ggtitle(paste(which, ": pseudotime distribution by group")) +
    theme_pub(10) + theme(legend.position = "none")
  save_fig3(p6, paste0(which, "_monocle3_pseudotime_group"), w = 5.6, h = 5.2)
  if (length(unique(df$group)) == 2) {
    wt <- wilcox.test(pseudotime ~ group, data = df)
    cat(which, "组间拟时序 wilcox p =", signif(wt$p.value, 4), "\n")
    write.csv(data.frame(dataset = which, statistic = unname(wt$statistic),
                         p.value = wt$p.value), row.names = FALSE,
              file.path(OUT, paste0("monocle3_pseudotime_wilcox_", which, ".csv")))
  }

  # ---- CSV: pseudotime ----
  out <- data.frame(cell = colnames(cds), pseudotime = as.numeric(pt),
                    celltype = as.character(colData(cds)[[ct_col]]),
                    group = as.character(colData(cds)[[grp_col]]),
                    cluster = tryCatch(as.character(clusters(cds)),
                                       error = function(e) rep(NA, ncol(cds))))
  write.csv(out, file.path(OUT, paste0("monocle3_pseudotime_", which, ".csv")), row.names = FALSE)

  # ---- graph_test (Moran's I)，已有 CSV 且非全 FAIL 则复用 ----
  gt_csv <- file.path(OUT, paste0("monocle3_graph_test_", which, ".csv"))
  reuse <- FALSE
  if (file.exists(gt_csv)) {
    pgd <- read.csv(gt_csv, stringsAsFactors = FALSE)
    fail_frac <- if ("status" %in% names(pgd)) mean(pgd$status == "FAIL", na.rm = TRUE) else 1
    if (fail_frac > 0.5) {
      cat(which, "既有 graph_test 全 FAIL（旧 bug），删除并重算\n")
      file.remove(gt_csv)
    } else {
      reuse <- TRUE
      cat(which, "graph_test 复用既有 CSV\n")
    }
  }
  if (!reuse) {
    t1 <- Sys.time()
    pg <- graph_test(cds, neighbor_graph = "principal_graph", cores = 1)
    pgd <- as.data.frame(pg)
    pgd$gene <- rownames(pgd)
    write.csv(pgd, gt_csv, row.names = FALSE)
    cat(which, "graph_test 耗时:", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), "min\n")
  }
  cat(which, "graph_test:", nrow(pgd), "genes |",
      sum(pgd$q_value < 0.05, na.rm = TRUE), "significant (q<0.05)\n")
  # hub 基因 Moran's I 汇总
  hub_tab <- pgd[pgd$gene %in% HUB, c("gene", "morans_I", "morans_test_statistic", "q_value")]
  hub_tab <- hub_tab[order(-hub_tab$morans_I), ]
  write.csv(hub_tab, file.path(OUT, paste0("monocle3_hub_moransI_", which, ".csv")), row.names = FALSE)
  print(hub_tab)

  # ---- 图 7: top 拟时序基因 + hub Moran's I 条形图 ----
  topg <- head(pgd[order(pgd$q_value, -pgd$morans_I), "gene"], 12)
  topg <- topg[!is.na(topg)]
  p7a <- plot_cells(cds, genes = topg, label_cell_groups = FALSE, cell_size = 0.55,
                    graph_label_size = 2) +
    ggtitle(paste(which, ": top autogenes (Moran's I)")) + theme_pub(9.5)
  save_fig3(p7a, paste0(which, "_monocle3_top_autogenes"), w = 12, h = 8.5)

  mi <- hub_tab[!is.na(hub_tab$morans_I), ]
  if (nrow(mi) > 0) {
    mi$gene <- factor(mi$gene, levels = mi$gene[order(mi$morans_I)])
    mi$sig <- ifelse(mi$q_value < 0.05, "q<0.05", "ns")
    p7b <- ggplot(mi, aes(x = morans_I, y = gene, fill = sig)) +
      geom_col(width = 0.65, color = "grey25", linewidth = 0.25) +
      scale_fill_manual(values = c("q<0.05" = "#6A5C95", ns = "#CFE3E6"), name = NULL) +
      ggtitle(paste(which, ": hub gene Moran's I along trajectory")) +
      xlab("Moran's I") + ylab(NULL) + theme_pub(10)
    save_fig3(p7b, paste0(which, "_monocle3_hub_moransI_bar"), w = 6, h = 4.4)
  } else {
    cat(which, "hub 基因 Moran's I 全部缺失，跳过条形图\n")
  }
  cat(which, "完成\n")
}

# GSE151263: root=Mono; GSE167363: root=ClassicalMono
run_monocle3("GSE151263", "Mono",
             list(group = "group", celltype = "celltype",
                  gcols = c(ARDS = "#8FA2D6", Sepsis = "#D99AAE")))
run_monocle3("GSE167363", "ClassicalMono",
             list(group = "group", celltype = "celltype",
                  gcols = c(Sepsis = "#D99AAE", HC = "#A9CBB8")))
cat("\n===== Monocle3 全部完成 =====\n")
