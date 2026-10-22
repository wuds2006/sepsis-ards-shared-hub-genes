# ============================================================
# Step 3n-3 —— 模块10a：分组分层 scTenifoldKnk 双虚拟敲除
# GSE151263 (ARDS vs Sepsis) 与 GSE167363 (Sepsis vs HC)：
#   1. 每组内独立构建 GRN 并对可敲 hub 基因虚拟敲除
#   2. 比较敲除后受影响基因（logFC 相关 / 显著基因 Jaccard）
#   3. 受影响基因的 GO BP 富集（clusterProfiler）并做组间通路比较图
# 输出: 08_singlecell/knockout_stratified/
# ============================================================
suppressPackageStartupMessages({
  library(scTenifoldKnk); library(Matrix); library(ggplot2); library(patchwork); library(ragg)
  library(clusterProfiler); library(org.Hs.eg.db)
})
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/common.R")
set.seed(42)

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT <- file.path(BASE, "08_singlecell/knockout_stratified")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
HUB <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")

save_fig3 <- function(p, name, w = 7, h = 6, dir = ".", line_art = FALSE) {
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf"))
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300))
  tdpi <- if (line_art) 1000 else 500
  agg_tiff(file.path(dir, paste0(name, ".tiff")), width = w, height = h,
           units = "in", res = tdpi, compression = "lzw")
  print(p); invisible(dev.off())
  invisible(p)
}

load_counts <- function(which) {
  if (which == "GSE151263") {
    cd <- file.path(BASE, "08_singlecell/cellchat_input")
  } else {
    cd <- file.path(BASE, "08_singlecell/GSE167363/input")
  }
  cnt <- readMM(file.path(cd, "counts.mtx"))
  genes <- readLines(file.path(cd, "genes.tsv"))
  genes <- sub("\t.*$", "", genes)
  rownames(cnt) <- genes
  meta <- read.csv(file.path(cd, "meta.csv"), stringsAsFactors = FALSE)
  cell_col <- if ("cell" %in% names(meta)) "cell" else names(meta)[1]
  rownames(meta) <- meta[[cell_col]]
  # counts.mtx 列顺序 = meta 行顺序（step3b/3e 导出时一致），此处按 meta 顺序补列名
  stopifnot(ncol(cnt) == nrow(meta))
  colnames(cnt) <- rownames(meta)
  list(cnt = cnt, meta = meta)
}

run_dataset <- function(which, groups, gene_candidates) {
  cat("\n##########", which, "##########\n")
  dat <- load_counts(which)
  cnt <- dat$cnt; meta <- dat$meta
  results <- list()   # results[[group]][[gene]] = diffRegulation df
  feas <- list()
  for (grp in groups) {
    cells_keep <- rownames(meta)[meta$group == grp]
    cells_keep <- intersect(cells_keep, colnames(cnt))
    cat(sprintf("\n--- %s: %s 组 %d 细胞 ---\n", which, grp, length(cells_keep)))
    cg <- cnt[, cells_keep]
    # top 1500 高表达基因 + hub
    gene_mean <- rowMeans(cg)
    top_genes <- names(sort(gene_mean, decreasing = TRUE))[1:1500]
    cg <- cg[union(top_genes, HUB), ]
    pct <- sapply(HUB, function(g) if (g %in% rownames(cg)) mean(cg[g, ] > 0) else 0)
    knockable <- gene_candidates
    set.seed(42)
    keep <- sample(seq_len(ncol(cg)), min(2500, ncol(cg)))
    cg <- cg[, keep]
    cg <- cg[rowSums(cg) > 0, ]
    feas[[grp]] <- data.frame(dataset = which, group = grp, gene = HUB,
                              pct_expressing = round(unname(pct), 4))
    results[[grp]] <- list()
    for (g in knockable) {
      if (!(g %in% rownames(cg))) { cat(g, "不在矩阵, 跳过\n"); next }
      csv_out <- file.path(OUT, sprintf("KO_%s_%s_%s_diffRegulation.csv", which, grp, g))
      if (file.exists(csv_out)) {
        cat(sprintf("跳过 %s %s (已存在)\n", grp, g))
        results[[grp]][[g]] <- read.csv(csv_out, stringsAsFactors = FALSE)
        next
      }
      cat(sprintf("\n===== %s %s 组敲除 %s =====\n", which, grp, g))
      t0 <- Sys.time()
      res <- tryCatch(
        scTenifoldKnk(countMatrix = cg, gKO = g, nc_nNet = 2, nc_nCells = 200, nCores = 1),
        error = function(e) { cat("敲除失败:", conditionMessage(e), "\n"); NULL })
      if (is.null(res)) next
      cat("耗时:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "分钟\n")
      dr <- res$diffRegulation
      dr <- dr[order(dr$distance, decreasing = TRUE), ]
      write.csv(dr, csv_out, row.names = FALSE)
      results[[grp]][[g]] <- dr
      cat(sprintf("受影响基因: %d | 显著(FDR<0.05): %d | Top5: %s\n",
                  nrow(dr), sum(dr$p.adj < 0.05), paste(head(dr$gene, 5), collapse = ", ")))
    }
  }
  list(which = which, groups = groups, results = results, feas = do.call(rbind, feas))
}

# ---- GO 富集（受影响基因 top 150 by distance）----
go_enrich <- function(dr, universe) {
  genes <- head(dr$gene, 150)
  eg <- mapIds(org.Hs.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL")
  eg <- eg[!is.na(eg)]
  uni <- mapIds(org.Hs.eg.db, keys = universe, column = "ENTREZID", keytype = "SYMBOL")
  uni <- uni[!is.na(uni)]
  if (length(eg) < 8) return(NULL)
  ego <- tryCatch(
    enrichGO(gene = eg, universe = uni, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
             ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.1, readable = TRUE),
    error = function(e) NULL)
  ego
}

# =============== 运行 ===============
pbmc <- run_dataset("GSE151263", c("ARDS", "Sepsis"), c("UPP1", "EXOSC4"))
cd14 <- run_dataset("GSE167363", c("Sepsis", "HC"), c("UPP1", "EXOSC4", "VNN1", "CDK5RAP2"))

# ---- 敲除可行性表 ----
feas_all <- rbind(pbmc$feas, cd14$feas)
write.csv(feas_all, file.path(OUT, "KO_stratified_feasibility.csv"), row.names = FALSE)

# ---- 组间比较 + 富集 ----
comp_rows <- list()
enrich_list <- list()
datasets <- list(pbmc = pbmc, cd14 = cd14)
for (ds in datasets) {
  which <- ds$which; groups <- ds$groups
  for (g in names(ds$results[[groups[1]]])) {
    a <- ds$results[[groups[1]]][[g]]
    b <- ds$results[[groups[2]]][[g]]
    if (is.null(a) || is.null(b)) next
    common <- intersect(a$gene, b$gene)
    ct <- cor.test(a$FC[match(common, a$gene)], b$FC[match(common, b$gene)], method = "spearman")
    sigA <- a$gene[a$p.adj < 0.05]; sigB <- b$gene[b$p.adj < 0.05]
    jac <- length(intersect(sigA, sigB)) / length(union(sigA, sigB))
    comp_rows[[length(comp_rows) + 1]] <- data.frame(
      dataset = which, gene = g,
      n_sig_case = length(sigA), n_sig_ctrl = length(sigB),
      n_common_genes = length(common),
      spearman_rho = round(unname(ct$estimate), 3),
      spearman_p = signif(ct$p.value, 3),
      jaccard_sig = round(jac, 3))
    # 富集（组间通路差异）
    for (grp in groups) {
      dr <- ds$results[[grp]][[g]]
      ego <- go_enrich(dr, unique(c(dr$gene, HUB)))
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      edf <- as.data.frame(ego)
      edf$group <- grp; edf$dataset <- which; edf$geneKO <- g
      edf$GeneRatioS <- vapply(strsplit(as.character(edf$GeneRatio), "/"),
                               function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
      enrich_list[[length(enrich_list) + 1]] <- edf
      }
    }
  }
}
comp <- do.call(rbind, comp_rows)
write.csv(comp, file.path(OUT, "KO_stratified_group_comparison.csv"), row.names = FALSE)
print(comp)

if (length(enrich_list)) {
  enrich_all <- do.call(rbind, enrich_list)
  write.csv(enrich_all, file.path(OUT, "KO_stratified_GO_BP.csv"), row.names = FALSE)

  # ---- 通路比较图（每个数据集 x 敲除基因：两组 top 8 GO term 分面）----
  for (ds in datasets) {
    which <- ds$which
    sub <- enrich_all[enrich_all$dataset == which, ]
    for (g in unique(sub$geneKO)) {
      sg <- sub[sub$geneKO == g, ]
      if (nrow(sg) == 0) next
      top_terms <- unique(unlist(lapply(split(sg, sg$group), function(d)
        head(d[order(d$p.adjust), "Description"], 8))))
      sg <- sg[sg$Description %in% top_terms, ]
      sg$Description <- factor(sg$Description, levels = rev(top_terms))
      p <- ggplot(sg, aes(x = group, y = Description, size = GeneRatioS, color = p.adjust)) +
        geom_point() +
        scale_color_gradientn(colors = c("#5B608C", "#8FA2D6", "#CFE3E6")) +
        facet_grid(. ~ geneKO, drop = TRUE) +
        ggtitle(paste0(which, ": KO-perturbed gene GO BP enrichment (top 8 per group)")) +
        xlab(NULL) + ylab(NULL) +
        theme_pub(10) + theme(axis.text.y = element_text(size = 7.5),
                              strip.text = element_text(face = "italic", size = 9))
      save_fig3(p, paste0(which, "_KO_stratified_GO_BP_", g), w = 9.5,
                h = max(3.2, 0.42 * length(top_terms) + 1.6), dir = OUT)
    }
  }
}

# ---- top10 受影响基因组间对比图（每基因 2 组并排）----
for (ds in datasets) {
  which <- ds$which; groups <- ds$groups
  for (g in names(ds$results[[groups[1]]])) {
    a <- ds$results[[groups[1]]][[g]]; b <- ds$results[[groups[2]]][[g]]
    if (is.null(a) || is.null(b)) next
    ta <- head(a, 10); ta$group <- groups[1]; ta$grp_f <- factor(groups[1], levels = groups)
    tb <- head(b, 10); tb$group <- groups[2]; tb$grp_f <- factor(groups[2], levels = groups)
    tt <- rbind(ta, tb)
    tt$gene <- factor(tt$gene, levels = rev(unique(tt$gene)))
    tt$dir <- ifelse(tt$distance > 0, "up", "down")
    p <- ggplot(tt, aes(x = grp_f, y = gene)) +
      geom_point(aes(size = abs(distance), fill = dir), shape = 21, color = "grey20", stroke = 0.4) +
      scale_fill_manual(values = c(up = "#D99AAE", down = "#8FA2D6"), name = "Direction") +
      scale_size_continuous(range = c(1.5, 6), name = "|dRegulation|") +
      ggtitle(paste0(which, ": scTenifoldKnk KO ", g, " — top 10 perturbed genes per group")) +
      xlab(NULL) + ylab(NULL) +
      theme_pub(10) + theme(axis.text.y = element_text(size = 7.5, face = "italic"))
    save_fig3(p, paste0(which, "_KO_stratified_top10_", g), w = 6.5, h = 5.2, dir = OUT)
  }
}

cat("\n===== Step3n KO stratified 完成 =====\n")
