# ============================================================
# Step 3n-3c —— 分组 scTenifoldKnk 结果的组间比较 + GO 富集 + 图件
# 读取 knockout_stratified/KO_{ds}_{grp}_{gene}_diffRegulation.csv 全量结果
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2); library(ragg); library(clusterProfiler); library(org.Hs.eg.db)
})
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/common.R")
set.seed(42)
BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT <- file.path(BASE, "08_singlecell/knockout_stratified")
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

DS <- list(
  list(which = "GSE151263", groups = c("ARDS", "Sepsis"), genes = c("UPP1", "EXOSC4")),
  list(which = "GSE167363", groups = c("Sepsis", "HC"), genes = c("UPP1", "EXOSC4", "VNN1", "CDK5RAP2"))
)

results <- list()
for (ds in DS) {
  for (grp in ds$groups) {
    results[[ds$which]][[grp]] <- list()
    for (g in ds$genes) {
      f <- file.path(OUT, sprintf("KO_%s_%s_%s_diffRegulation.csv", ds$which, grp, g))
      if (file.exists(f)) {
        dr <- read.csv(f, stringsAsFactors = FALSE)
        results[[ds$which]][[grp]][[g]] <- dr
      }
    }
  }
}

# ---- 组间比较表 ----
comp_rows <- list()
enrich_list <- list()
for (ds in DS) {
  which <- ds$which; groups <- ds$groups
  for (g in ds$genes) {
    a <- results[[which]][[groups[1]]][[g]]
    b <- results[[which]][[groups[2]]][[g]]
    if (is.null(a) || is.null(b)) {
      cat(sprintf("%s %s: 缺一组结果, 跳过比较\n", which, g)); next
    }
    common <- intersect(a$gene, b$gene)
    ct <- suppressWarnings(cor.test(a$FC[match(common, a$gene)],
                                    b$FC[match(common, b$gene)], method = "spearman"))
    sigA <- a$gene[a$p.adj < 0.05]; sigB <- b$gene[b$p.adj < 0.05]
    jac <- length(intersect(sigA, sigB)) / max(length(union(sigA, sigB)), 1)
    topA <- head(a$gene, 150); topB <- head(b$gene, 150)
    jac_top <- length(intersect(topA, topB)) / length(union(topA, topB))
    comp_rows[[length(comp_rows) + 1]] <- data.frame(
      dataset = which, gene = g,
      n_affected_case = nrow(a), n_affected_ctrl = nrow(b),
      n_common_genes = length(common),
      spearman_rho = round(unname(ct$estimate), 3),
      spearman_p = signif(ct$p.value, 3),
      jaccard_top150 = round(jac_top, 3),
      jaccard_sig = round(jac, 3))
    # GO 富集（每组 top150 受影响基因）
    for (grp in groups) {
      dr <- results[[which]][[grp]][[g]]
      genes <- head(dr$gene, 150)
      eg <- mapIds(org.Hs.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL",
                   multiVals = "first")
      eg <- eg[!is.na(eg)]
      if (length(eg) < 8) next
      ego <- tryCatch(
        enrichGO(gene = eg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                 ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.1, readable = TRUE),
        error = function(e) NULL)
      if (is.null(ego)) next
      edf <- as.data.frame(ego)
      if (nrow(edf) == 0) next
      edf$group <- grp; edf$dataset <- which; edf$geneKO <- g
      edf$GeneRatioS <- vapply(strsplit(as.character(edf$GeneRatio), "/"),
                               function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
      enrich_list[[length(enrich_list) + 1]] <- edf
    }
  }
}
comp <- do.call(rbind, comp_rows)
write.csv(comp, file.path(OUT, "KO_stratified_group_comparison.csv"), row.names = FALSE)
print(comp)

if (length(enrich_list)) {
  enrich_all <- do.call(rbind, enrich_list)
  write.csv(enrich_all, file.path(OUT, "KO_stratified_GO_BP.csv"), row.names = FALSE)
  for (ds in DS) {
    which <- ds$which
    sub <- enrich_all[enrich_all$dataset == which, ]
    for (g in unique(sub$geneKO)) {
      sg <- sub[sub$geneKO == g, ]
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

# ---- top10 受影响基因分组对比图 ----
for (ds in DS) {
  which <- ds$which; groups <- ds$groups
  for (g in names(results[[which]][[groups[1]]])) {
    a <- results[[which]][[groups[1]]][[g]]
    b <- results[[which]][[groups[2]]][[g]]
    if (is.null(a) || is.null(b)) next
    ta <- head(a, 10); ta$grp_f <- factor(groups[1], levels = groups)
    tb <- head(b, 10); tb$grp_f <- factor(groups[2], levels = groups)
    tt <- rbind(ta, tb)
    tt$gene <- factor(tt$gene, levels = rev(unique(tt$gene)))
    tt$dir <- ifelse(tt$distance > 0, "up", "down")
    p <- ggplot(tt, aes(x = grp_f, y = gene)) +
      geom_point(aes(size = abs(distance), fill = dir), shape = 21, color = "grey20", stroke = 0.4) +
      scale_fill_manual(values = c(up = "#D99AAE", down = "#8FA2D6"), name = "Direction") +
      scale_size_continuous(range = c(1.5, 6), name = "|dRegulation|") +
      ggtitle(paste0(which, ": scTenifoldKnk KO ", g, " - top 10 perturbed genes per group")) +
      xlab(NULL) + ylab(NULL) +
      theme_pub(10) + theme(axis.text.y = element_text(size = 7.5, face = "italic"))
    save_fig3(p, paste0(which, "_KO_stratified_top10_", g), w = 6.5, h = 5.2, dir = OUT)
  }
}

# ---- 组内检出率与可敲性（含不可敲记录）----
feas <- list()
for (ds in DS) {
  which <- ds$which
  for (grp in ds$groups) {
    for (g in HUB) {
      f_ok <- file.path(OUT, sprintf("KO_%s_%s_%s_diffRegulation.csv", which, grp, g))
      f_no <- file.path(OUT, sprintf("KO_%s_%s_%s_NOTKNOCKABLE.csv", which, grp, g))
      feas[[length(feas) + 1]] <- data.frame(
        dataset = which, group = grp, gene = g,
        status = if (file.exists(f_ok)) "knocked_out" else
          (if (file.exists(f_no)) "not_knockable_in_group" else "low_expression_skipped"))
    }
  }
}
write.csv(do.call(rbind, feas), file.path(OUT, "KO_stratified_feasibility.csv"), row.names = FALSE)
cat("\n===== 比较与富集完成 =====\n")
