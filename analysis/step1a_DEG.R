# ============================================================
# Step 1a —— 差异表达分析（limma）
# GSE76293(ARDS) 与 GSE69528(Sepsis) 从 GEO 原始探针级矩阵重跑
# 阈值: |log2FC| > 1 且 adj.P.Val < 0.05
# ============================================================
suppressPackageStartupMessages({
  library(limma)
})
source("scripts/common.R")

BASE    <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
PAPER   <- "/Users/mac/Desktop/2026年撰写的文章/2026论文撰写一Acute Respiratory Distress Syndrome and Sepsis"
OUT     <- file.path(BASE, "01_DEG")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# 注释映射
gpl570_col <- c(probe = "ID", gene = "Gene Symbol")
gpl10558_col <- c(probe = "ID", gene = "Symbol")

load_annot <- function(path, probe_col, gene_col) {
  a <- read.delim(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE,
                  comment.char = "#", quote = "", sep = "\t")
  a <- a[, c(probe_col, gene_col)]
  colnames(a) <- c("probe", "gene")
  a$gene <- sub(" ///.*", "", a$gene)     # 多基因取首个
  a <- a[a$gene != "" & !is.na(a$gene), ]
  a
}

collapse_probes <- function(expr, annot) {
  # 只保留有注释的探针
  common <- intersect(rownames(expr), annot$probe)
  expr <- expr[common, , drop = FALSE]
  gmap <- setNames(annot$gene[match(common, annot$probe)], common)
  # 每个基因取平均表达最高的探针
  gmean <- rowMeans(expr, na.rm = TRUE)
  df <- data.frame(probe = common, gene = gmap, mean = gmean, stringsAsFactors = FALSE)
  df <- df[order(df$mean, decreasing = TRUE), ]
  df <- df[!duplicated(df$gene), ]
  out <- expr[df$probe, , drop = FALSE]
  rownames(out) <- df$gene[match(df$probe, df$probe)]
  out
}

run_deg <- function(ds, disease, control, annot_file, annot_cols) {
  cat("\n##########", ds, "##########\n")
  expr  <- read_expr(file.path(OUT, paste0(ds, "_expr.txt")))
  grp   <- read.csv(file.path(OUT, paste0(ds, "_groups.csv")), stringsAsFactors = FALSE)

  # 过滤相关样本
  keep <- grp$group %in% c(disease, control)
  expr <- expr[, grp$sample[keep], drop = FALSE]
  grp  <- grp[keep, ]
  grp$group <- factor(grp$group, levels = c(control, disease))
  cat(sprintf("样本: %s=%d, %s=%d\n", disease, sum(grp$group==disease),
              control, sum(grp$group==control)))

  # 探针 -> 基因
  annot <- load_annot(annot_file, annot_cols["probe"], annot_cols["gene"])
  gene_expr <- collapse_probes(expr, annot)
  cat(sprintf("基因数(去重后): %d\n", nrow(gene_expr)))

  # limma
  design <- model.matrix(~0 + grp$group)
  colnames(design) <- levels(grp$group)
  cm <- makeContrasts(contrasts = paste(disease, control, sep = "-"), levels = design)
  fit <- lmFit(gene_expr, design)
  fit <- contrasts.fit(fit, cm)
  fit <- eBayes(fit, trend = FALSE)
  res <- topTable(fit, coef = 1, number = Inf, sort.by = "P")
  res$gene <- rownames(res)
  res$regulation <- ifelse(res$logFC > 1 & res$adj.P.Val < 0.05, "Up",
                    ifelse(res$logFC < -1 & res$adj.P.Val < 0.05, "Down", "NS"))

  deg <- res[res$regulation != "NS", ]
  n_up <- sum(res$regulation == "Up")
  n_dn <- sum(res$regulation == "Down")
  cat(sprintf("DEG: 上调 %d, 下调 %d, 合计 %d\n", n_up, n_dn, nrow(deg)))

  write.csv(res, file.path(OUT, paste0(ds, "_DEG_all.csv")), row.names = FALSE)
  write.csv(deg, file.path(OUT, paste0(ds, "_DEG.csv")), row.names = FALSE)

  # 火山图
  vp <- ggplot(res, aes(x = logFC, y = -log10(P.Value), color = regulation)) +
    geom_point(size = 0.7, alpha = 0.7) +
    scale_color_manual(values = c(Up = PAL_UP, Down = PAL_DOWN, NS = PAL_NS),
                       breaks = c("Up", "Down", "NS"),
                       labels = c(paste0("Up (", n_up, ")"),
                                  paste0("Down (", n_dn, ")"), "NS")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    labs(title = paste(ds, ":", disease, "vs", control),
         x = "log2 Fold Change", y = "-log10(P value)") +
    theme_pub()
  save_fig(vp, paste0(ds, "_volcano"), w = 7, h = 6, dir = OUT)

  # 保存基因表达矩阵（供 WGCNA 与后续使用）
  saveRDS(gene_expr, file.path(OUT, paste0(ds, "_gene_expr.rds")))
  saveRDS(grp, file.path(OUT, paste0(ds, "_grp.rds")))

  list(res = res, deg = deg, gene_expr = gene_expr, grp = grp)
}

# 运行
r1 <- run_deg("GSE76293", "ARDS", "Control",
              file.path(PAPER, "ARDS/GSE76293/GPL570-55999.txt"), gpl570_col)
r2 <- run_deg("GSE69528", "Sepsis", "Control",
              file.path(PAPER, "脓毒症/GSE69528/GPL10558-50081.txt"), gpl10558_col)

cat("\n===== Step1a DEG 完成 =====\n")
cat("GSE76293 DEG 数:", nrow(r1$deg), "| GSE69528 DEG 数:", nrow(r2$deg), "\n")
