# ============================================================
# Step 55 —— WGCNA 跨发现队列模块保留性（审稿详细评语：模块稳定性）
# 参考 = GSE76293 (ARDS)，测试 = GSE69528 (Sepsis)，及反向
# 关注疾病相关模块（ARDS turquoise / Sepsis blue）的 Zsummary
# ============================================================
suppressPackageStartupMessages({ library(WGCNA) })
BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "02_WGCNA", "preservation")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
options(stringsAsFactors = FALSE)

load_set <- function(ds) {
  expr <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_gene_expr.rds")))
  mods <- read.csv(file.path(BASE, "02_WGCNA", paste0(ds, "_gene_module.csv")))
  mads <- apply(expr, 1, mad, na.rm = TRUE)
  expr <- expr[order(mads, decreasing = TRUE), , drop = FALSE]
  expr <- expr[1:min(3000, nrow(expr)), , drop = FALSE]   # 与 step1b 相同口径
  list(expr = expr, mods = mods)
}
a <- load_set("GSE76293")   # ARDS
s <- load_set("GSE69528")   # Sepsis

run_pres <- function(ref, test, ref_name, test_name) {
  # 仅用参考网络中有模块标签且两队列均有表达的基因
  common <- intersect(intersect(rownames(ref$expr), rownames(test$expr)), ref$mods$gene)
  mref <- ref$mods$module[match(common, ref$mods$gene)]
  names(mref) <- common
  stopifnot(!any(is.na(mref)))
  ref_expr <- ref$expr[common, , drop = FALSE]
  test_expr <- test$expr[common, , drop = FALSE]
  # WGCNA 内部 .checkExpr 访问 multiExpr[[set]]$data —— 组件名必须是 data
  multi <- list(ref = list(data = t(ref_expr)),
                test = list(data = t(test_expr)))
  # 基因顺序一致（列为基因）
  multi$test$data <- multi$test$data[, colnames(multi$ref$data)]
  # 先用 perm=0 快速调用拿 quality 行名 = 模块真实顺序（WGCNA 会追加 gold 占位模块）
  mp0 <- modulePreservation(multiData = multi, multiColor = list(ref = mref, test = mref),
                            dataIsExpr = TRUE, networkType = "unsigned", corFnc = "cor",
                            referenceNetworks = 1, nPermutations = 0,
                            maxModuleSize = 5000, savePermutedStatistics = FALSE,
                            plotInterpolation = FALSE, verbose = 0)
  mods <- rownames(mp0$observed[[1]]$quality[[2]])
  stopifnot(!is.null(mods), !any(is.na(mods)))
  mp <- modulePreservation(multiData = multi, multiColor = list(ref = mref, test = mref),
                           dataIsExpr = TRUE, networkType = "unsigned",
                           corFnc = "cor", referenceNetworks = 1,
                           nPermutations = 200, randomSeed = 2026,
                           maxModuleSize = 5000, savePermutedStatistics = FALSE,
                           plotInterpolation = FALSE, verbose = 0, indent = 0)
  # 输出结构（WGCNA 1.9+）：mp$preservation$Z[[1]][[2]]$moduleSize / $Zsummary.pres
  zstat <- mp$preservation$Z[[1]][[2]]
  ms <- as.numeric(zstat$moduleSize)
  stopifnot(length(mods) == length(ms), length(ms) == length(zstat$Zsummary.pres))
  res <- data.frame(module = mods,
                    moduleSize = ms,
                    Zsummary = round(as.numeric(zstat$Zsummary.pres), 2),
                    row.names = NULL)
  res <- res[order(-res$Zsummary), ]
  write.csv(res, file.path(OUT, sprintf("modulePreservation_%s_in_%s.csv",
                                        ref_name, test_name)), row.names = FALSE)
  cat(sprintf("\n===== %s 模块在 %s 中的保留性 =====\n", ref_name, test_name))
  print(res)
  res
}

r1 <- run_pres(a, s, "GSE76293", "GSE69528")
r2 <- run_pres(s, a, "GSE69528", "GSE76293")

cat("\n疾病相关模块要点：\n")
cat("ARDS turquoise 在 Sepsis 队列 Zsummary =",
    r1$Zsummary[r1$module == "turquoise"], "\n")
cat("Sepsis blue 在 ARDS 队列 Zsummary =",
    r2$Zsummary[r2$module == "blue"], "\n")
cat("\n=== Step 55 完成 ===\n")
