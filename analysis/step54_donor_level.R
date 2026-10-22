# ============================================================
# Step 54 —— 供体级单细胞统计（审稿 W4/问题4，R 部分）
# 1) GSE167363 每供体 UPP1-pseudotime 轨迹相关（Spearman）
# 2) 每供体 pseudobulk CEBPB/UPP1 表达（log2CPM）
# 3) pseudobulk DESeq2：GSE167363 (Sepsis vs HC) 与 GSE151263 (ARDS vs Sepsis)
# ============================================================
suppressPackageStartupMessages({ library(DESeq2); library(Matrix) })
BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "08_singlecell", "donor_level")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
HUB6 <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")

# ---------- GSE167363 装载 ----------
cd  <- file.path(BASE, "08_singlecell", "GSE167363", "input")
cnt <- readMM(file.path(cd, "counts.mtx"))            # genes x cells
genes <- sub("\t.*", "", readLines(file.path(cd, "genes.tsv")))
meta <- read.csv(file.path(cd, "meta.csv"), row.names = 1)
stopifnot(ncol(cnt) == nrow(meta))
rownames(cnt) <- genes; colnames(cnt) <- rownames(meta)
pt <- read.csv(file.path(BASE, "08_singlecell", "monocle3", "monocle3_pseudotime_GSE167363.csv"),
               row.names = 1)
lib <- colSums(cnt)
lcpm <- log2(t(t(cnt) / pmax(lib, 1) * 1e6) + 1)

# ---------- 1) 每供体轨迹相关 + 2) 每供体 pseudobulk ----------
samples <- c("GSM5102900_HC1", "GSM5102901_HC2",
             "GSM5511353_S2_T0", "GSM5511354_S2_T6", "GSM5511355_S3_T0")
donor_of <- meta$sample
donor_of[grepl("S2_T", meta$sample)] <- "S2"          # 同一患者两时点合并（供体级口径）
rows <- list()
for (smp in samples) {
  cells <- colnames(cnt)[meta$sample == smp]
  cells_pt <- intersect(cells, rownames(pt))
  upp1 <- lcpm["UPP1", cells_pt]
  pv <- pt[cells_pt, "pseudotime"]
  rho <- suppressWarnings(cor(upp1, pv, method = "spearman"))
  ct_p <- tryCatch(suppressWarnings(cor.test(upp1, pv, method = "spearman")$p.value),
                   error = function(e) NA)
  pb <- rowSums(cnt[, cells])
  rows[[length(rows) + 1]] <- data.frame(
    donor = sub("^GSM[0-9]+_", "", smp), sample = smp,
    group = unique(meta$group[meta$sample == smp]), n_cells = length(cells),
    UPP1_pseudotime_rho = round(rho, 3), rho_P = signif(ct_p, 3),
    CEBPB_pb_log2CPM = round(log2(sum(pb["CEBPB"]) / sum(pb) * 1e6 + 1), 3),
    UPP1_pb_log2CPM  = round(log2(sum(pb["UPP1"]) / sum(pb) * 1e6 + 1), 3))
}
don_tab <- do.call(rbind, rows)
write.csv(don_tab, file.path(OUT, "GSE167363_donor_level_UPP1_CEBPB.csv"), row.names = FALSE)
print(don_tab)

# ---------- 3) pseudobulk DESeq2 ----------
make_pb <- function(cnt, meta, donor_id) {
  cts <- sort(unique(meta$celltype))
  dmap <- tapply(meta$group, donor_id, function(x) names(which.max(table(x))))
  out <- list()
  for (ctu in cts) for (d in unique(donor_id)) {
    cells <- colnames(cnt)[meta$celltype == ctu & donor_id == d]
    if (length(cells) >= 20) out[[paste(ctu, d, sep = "__")]] <- rowSums(cnt[, cells, drop = FALSE])
  }
  pb <- do.call(cbind, out)
  ctv <- sub("__.*", "", colnames(pb)); don <- sub(".*__", "", colnames(pb))
  list(pb = pb, ct = ctv, donor = don, group = dmap[don])
}
pb_de <- function(pbo, g1, g0, tag) {
  res_all <- list()
  for (ctu in unique(pbo$ct)) {
    idx <- which(pbo$ct == ctu)
    m <- pbo$pb[, idx, drop = FALSE]; g <- pbo$group[idx]
    if (length(unique(g)) < 2 || min(table(g)) < 2) next
    dds <- tryCatch({
      coldata <- data.frame(group = factor(g))
      rownames(coldata) <- colnames(m)
      dd <- DESeqDataSetFromMatrix(round(m), coldata, ~ group)
      DESeq(dd, quiet = TRUE)
    }, error = function(e) { cat(sprintf("  [DESeq2 %s/%s ERROR] %s\n", tag, ctu, conditionMessage(e))); NULL })
    if (is.null(dds)) next
    r <- as.data.frame(results(dds, contrast = c("group", g1, g0)))
    hit <- HUB6[HUB6 %in% rownames(r)]
    if (length(hit)) res_all[[length(res_all) + 1]] <- cbind(celltype = ctu, gene = hit, r[hit, ])
  }
  if (length(res_all)) {
    x <- do.call(rbind, res_all)
    write.csv(x[, c("celltype", "gene", "baseMean", "log2FoldChange", "pvalue", "padj")],
              file.path(OUT, sprintf("pseudobulk_DESeq2_%s.csv", tag)), row.names = FALSE)
    cat(sprintf("\n[%s] pseudobulk DESeq2: %d celltype-gene 组合\n", tag, nrow(x)))
    print(x[, c("celltype", "gene", "log2FoldChange", "padj")])
  } else cat(sprintf("\n[%s] 无可用 pseudobulk 组合\n", tag))
}

donor_id_167 <- donor_of
pb1 <- make_pb(cnt, meta, donor_id_167)
pb_de(pb1, "Sepsis", "HC", "GSE167363_Sepsis_vs_HC")

# ---------- GSE151263 (ARDS vs Sepsis) ----------
cd2 <- file.path(BASE, "08_singlecell", "cellchat_input")
cnt2 <- readMM(file.path(cd2, "counts.mtx"))
genes2 <- sub("\t.*", "", readLines(file.path(cd2, "genes.tsv")))
meta2 <- read.csv(file.path(cd2, "meta.csv"), row.names = 1)
stopifnot(ncol(cnt2) == nrow(meta2))
rownames(cnt2) <- genes2; colnames(cnt2) <- rownames(meta2)
pb2 <- make_pb(cnt2, meta2, meta2$sample)
pb_de(pb2, "ARDS", "Sepsis", "GSE151263_ARDS_vs_Sepsis")

cat("\n=== Step 54 完成 ===\n")
