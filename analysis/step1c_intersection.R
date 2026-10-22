# ============================================================
# Step 1c —— DEG ∩ WGCNA 交集基因 CGs + 双病种共享基因 + GO/KEGG
# 流程: ARDS_CGs = DEG_ARDS ∩ WGCNA_ARDS
#       Sepsis_CGs = DEG_Sepsis ∩ WGCNA_Sepsis
#       共享 CGs = ARDS_CGs ∩ Sepsis_CGs  (共同疾病基因)
# ============================================================
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(VennDiagram)
  library(ggplot2)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "03_intersection")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# 1. 读 DEG
deg_ards    <- read.csv(file.path(BASE, "01_DEG", "GSE76293_DEG.csv"), stringsAsFactors = FALSE)
deg_sepsis  <- read.csv(file.path(BASE, "01_DEG", "GSE69528_DEG.csv"), stringsAsFactors = FALSE)

# 2. 读 WGCNA 疾病模块基因
wgcna_ards    <- readLines(file.path(BASE, "02_WGCNA", "GSE76293_disease_module_genes.txt"))
wgcna_sepsis  <- readLines(file.path(BASE, "02_WGCNA", "GSE69528_disease_module_genes.txt"))

deg_ards_g   <- deg_ards$gene
deg_sepsis_g <- deg_sepsis$gene

# 3. 交集
ards_cgs    <- intersect(deg_ards_g, wgcna_ards)
sepsis_cgs  <- intersect(deg_sepsis_g, wgcna_sepsis)
shared_cgs  <- intersect(ards_cgs, sepsis_cgs)

cat("ARDS DEG:", length(deg_ards_g), "| ARDS WGCNA:", length(wgcna_ards),
    "| ARDS_CGs:", length(ards_cgs), "\n")
cat("Sepsis DEG:", length(deg_sepsis_g), "| Sepsis WGCNA:", length(wgcna_sepsis),
    "| Sepsis_CGs:", length(sepsis_cgs), "\n")
cat("共享 CGs (共同疾病基因):", length(shared_cgs), "\n")

# 保存
writeLines(ards_cgs,   file.path(OUT, "ARDS_CGs.txt"))
writeLines(sepsis_cgs, file.path(OUT, "Sepsis_CGs.txt"))
writeLines(shared_cgs, file.path(OUT, "Shared_CGs.txt"))
write.csv(data.frame(gene = shared_cgs), file.path(OUT, "Shared_CGs.csv"), row.names = FALSE)

# 4. 韦恩图 (VennDiagram)
venn_list <- list(ARDS_DEG = deg_ards_g, ARDS_WGCNA = wgcna_ards,
                  Sepsis_DEG = deg_sepsis_g, Sepsis_WGCNA = wgcna_sepsis)
pdf(file.path(OUT, "venn_4sets.pdf"), width = 8, height = 8)
grid.newpage()
draw.quad.venn(area1 = length(deg_ards_g), area2 = length(wgcna_ards),
               area3 = length(deg_sepsis_g), area4 = length(wgcna_sepsis),
               n12 = length(intersect(deg_ards_g, wgcna_ards)),
               n13 = length(intersect(deg_ards_g, deg_sepsis_g)),
               n14 = length(intersect(deg_ards_g, wgcna_sepsis)),
               n23 = length(intersect(wgcna_ards, deg_sepsis_g)),
               n24 = length(intersect(wgcna_ards, wgcna_sepsis)),
               n34 = length(intersect(deg_sepsis_g, wgcna_sepsis)),
               n123 = length(Reduce(intersect, list(deg_ards_g, wgcna_ards, deg_sepsis_g))),
               n124 = length(Reduce(intersect, list(deg_ards_g, wgcna_ards, wgcna_sepsis))),
               n134 = length(Reduce(intersect, list(deg_ards_g, deg_sepsis_g, wgcna_sepsis))),
               n234 = length(Reduce(intersect, list(wgcna_ards, deg_sepsis_g, wgcna_sepsis))),
               n1234 = length(shared_cgs),
               category = c("ARDS DEG", "ARDS WGCNA", "Sepsis DEG", "Sepsis WGCNA"),
               fill = c(PAL["peach"], PAL["green"], PAL["periwinkle"], PAL["lightblue"]),
               alpha = 0.5, cat.cex = 1.2, cex = 1.4)
dev.off()

# 5. GO/KEGG 富集 (对共享 CGs)
enrich_and_save <- function(genes, prefix, title) {
  if (length(genes) == 0) { cat(prefix, "无基因，跳过\n"); return(NULL) }
  eg <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (nrow(eg) == 0) { cat(prefix, "无映射，跳过\n"); return(NULL) }
  ego <- enrichGO(gene = eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  ekegg <- enrichKEGG(gene = eg$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    write.csv(as.data.frame(ego), file.path(OUT, paste0(prefix, "_GO_BP.csv")), row.names = FALSE)
    p <- dotplot(ego, showCategory = 15) + ggtitle(title) + theme_pub()
    save_fig(p, paste0(prefix, "_GO_dot"), w = 8, h = 7, dir = OUT)
  }
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    write.csv(as.data.frame(ekegg), file.path(OUT, paste0(prefix, "_KEGG.csv")), row.names = FALSE)
  }
  cat(prefix, "富集完成: GO", if(!is.null(ego)) nrow(as.data.frame(ego)) else 0,
      "条, KEGG", if(!is.null(ekegg)) nrow(as.data.frame(ekegg)) else 0, "条\n")
  list(ego = ego, ekegg = ekegg)
}

enrich_and_save(shared_cgs, "Shared_CGs", "Shared CGs: GO-BP")

cat("\n===== Step1c 完成 =====\n")
cat("共享 CGs 数:", length(shared_cgs), "\n")
cat("共享 CGs 列表:", paste(shared_cgs, collapse = ", "), "\n")
