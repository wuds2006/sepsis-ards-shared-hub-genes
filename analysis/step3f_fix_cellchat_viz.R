# ============================================================
# Step 3f-fix —— CellChat 可视化（从已保存 rds 重载）
# ============================================================
suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(patchwork)
  library(ragg)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "08_singlecell/GSE167363/cellchat")

save_fig3 <- function(p, name, w = 7, h = 6, dir = ".", line_art = FALSE) {
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf"))
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300))
  tdpi <- if (line_art) 1000 else 500
  agg_tiff(file.path(dir, paste0(name, ".tiff")), width = w, height = h,
           units = "in", res = tdpi, compression = "lzw")
  print(p); invisible(dev.off())
  invisible(p)
}

cc_HC  <- readRDS(file.path(OUT, "CellChat_HC.rds"))
cc_Sep <- readRDS(file.path(OUT, "CellChat_Sepsis.rds"))

CT_COLS <- c(
  "ClassicalMono"      = "#6A5C95",
  "IntermediateMono"   = "#D99AAE",
  "NonclassicalMono"   = "#E7A889",
  "InflammatoryMono"   = "#5B608C",
  "ISGMono"            = "#8FA2D6",
  "MonoDC"             = "#A9CBB8",
  "pDC"                = "#CFE3E6",
  "Cycling"            = "#A9CBB8"
)

# ---------- 圈图（count/weight, HC/Sepsis）----------
for (grp in c("HC", "Sepsis")) {
  cc <- if (grp == "HC") cc_HC else cc_Sep
  gs <- table(cc@idents)
  cols <- CT_COLS[intersect(names(gs), names(CT_COLS))]
  for (s in c("count", "weight")) {
    png(file.path(OUT, paste0("net_", s, "_", grp, ".png")), width = 2000, height = 1800, res = 300)
    par(mar = c(2, 2, 2, 1))
    netVisual_circle(cc@net[[s]], vertex.weight = as.numeric(gs), weight.scale = TRUE,
                     label.edge = FALSE, title.name = paste(toupper(s), "-", grp),
                     vertex.label.cex = 0.8, color.use = cols)
    dev.off()
  }
}

# ---------- 合并比较 ----------
cat("合并比较 HC vs Sepsis...\n")
object.list <- list(HC = cc_HC, Sepsis = cc_Sep)
cellchat_merged <- mergeCellChat(object.list, add.names = c("HC", "Sepsis"), cell.prefix = TRUE)
saveRDS(cellchat_merged, file.path(OUT, "CellChat_merged.rds"))

tryCatch({
  gg1 <- netVisual_heatmap(cellchat_merged)
  save_fig3(gg1, "diff_heatmap", w = 8, h = 6, dir = OUT)
}, error = function(e) cat("diff heatmap 失败:", conditionMessage(e), "\n"))

tryCatch({
  gg2 <- rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
  save_fig3(gg2, "rankNet_pathways", w = 7, h = 8, dir = OUT)
}, error = function(e) cat("rankNet 失败:", conditionMessage(e), "\n"))

tryCatch({
  pb <- netVisual_bubble(cellchat_merged, comparison = c(1, 2), angle.x = 45)
  save_fig3(pb, "diff_LR_bubble", w = 13, h = 10, dir = OUT)
}, error = function(e) cat("bubble 图失败:", conditionMessage(e), "\n"))

# ---------- 关键炎症通路两组对比圈图 ----------
for (pw in c("MIF", "CCL", "CXCL", "IL1", "IL6", "TNF")) {
  ok <- tryCatch({
    p1 <- netVisual_aggregate(cc_HC, signaling = pw, layout = "circle",
                              color.use = CT_COLS[intersect(names(table(cc_HC@idents)), names(CT_COLS))])
    p2 <- netVisual_aggregate(cc_Sep, signaling = pw, layout = "circle",
                              color.use = CT_COLS[intersect(names(table(cc_Sep@idents)), names(CT_COLS))])
    save_fig3(p1 + p2, paste0("pathway_", pw, "_HCvsSepsis"), w = 12, h = 6, dir = OUT)
    TRUE
  }, error = function(e) { cat(pw, "通路图跳过:", conditionMessage(e), "\n"); FALSE })
}

# ---------- 通路活性汇总表 ----------
tab <- data.frame()
for (grp in c("HC", "Sepsis")) {
  cc <- if (grp == "HC") cc_HC else cc_Sep
  pn <- tryCatch(cc@netPVE$prob, error = function(e) NULL)
  # 通路总强度（per pathway, sum over LR pairs）
  df <- tryCatch(subsetCommunication(cc, slot.name = "netP"), error = function(e) NULL)
  if (!is.null(df)) {
    agg <- aggregate(prob ~ pathway, data = df, FUN = sum)
    agg$group <- grp
    tab <- rbind(tab, agg)
  }
}
if (nrow(tab)) {
  pw_wide <- reshape(tab, idvar = "pathway", timevar = "group", direction = "wide")
  colnames(pw_wide) <- gsub("prob\\.", "", colnames(pw_wide))
  pw_wide$log2FC_Sepsis_HC <- log2((pw_wide$Sepsis + 1e-10) / (pw_wide$HC + 1e-10))
  pw_wide <- pw_wide[order(pw_wide$log2FC_Sepsis_HC, decreasing = TRUE), ]
  write.csv(pw_wide, file.path(OUT, "pathway_strength_HC_vs_Sepsis.csv"), row.names = FALSE)
  cat("\n通路强度变化 Top10 (Sepsis vs HC):\n")
  print(head(pw_wide, 10))

  # 通路差异条形图
  tp <- head(rbind(pw_wide[order(pw_wide$log2FC_Sepsis_HC, decreasing = TRUE), ][1:15, ],
                   pw_wide[order(pw_wide$log2FC_Sepsis_HC), ][1:15, ]), 25)
  tp$dir <- ifelse(tp$log2FC_Sepsis_HC > 0, "Sepsis_enriched", "HC_enriched")
  tp$pathway <- factor(tp$pathway, levels = rev(tp$pathway))
  pbar <- ggplot(tp, aes(x = log2FC_Sepsis_HC, y = pathway, fill = dir)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(Sepsis_enriched = "#D99AAE", HC_enriched = "#8FA2D6")) +
    geom_vline(xintercept = 0, color = "#5B608C", linewidth = 0.5) +
    ggtitle("Signaling Pathway Strength: Sepsis vs HC (log2FC)") +
    xlab("log2 FC (Sepsis / HC)") + ylab(NULL) + theme_pub() + NoLegend()
  save_fig3(pbar, "pathway_diff_bar", w = 8, h = 7, dir = OUT, line_art = TRUE)
}

cat("\n===== Step3f-fix CellChat 可视化完成 =====\n")
