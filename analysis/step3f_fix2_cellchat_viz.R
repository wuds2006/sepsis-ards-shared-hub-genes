# ============================================================
# Step 3f-fix2 —— 修正 diff heatmap (ComplexHeatmap) + 通路强度对比
# ============================================================
suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(ragg)
  library(grid)
  library(ComplexHeatmap)
})
draw <- ComplexHeatmap::draw
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "08_singlecell/GSE167363/cellchat")

cc_HC  <- readRDS(file.path(OUT, "CellChat_HC.rds"))
cc_Sep <- readRDS(file.path(OUT, "CellChat_Sepsis.rds"))
cellchat_merged <- readRDS(file.path(OUT, "CellChat_merged.rds"))

# ---------- diff heatmap（ComplexHeatmap 对象）----------
save_heatmap3 <- function(hm, name, w = 8, h = 6, dir = ".", res = 500) {
  # ComplexHeatmap -> 设备输出
  pdf(file.path(dir, paste0(name, ".pdf")), width = w, height = h)
  draw(hm)
  invisible(dev.off())
  png(file.path(dir, paste0(name, ".png")), width = w, height = h, units = "in", res = 300)
  draw(hm)
  invisible(dev.off())
  agg_tiff(file.path(dir, paste0(name, ".tiff")), width = w, height = h,
           units = "in", res = res, compression = "lzw")
  draw(hm)
  invisible(dev.off())
}

tryCatch({
  hm <- netVisual_heatmap(cellchat_merged)
  save_heatmap3(hm, "diff_heatmap", dir = OUT)
  cat("diff_heatmap 重出完成\n")
}, error = function(e) cat("diff heatmap 失败:", conditionMessage(e), "\n"))

# 交互数量差 heatmap
tryCatch({
  hm2 <- netVisual_heatmap(cellchat_merged, measure = "weight")
  save_heatmap3(hm2, "diff_heatmap_weight", dir = OUT)
}, error = function(e) cat("weight heatmap 失败:", conditionMessage(e), "\n"))

# ---------- 通路强度对比（column 为 pathway_name）----------
tab <- data.frame()
for (grp in c("HC", "Sepsis")) {
  cc <- if (grp == "HC") cc_HC else cc_Sep
  df <- subsetCommunication(cc)
  agg <- aggregate(prob ~ pathway_name, data = df, FUN = sum)
  agg$group <- grp
  tab <- rbind(tab, agg)
}
pw_wide <- reshape(tab, idvar = "pathway_name", timevar = "group", direction = "wide")
colnames(pw_wide) <- gsub("prob\\.", "", colnames(pw_wide))
pw_wide$HC[is.na(pw_wide$HC)] <- 0
pw_wide$Sepsis[is.na(pw_wide$Sepsis)] <- 0
pw_wide$log2FC_Sepsis_HC <- log2((pw_wide$Sepsis + 1e-10) / (pw_wide$HC + 1e-10))
pw_wide <- pw_wide[order(pw_wide$log2FC_Sepsis_HC, decreasing = TRUE), ]
write.csv(pw_wide, file.path(OUT, "pathway_strength_HC_vs_Sepsis.csv"), row.names = FALSE)
cat("\n通路强度变化 Top10 (Sepsis 上调):\n"); print(head(pw_wide, 10))
cat("\n通路强度变化 Top10 (HC 上调):\n"); print(tail(pw_wide, 10))

tp <- rbind(head(pw_wide, 15), tail(pw_wide, 15))
tp <- tp[!duplicated(tp$pathway_name), ]
tp$dir <- ifelse(tp$log2FC_Sepsis_HC > 0, "Sepsis_enriched", "HC_enriched")
tp$pathway_name <- factor(tp$pathway_name, levels = rev(tp$pathway_name))
pbar <- ggplot(tp, aes(x = log2FC_Sepsis_HC, y = pathway_name, fill = dir)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(Sepsis_enriched = "#D99AAE", HC_enriched = "#8FA2D6")) +
  geom_vline(xintercept = 0, color = "#5B608C", linewidth = 0.5) +
  ggtitle("Signaling Pathway Strength: Sepsis vs HC (log2FC)") +
  xlab("log2 FC (Sepsis / HC)") + ylab(NULL) + theme_pub() + theme(legend.position = "none")

suppressWarnings(ggsave(file.path(OUT, "pathway_diff_bar.pdf"), pbar, width = 8, height = 8, device = "pdf"))
suppressWarnings(ggsave(file.path(OUT, "pathway_diff_bar.png"), pbar, width = 8, height = 8, dpi = 300))
agg_tiff(file.path(OUT, "pathway_diff_bar.tiff"), width = 8, height = 8,
         units = "in", res = 1000, compression = "lzw")
print(pbar); invisible(dev.off())

# ---------- 关键通路圈图（单组, base 图形设备）----------
CT_COLS <- c(
  "ClassicalMono"      = "#6A5C95",
  "IntermediateMono"   = "#D99AAE",
  "NonclassicalMono"   = "#E7A889",
  "InflammatoryMono"   = "#5B608C",
  "ISGMono"            = "#8FA2D6",
  "MonoDC"             = "#A9CBB8",
  "pDC"                = "#CFE3E6"
)
for (pw in c("MIF", "CCL", "TNF")) {
  for (grp in c("HC", "Sepsis")) {
    cc <- if (grp == "HC") cc_HC else cc_Sep
    gs <- table(cc@idents)
    cols <- CT_COLS[intersect(names(gs), names(CT_COLS))]
    tryCatch({
      png(file.path(OUT, paste0("pathway_", pw, "_", grp, ".png")), width = 1800, height = 1600, res = 300)
      netVisual_aggregate(cc, signaling = pw, layout = "circle", color.use = cols)
      dev.off()
    }, error = function(e) { dev.off(); cat(pw, grp, "跳过:", conditionMessage(e), "\n") })
  }
}

cat("\n===== Step3f-fix2 完成 =====\n")
