# ============================================================
# Step 3m —— 两数据集 CellChat 补件与格式规范化
# PBMC (GSE151263): 补 net_count/net_weight、pathway 网络图(MIF/MHC-I/GALECTIN)、
#                   pathway_diff_bar、pathway_strength.csv、diff_heatmap_weight
# CD14+ (GSE167363): 重出 net_count/net_weight、pathway_CCL/MIF/TNF 为三格式
# 三格式: pdf(矢量) + png(300dpi RGB) + tiff(500dpi LZW RGB)
# ============================================================
suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(ragg)
  library(grid)
})
setwd("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis")
BASE <- getwd()
NET_COLS <- c("#6A5C95", "#D99AAE", "#E7A889", "#A9CBB8", "#CFE3E6", "#8FA2D6", "#5B608C")

# ---- 三格式设备封装（CellChat base 图形用，传匿名函数闭包）----
save3_base <- function(plot_fn, name, dir, w = 7, h = 6.5, tiff_dpi = 500) {
  pdf(file.path(dir, paste0(name, ".pdf")), width = w, height = h, family = "Helvetica",
      useDingbats = FALSE)
  plot_fn(); dev.off()
  agg_png(file.path(dir, paste0(name, ".png")), width = w, height = h, units = "in",
          res = 300, background = "white")
  plot_fn(); dev.off()
  agg_tiff(file.path(dir, paste0(name, ".tiff")), width = w, height = h, units = "in",
           res = tiff_dpi, background = "white", compression = "lzw")
  plot_fn(); dev.off()
  cat("  OK:", name, "\n")
}
# ---- ggplot 对象三格式 ----
save3_gg <- function(p, name, dir, w = 8, h = 6, tiff_dpi = 500) {
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h,
                          device = "pdf"))
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h,
                          dpi = 300, bg = "white"))
  agg_tiff(file.path(dir, paste0(name, ".tiff")), width = w, height = h, units = "in",
           res = tiff_dpi, background = "white", compression = "lzw")
  print(p); dev.off()
  cat("  OK:", name, "\n")
}

process_dataset <- function(out_dir, g1, g2, pathways, ds_label) {
  cat(sprintf("\n========== %s (%s vs %s) ==========\n", ds_label, g1, g2))
  cc1 <- readRDS(file.path(out_dir, paste0("CellChat_", g1, ".rds")))
  cc2 <- readRDS(file.path(out_dir, paste0("CellChat_", g2, ".rds")))
  ccm <- readRDS(file.path(out_dir, "CellChat_merged.rds"))
  gs1 <- table(cc1@idents); gs2 <- table(cc2@idents)

  # ---- 1) net count / weight 圈图（每组）----
  for (s in c("count", "weight")) {
    for (grp in c(g1, g2)) {
      cc <- if (grp == g1) cc1 else cc2
      gs <- if (grp == g1) gs1 else gs2
      ttl <- paste(if (s == "count") "Number of interactions" else "Interaction strength",
                   "-", grp, paste0("(", ds_label, ")"))
      save3_base(function() {
        par(mar = c(1, 1, 2.5, 1))
        netVisual_circle(cc@net[[s]], vertex.weight = gs, weight.scale = TRUE,
                         label.edge = FALSE, title.name = ttl,
                         vertex.label.cex = 0.8,
                         color.use = setNames(NET_COLS[seq_along(gs)], names(gs)))
      }, paste0("net_", s, "_", grp), out_dir)
    }
  }

  # ---- 2) 关键通路网络图（每组）----
  for (pw in pathways) {
    for (grp in c(g1, g2)) {
      cc <- if (grp == g1) cc1 else cc2
      gs <- if (grp == g1) gs1 else gs2
      if (!(pw %in% dimnames(cc@netP$prob)[[3]])) {
        cat("  跳过(通路不存在):", pw, grp, "\n"); next
      }
      save3_base(function() {
        par(mar = c(1, 1, 2.5, 1))
        netVisual_aggregate(cc, signaling = pw,
                            signaling.name = paste(pw, "-", grp,
                                                   paste0("(", ds_label, ")")),
                            layout = "circle",
                            vertex.weight = gs, weight.scale = TRUE,
                            vertex.label.cex = 0.8,
                            color.use = setNames(NET_COLS[seq_along(gs)], names(gs)))
      }, paste0("pathway_", pw, "_", grp), out_dir)
    }
  }

  # ---- 3) 通路总强度对比条形图 + CSV ----
  # netP$prob 为 3D 数组(来源细胞 × 目标细胞 × 通路), 第 3 维为通路
  sp1 <- apply(cc1@netP$prob, 3, sum); sp2 <- apply(cc2@netP$prob, 3, sum)
  pw_all <- union(names(sp1), names(sp2))
  df <- data.frame(
    pathway = pw_all,
    g1 = as.numeric(sp1[pw_all]), g2 = as.numeric(sp2[pw_all])
  )
  df$g1[is.na(df$g1)] <- 0; df$g2[is.na(df$g2)] <- 0
  colnames(df)[2:3] <- c(g1, g2)
  df <- df[order(df[[3]] - df[[2]], decreasing = TRUE), ]
  write.csv(df, file.path(out_dir, "pathway_strength_comparison.csv"), row.names = FALSE)
  top <- df[seq_len(min(20, nrow(df))), ]
  top$pathway <- factor(top$pathway, levels = rev(top$pathway))
  dl <- data.frame(pathway = rep(top$pathway, 2),
                   strength = c(top[[2]], top[[3]]),
                   group = rep(c(g1, g2), each = nrow(top)))
  dl$group <- factor(dl$group, levels = c(g1, g2))
  p <- ggplot(dl, aes(x = pathway, y = strength, fill = group)) +
    geom_col(position = position_dodge(0.7), width = 0.65, color = "black", linewidth = 0.2) +
    coord_flip() +
    scale_fill_manual(values = setNames(c("#8FA2D6", "#D99AAE"), c(g1, g2))) +
    labs(x = NULL, y = "Total interaction strength",
         title = paste0("Signaling pathway strength: ", g1, " vs ", g2, " (", ds_label, ")")) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11),
          panel.grid.minor = element_blank(), legend.position = "top",
          legend.title = element_blank())
  save3_gg(p, "pathway_diff_bar", out_dir, w = 7.5, h = 6.5)

  # ---- 4) 差异热图(weight) ----
  tryCatch({
    ht <- netVisual_heatmap(ccm, measure = "weight",
                            color.heatmap = c("#8FA2D6", "#D99AAE"))
    save3_base(function() draw(ht), "diff_heatmap_weight", out_dir, w = 8, h = 6)
  }, error = function(e) cat("  diff_heatmap_weight 失败:", conditionMessage(e), "\n"))
}

# PBMC: top 通路 MIF / MHC-I / GALECTIN（两组合计前三）
process_dataset(file.path(BASE, "08_singlecell/cellchat"),
                "ARDS", "Sepsis", c("MIF", "MHC-I", "GALECTIN"), "GSE151263 PBMC")
# CD14+: CCL / MIF / TNF（与既有 png 一致）
process_dataset(file.path(BASE, "08_singlecell/GSE167363/cellchat"),
                "HC", "Sepsis", c("CCL", "MIF", "TNF"), "GSE167363 CD14+")

cat("\n===== Step3m 完成 =====\n")
