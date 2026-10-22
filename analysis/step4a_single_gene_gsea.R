# step4a_single_gene_gsea.R ----------------------------------------------------
# Single-gene GSEA (fgsea, MSigDB Hallmark) for the six hub genes in the two
# discovery cohorts (GSE69528 sepsis whole blood; GSE76293 ARDS neutrophils).
# Rank metric: Pearson correlation of each gene with the hub gene across samples.
# Output: full fgsea tables + Figure S8 (running-ES curves + rank metric).
# ------------------------------------------------------------------------------
suppressMessages({
  library(fgsea); library(ggplot2); library(data.table); library(patchwork); library(ragg)
})

set.seed(42)
root <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
outdir <- file.path(root, "12_single_gene_gsea")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

hub <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1")
cohorts <- c(GSE69528 = "GSE69528 (sepsis whole blood)",
             GSE76293 = "GSE76293 (ARDS neutrophils)")

gmt_file <- file.path(root, "00_rawdata/hallmark/h.all.v2024.1.Hs.symbols.gmt")
pathways <- gmtPathways(gmt_file)
cat("Hallmark sets:", length(pathways), "\n")

# Nature 7-colour palette ------------------------------------------------------
pal <- c("#6A5C95", "#D99AAE", "#E7A889", "#A9CBB8", "#CFE3E6", "#8FA2D6", "#5B608C")

run_es <- function(stats, geneset) {
  inset <- names(stats) %in% geneset
  N <- length(stats); Nh <- sum(inset)
  if (Nh == 0) return(rep(0, N))
  hit <- ifelse(inset, abs(stats) / sum(abs(stats)[inset]), 0)
  miss <- ifelse(inset, 0, -1 / (N - Nh))
  cumsum(hit + miss)
}

all_res <- list(); plots <- list(); panel_letters <- LETTERS[1:12]
pi <- 0

for (ds in names(cohorts)) {
  m <- readRDS(file.path(root, "01_DEG", paste0(ds, "_gene_expr.rds")))
  m <- as.matrix(m)
  # keep uniquely named rows, most-variance first if duplicated
  rn <- rownames(m); dup <- rn[duplicated(rn)]
  if (length(dup)) m <- m[!rn %in% dup | rn %in% names(which.max(sapply(split(rowSds <- apply(m, 1, sd), rn), max))), ]
  stopifnot(all(hub %in% rownames(m)))

  for (g in hub) {
    stats <- suppressWarnings(cor(t(m), m[g, ], method = "pearson"))
    stats <- as.numeric(stats); names(stats) <- rownames(m)
    stats <- stats[!is.na(stats)]
    stats <- sort(stats[names(stats) != g], decreasing = TRUE)

    fg <- fgsea(pathways = pathways, stats = stats, minSize = 10, maxSize = 500)
    fg <- as.data.frame(fg)
    fg$cohort <- ds; fg$hub_gene <- g; fg$rank_metric <- "Pearson r"
    all_res[[paste(ds, g)]] <- fg

    top2 <- fg[order(-abs(fg$NES)), ][1:2, ]
    cat(sprintf("%s %s | top2: %s (NES=%.2f, p=%.2g) | %s (NES=%.2f, p=%.2g)\n",
                ds, g, top2$pathway[1], top2$NES[1], top2$pval[1],
                top2$pathway[2], top2$NES[2], top2$pval[2]))

    pi <- pi + 1
    es_df <- do.call(rbind, lapply(1:2, function(k) {
      es <- run_es(stats, pathways[[top2$pathway[k]]])
      data.frame(rank = seq_along(es), ES = es,
                 pathway = top2$pathway[k], NES = top2$NES[k], col = pal[k])
    }))
    met_df <- data.frame(rank = seq_along(stats), r = stats,
                         hit = names(stats) %in% unlist(pathways[top2$pathway]))
    met_df$col <- ifelse(met_df$hit, "#C0392B", "grey70")

    neg <- top2$NES[1] < 0  # negative enrichment: curves dip, top corners free
    txt_y <- if (neg) max(es_df$ES) * 0.96 else min(es_df$ES) * 0.96
    txt_x <- max(es_df$rank) * 0.98
    p_es <- ggplot(es_df, aes(rank, ES, colour = col)) +
      geom_line(linewidth = 0.45) +
      scale_colour_identity() +
      geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.25) +
      ggtitle(paste0(LETTERS[pi], "  ", g, " \u2014 ", ds)) +
      annotate("text", x = txt_x, y = txt_y, hjust = 1,
               vjust = if (neg) 1 else 0, size = 2.1, lineheight = 0.95,
               label = paste0(gsub("HALLMARK_", "", top2$pathway),
                              " (NES ", sprintf("%.2f", top2$NES), ")",
                              collapse = "\n")) +
      labs(x = NULL, y = "Running ES") +
      theme_classic(base_size = 8) +
      theme(legend.position = "none", plot.margin = margin(2, 2, 0, 2),
            plot.title = element_text(size = 7.5, face = "bold"))

    p_met <- ggplot(met_df, aes(rank, r)) +
      geom_area(aes(fill = col), colour = NA, alpha = 0.9) +
      scale_fill_identity() +
      labs(x = "Gene rank", y = "Pearson r") +
      theme_classic(base_size = 8) +
      theme(legend.position = "none", plot.margin = margin(0, 2, 2, 2),
            axis.text = element_text(size = 6), axis.title = element_text(size = 7))

    plots[[panel_letters[pi]]] <- p_es / p_met + plot_layout(heights = c(2, 1))
  }
}

# Save full results ------------------------------------------------------------
res <- rbindlist(all_res, fill = TRUE)
res <- res[order(cohort, hub_gene, -abs(NES))]
fwrite(res, file.path(outdir, "single_gene_gsea_all_results.csv"))
top_tbl <- res[, .SD[1:2], by = .(cohort, hub_gene)]
fwrite(top_tbl, file.path(outdir, "single_gene_gsea_top2_pathways.csv"))

# Figure S8: 4 cols x 3 rows ---------------------------------------------------
fig <- wrap_plots(plots, ncol = 4) +
  plot_annotation(theme = theme(plot.margin = margin(4, 4, 4, 4)))

ggsave(file.path(outdir, "FigureS8_single_gene_gsea.pdf"), fig,
       width = 13.5, height = 10.5)
ggsave(file.path(outdir, "FigureS8_single_gene_gsea.png"), fig,
       width = 13.5, height = 10.5, dpi = 300, bg = "white")
agg_tiff(file.path(outdir, "FigureS8_single_gene_gsea.tiff"), width = 4050,
         height = 3150, units = "px", res = 300, compression = "lzw", background = "white")
print(fig); invisible(dev.off())
cat("Done. Outputs in", outdir, "\n")
