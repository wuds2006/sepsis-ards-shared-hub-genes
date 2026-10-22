# ============================================================
# Step 2d2: 4 个独立箱式图(Wilcoxon + 星号) + ROC 标注 AUC/CI
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2); library(pROC); library(dplyr)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "AJRCCM_submission", "07_figure_exports", "Figure2_expression_ROC")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

hub <- c("CDK5RAP2","CD177","GRB10","UPP1","EXOSC4","VNN1")
PAL7 <- c("#6A5C95","#D99AAE","#E7A889","#A9CBB8","#CFE3E6","#8FA2D6","#5B608C")

# ---- 数据读取（复用 step2d）----
load_gene_expr <- function(ds) {
  list(expr = readRDS(file.path(BASE, "01_DEG", paste0(ds, "_gene_expr.rds"))),
       grp  = readRDS(file.path(BASE, "01_DEG", paste0(ds, "_grp.rds")))$group)
}
d1 <- load_gene_expr("GSE76293"); d2 <- load_gene_expr("GSE69528")

dd <- file.path(BASE, "00_rawdata/validation/GSE243066")
files <- list.files(dd, pattern="raw_count.txt.gz$", full.names=TRUE)
mat <- NULL
for (f in files) { x <- read.delim(gzfile(f), header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
  colnames(x)[2] <- basename(f); mat <- if (is.null(mat)) x else merge(mat, x, by="Gene_symbol", all=TRUE) }
rownames(mat) <- mat$Gene_symbol; mat <- as.matrix(mat[, -1, drop=FALSE]); mat[is.na(mat)] <- 0; storage.mode(mat) <- "double"
sn <- colnames(mat); g3 <- ifelse(grepl("_A[0-9]", sn), "ARDS", ifelse(grepl("_H[0-9]", sn), "Control", "Exclude"))
cpm <- t(t(mat) / colSums(mat) * 1e6); lc <- log2(cpm + 1)
keep <- g3 %in% c("ARDS","Control"); lc <- lc[, keep, drop=FALSE]; g3 <- g3[keep]
d3 <- list(expr = lc, grp = g3)

DD <- file.path(BASE, "00_rawdata/validation/GSE100159")
lines <- readLines(file.path(DD, "GPL6884.annot")); b <- which(lines == "!platform_table_begin")
hdr <- strsplit(lines[b+1], "\t")[[1]]
ann <- read.delim(file.path(DD, "GPL6884.annot"), header=FALSE, stringsAsFactors=FALSE,
                  skip=b, quote="", comment.char="", check.names=FALSE)
colnames(ann) <- hdr; ann <- ann[ann[["Gene symbol"]] != "", c("ID","Gene symbol")]; colnames(ann) <- c("probe","gene")
sm <- readLines(file.path(DD, "GSE100159_series_matrix.txt"))
tb <- which(sm == "!series_matrix_table_begin"); te <- which(sm == "!series_matrix_table_end")
hdr2 <- gsub('"','',strsplit(sm[tb+1],"\t")[[1]])
m4 <- read.delim(textConnection(sm[(tb+2):(te-1)]), header=FALSE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
colnames(m4) <- hdr2; rownames(m4) <- m4[,1]; e4 <- as.matrix(m4[, -1, drop=FALSE]); storage.mode(e4) <- "double"
i <- which(grepl("sample class:", sm, fixed=TRUE))[1]
cls <- strsplit(gsub("!Sample_characteristics_ch1\t|\"","",sm[i]),"\t")[[1]]; cls <- cls[cls != ""][1:ncol(e4)]
g4 <- ifelse(cls %in% c("sample class: melioidosis","sample class: other infection"), "Sepsis",
      ifelse(cls == "sample class: Control", "Control", "Exclude"))
keep <- g4 %in% c("Sepsis","Control"); e4 <- e4[, keep, drop=FALSE]; g4 <- g4[keep]
common <- intersect(rownames(e4), ann$probe); e4 <- e4[common,,drop=FALSE]
gmap <- setNames(ann$gene[match(common, ann$probe)], common)
gmean <- rowMeans(e4); df <- data.frame(probe=common, gene=gmap, mean=gmean, stringsAsFactors=FALSE)
df <- df[order(df$mean, decreasing=TRUE), ]; df <- df[!duplicated(df$gene), ]
e4 <- e4[df$probe,,drop=FALSE]; rownames(e4) <- df$gene[match(df$probe, df$probe)]
d4 <- list(expr = e4, grp = g4)

datasets <- list(
  list(name="GSE76293", title="GSE76293 (ARDS discovery)", disease="ARDS", d=d1),
  list(name="GSE243066", title="GSE243066 (ARDS validation)", disease="ARDS", d=d3),
  list(name="GSE69528", title="GSE69528 (Sepsis discovery)", disease="Sepsis", d=d2),
  list(name="GSE100159", title="GSE100159 (Sepsis validation)", disease="Sepsis", d=d4)
)

# ---- 星号函数 ----
sig_star <- function(p) ifelse(p < 0.0001, "****", ifelse(p < 0.001, "***",
               ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

# ---- 4 个独立箱式图 ----
for (z in datasets) {
  g <- intersect(hub, rownames(z$d$expr))
  long <- do.call(rbind, lapply(g, function(gene){
    data.frame(gene=gene, expr=z$d$expr[gene, ], group=z$d$grp, stringsAsFactors=FALSE)
  }))
  long$gene <- factor(long$gene, levels=hub)
  long$group <- factor(long$group, levels=c("Control", z$disease))

  # Wilcoxon 检验 + 星号
  stats <- do.call(rbind, lapply(g, function(gene){
    sub <- long[long$gene == gene, ]
    x <- sub$expr[sub$group == z$disease]; y <- sub$expr[sub$group == "Control"]
    p <- wilcox.test(x, y)$p.value
    data.frame(gene=gene, p=p, star=sig_star(p),
               ypos=max(sub$expr) * 1.08, stringsAsFactors=FALSE)
  }))
  stats$gene <- factor(stats$gene, levels=hub)

  p <- ggplot(long, aes(x=gene, y=expr, fill=group)) +
    geom_boxplot(outlier.shape=NA, lwd=0.35, alpha=0.9, width=0.7, position=position_dodge(0.75)) +
    geom_jitter(aes(color=group), position=position_jitterdodge(jitter.width=0.15, dodge.width=0.75),
                size=0.5, alpha=0.55) +
    geom_text(data=stats, aes(x=gene, y=ypos, label=star), inherit.aes=FALSE,
              size=4, color="black") +
    scale_fill_manual(values=setNames(c("#8FA2D6", "#D99AAE"), c("Control", z$disease))) +
    scale_color_manual(values=setNames(c("#5B608C", "#C06C84"), c("Control", z$disease))) +
    labs(y="Expression level", x=NULL, title=z$title, fill=NULL) +
    theme_classic(base_size=10) +
    theme(legend.position="top", legend.text=element_text(size=9),
          axis.text.x=element_text(size=10, face="italic", angle=30, hjust=1),
          plot.title=element_text(size=12, face="bold", hjust=0.5),
          panel.grid=element_blank()) +
    guides(color="none")
  ggsave(file.path(OUT, paste0("Figure2_boxplot_", z$name, ".pdf")), p, width=8, height=5.5)
  ggsave(file.path(OUT, paste0("Figure2_boxplot_", z$name, ".png")), p, width=8, height=5.5, dpi=300)
  cat("箱式图", z$name, "已保存\n")
}

# ---- ROC 标注 AUC + 95% CI ----
roc_long <- do.call(rbind, lapply(datasets, function(z){
  g <- intersect(hub, rownames(z$d$expr)); yy <- as.numeric(z$d$grp != "Control")
  rr <- do.call(rbind, lapply(g, function(gene){
    ro <- roc(yy, z$d$expr[gene, ], quiet=TRUE); ci <- ci.auc(ro)
    data.frame(dataset=z$name, gene=gene, FPR=1-ro$specificities, TPR=ro$sensitivities,
               AUC=as.numeric(ro$auc), lo=ci[1], hi=ci[3])
  }))
  X <- t(z$d$expr[g, , drop=FALSE])
  if(length(g) >= 3){
    m <- glm(yy ~ X, family=binomial); sc <- predict(m, type="response")
    ro <- roc(yy, sc, quiet=TRUE); ci <- ci.auc(ro)
    rr <- rbind(rr, data.frame(dataset=z$name, gene="Combined", FPR=1-ro$specificities, TPR=ro$sensitivities,
                               AUC=as.numeric(ro$auc), lo=ci[1], hi=ci[3]))
  }
  rr
}))
roc_long$dataset <- factor(roc_long$dataset, levels=sapply(datasets, function(z) z$name))
roc_long$gene <- factor(roc_long$gene, levels=c(hub, "Combined"))
col7 <- c(setNames(PAL7[1:6], hub), Combined="#2A2A3A")

# 每个 facet 右下角标注 AUC(CI)
lab <- roc_long %>% group_by(dataset, gene) %>% summarise(AUC=first(AUC), lo=first(lo), hi=first(hi), .groups="drop") %>%
  arrange(dataset, desc(AUC)) %>%
  mutate(txt = sprintf("%s: %.3f (%.3f-%.3f)", gene, AUC, lo, hi)) %>%
  group_by(dataset) %>% summarise(txt = paste(txt, collapse="\n"), .groups="drop")
lab$dataset <- factor(lab$dataset, levels=sapply(datasets, function(z) z$name))

p_roc <- ggplot(roc_long, aes(x=FPR, y=TPR, color=gene)) +
  geom_line(lwd=0.55) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60", lwd=0.3) +
  geom_label(data=lab, aes(x=0.52, y=0.06, label=txt), inherit.aes=FALSE,
             size=2.3, hjust=0, vjust=0, linewidth=0, fill="white", alpha=0.85,
             color="black", lineheight=0.85) +
  facet_wrap(~dataset, ncol=2) +
  scale_color_manual(values=col7) +
  labs(x="1 - Specificity (FPR)", y="Sensitivity (TPR)", color=NULL) +
  theme_classic(base_size=9) +
  theme(strip.text=element_text(size=8, face="bold"),
        legend.text=element_text(size=7, face="italic"),
        legend.key.height=unit(0.55,"lines"),
        panel.grid=element_blank(), legend.position="right")
ggsave(file.path(OUT, "Figure2_ROC_AUC_CI.pdf"), p_roc, width=11, height=8.5)
ggsave(file.path(OUT, "Figure2_ROC_AUC_CI.png"), p_roc, width=11, height=8.5, dpi=300)
cat("ROC 图(标注 AUC/CI)已保存\n")

# 保存 AUC + CI 表
write.csv(roc_long %>% group_by(dataset, gene) %>% summarise(AUC=first(AUC), lo=first(lo), hi=first(hi), .groups="drop"),
          file.path(OUT, "Figure2_AUC_CI_summary.csv"), row.names=FALSE)
cat("\n完成\n")
