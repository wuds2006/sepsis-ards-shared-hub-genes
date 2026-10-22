# ============================================================
# Step 2d: 6 核心基因表达模式(箱式图) + ROC/AUC (4 数据集)
# 实验组: GSE76293(ARDS) + GSE69528(Sepsis)
# 验证组: GSE243066(ARDS) + GSE100159(Sepsis)
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2); library(pROC); library(patchwork); library(dplyr); library(scales)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "AJRCCM_submission", "07_figure_exports", "Figure2_expression_ROC")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

hub <- c("CDK5RAP2","CD177","GRB10","UPP1","EXOSC4","VNN1")
PAL7 <- c("#6A5C95","#D99AAE","#E7A889","#A9CBB8","#CFE3E6","#8FA2D6","#5B608C")

load_gene_expr <- function(ds) {
  expr <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_gene_expr.rds")))
  grp  <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_grp.rds")))
  list(expr = expr, grp = grp)
}

# ---- 1. GSE76293 (ARDS 实验组) ----
d1 <- load_gene_expr("GSE76293")

# ---- 2. GSE69528 (Sepsis 实验组) ----
d2 <- load_gene_expr("GSE69528")

# ---- 3. GSE243066 (ARDS 验证组, RNA-seq) ----
dd <- file.path(BASE, "00_rawdata/validation/GSE243066")
files <- list.files(dd, pattern="raw_count.txt.gz$", full.names=TRUE)
mat <- NULL
for (f in files) {
  x <- read.delim(gzfile(f), header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
  colnames(x)[2] <- basename(f)
  mat <- if (is.null(mat)) x else merge(mat, x, by="Gene_symbol", all=TRUE)
}
rownames(mat) <- mat$Gene_symbol; mat <- as.matrix(mat[, -1, drop=FALSE])
mat[is.na(mat)] <- 0; storage.mode(mat) <- "double"
sn <- colnames(mat)
g3 <- ifelse(grepl("_A[0-9]", sn), "ARDS", ifelse(grepl("_H[0-9]", sn), "Control", "Exclude"))
cpm <- t(t(mat) / colSums(mat) * 1e6); logcpm <- log2(cpm + 1)
keep <- g3 %in% c("ARDS","Control"); logcpm <- logcpm[, keep, drop=FALSE]; g3 <- g3[keep]
d3 <- list(expr = logcpm, grp = data.frame(group = g3))

# ---- 4. GSE100159 (Sepsis 验证组, GPL6884) ----
DD <- file.path(BASE, "00_rawdata/validation/GSE100159")
lines <- readLines(file.path(DD, "GPL6884.annot"))
b <- which(lines == "!platform_table_begin")
hdr <- strsplit(lines[b+1], "\t")[[1]]
ann <- read.delim(file.path(DD, "GPL6884.annot"), header=FALSE, stringsAsFactors=FALSE,
                  skip=b, quote="", comment.char="", check.names=FALSE)
colnames(ann) <- hdr
ann <- ann[ann[["Gene symbol"]] != "", c("ID","Gene symbol")]
colnames(ann) <- c("probe","gene")
sm <- readLines(file.path(DD, "GSE100159_series_matrix.txt"))
tb <- which(sm == "!series_matrix_table_begin"); te <- which(sm == "!series_matrix_table_end")
hdr2 <- gsub('"','',strsplit(sm[tb+1],"\t")[[1]])
m4 <- read.delim(textConnection(sm[(tb+2):(te-1)]), header=FALSE, sep="\t",
                 stringsAsFactors=FALSE, check.names=FALSE)
colnames(m4) <- hdr2; rownames(m4) <- m4[,1]
e4 <- as.matrix(m4[, -1, drop=FALSE]); storage.mode(e4) <- "double"
i <- which(grepl("sample class:", sm, fixed=TRUE))[1]
cls <- strsplit(gsub("!Sample_characteristics_ch1\t|\"","",sm[i]),"\t")[[1]]
cls <- cls[cls != ""][1:ncol(e4)]
g4 <- ifelse(cls %in% c("sample class: melioidosis","sample class: other infection"), "Sepsis",
      ifelse(cls == "sample class: Control", "Control", "Exclude"))
keep <- g4 %in% c("Sepsis","Control"); e4 <- e4[, keep, drop=FALSE]; g4 <- g4[keep]
common <- intersect(rownames(e4), ann$probe); e4 <- e4[common,,drop=FALSE]
gmap <- setNames(ann$gene[match(common, ann$probe)], common)
gmean <- rowMeans(e4)
df <- data.frame(probe=common, gene=gmap, mean=gmean, stringsAsFactors=FALSE)
df <- df[order(df$mean, decreasing=TRUE), ]; df <- df[!duplicated(df$gene), ]
e4 <- e4[df$probe,,drop=FALSE]; rownames(e4) <- df$gene[match(df$probe, df$probe)]
d4 <- list(expr = e4, grp = data.frame(group = g4))

# ---- 汇总为长格式(箱式图用) ----
datasets <- list(
  list(name="GSE76293\nARDS discovery", disease="ARDS", d=d1),
  list(name="GSE243066\nARDS validation", disease="ARDS", d=d3),
  list(name="GSE69528\nSepsis discovery", disease="Sepsis", d=d2),
  list(name="GSE100159\nSepsis validation", disease="Sepsis", d=d4)
)
long <- do.call(rbind, lapply(datasets, function(z){
  g <- intersect(hub, rownames(z$d$expr))
  do.call(rbind, lapply(g, function(gene){
    data.frame(dataset=z$name, gene=gene,
               expr=z$d$expr[gene, ],
               group=z$d$grp$group, stringsAsFactors=FALSE)
  }))
}))
long$gene <- factor(long$gene, levels=hub)
long$group <- factor(long$group, levels=c("Control","ARDS","Sepsis"))
long$dataset <- factor(long$dataset, levels=sapply(datasets, function(z) z$name))

# ---- 箱式图: 6 基因 x 4 数据集 ----
p_box <- ggplot(long, aes(x=group, y=expr, fill=group)) +
  geom_boxplot(outlier.shape=NA, lwd=0.3, alpha=0.85) +
  geom_jitter(width=0.18, size=0.35, alpha=0.5, color="grey25") +
  facet_grid(dataset ~ gene, scales="free_y") +
  scale_fill_manual(values=c(Control="#8FA2D6", ARDS="#D99AAE", Sepsis="#D99AAE")) +
  labs(y="Expression level", x=NULL) +
  theme_classic(base_size=8) +
  theme(legend.position="none",
        strip.text.y=element_text(size=7, face="bold"),
        strip.text.x=element_text(size=7, face="italic"),
        axis.text.x=element_text(size=6.5, angle=30, hjust=1),
        axis.text.y=element_text(size=6),
        panel.grid=element_blank())
ggsave(file.path(OUT, "Figure2_expression_boxplot.pdf"), p_box, width=13, height=9)
ggsave(file.path(OUT, "Figure2_expression_boxplot.png"), p_box, width=13, height=9, dpi=300)
cat("箱式图已保存\n")

# ---- ROC: 4 数据集 x (6 基因 + 联合签名) ----
roc_long <- do.call(rbind, lapply(datasets, function(z){
  g <- intersect(hub, rownames(z$d$expr))
  yy <- as.numeric(z$d$grp$group != "Control")
  # 单基因
  rr <- do.call(rbind, lapply(g, function(gene){
    ro <- roc(yy, z$d$expr[gene, ], quiet=TRUE)
    data.frame(dataset=z$name, gene=gene, FPR=1-ro$specificities, TPR=ro$sensitivities, AUC=as.numeric(ro$auc))
  }))
  # 联合签名
  X <- t(z$d$expr[g, , drop=FALSE])
  if(length(g) >= 3){
    m <- glm(yy ~ X, family=binomial)
    sc <- predict(m, type="response")
    ro <- roc(yy, sc, quiet=TRUE)
    rr <- rbind(rr, data.frame(dataset=z$name, gene="Combined", FPR=1-ro$specificities, TPR=ro$sensitivities, AUC=as.numeric(ro$auc)))
  }
  rr
}))
roc_long$dataset <- factor(roc_long$dataset, levels=sapply(datasets, function(z) z$name))
roc_long$gene <- factor(roc_long$gene, levels=c(hub, "Combined"))
col7 <- c(setNames(PAL7[1:6], hub), Combined="#2A2A3A")

p_roc <- ggplot(roc_long, aes(x=FPR, y=TPR, color=gene)) +
  geom_line(lwd=0.5) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60", lwd=0.3) +
  facet_wrap(~dataset, ncol=2) +
  scale_color_manual(values=col7) +
  labs(x="1 - Specificity (FPR)", y="Sensitivity (TPR)", color=NULL) +
  theme_classic(base_size=8) +
  theme(strip.text=element_text(size=7, face="bold"),
        legend.text=element_text(size=6.5, face="italic"),
        legend.key.height=unit(0.55,"lines"),
        panel.grid=element_blank(),
        legend.position="right")
ggsave(file.path(OUT, "Figure2_ROC.pdf"), p_roc, width=10, height=8)
ggsave(file.path(OUT, "Figure2_ROC.png"), p_roc, width=10, height=8, dpi=300)
cat("ROC 图已保存\n")

# ---- 输出 AUC 汇总表 ----
auc_sum <- roc_long %>% group_by(dataset, gene) %>% summarise(AUC=mean(AUC), .groups="drop") %>%
  arrange(dataset, desc(AUC))
write.csv(auc_sum, file.path(OUT, "Figure2_AUC_summary.csv"), row.names=FALSE)
cat("\n===== AUC 汇总 =====\n")
print(as.data.frame(auc_sum), row.names=FALSE)
cat("\n完成\n")
