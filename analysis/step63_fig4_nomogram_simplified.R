#!/usr/bin/env Rscript
# Step 63 (subject-reviewer revision): Figure 4 nomogram panels with simplified axes.
# Reviewer request (minor comments): "图4（列线图）中的刻度线数量过多，建议简化轴标签以提高可读性。"
# Change: nomogram(..., nint = 10 default) -> nint = 4 (fewer tick marks per gene axis);
#         layout, fonts, and the Disease Probability fun.at grid are unchanged.
# Outputs: main_panels/Fig4A_nomogram_ARDS.pdf, Fig4B_nomogram_Sepsis.pdf
source("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/ajrccm_theme.R")
suppressPackageStartupMessages({ library(rms) })

BASE <- PROJ
hub <- HUB6

load_expr <- function(ds) {
  expr <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_gene_expr.rds")))
  grp  <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_grp.rds")))
  list(expr = expr, grp = grp$group)
}
d76293 <- load_expr("GSE76293")
d69528 <- load_expr("GSE69528")

build_nom <- function(train, tag, panel_nom) {
  g <- intersect(hub, rownames(train$expr))
  tr <- t(train$expr[g, , drop = FALSE]); y_tr <- as.numeric(train$grp != "Control")
  m_tr <- colMeans(tr); s_tr <- apply(tr, 2, sd)
  tr_z <- sweep(sweep(tr, 2, m_tr), 2, s_tr, "/"); colnames(tr_z) <- g
  d_tr <- data.frame(y = y_tr, tr_z)
  dd <- datadist(d_tr); assign("dd", dd, envir = .GlobalEnv); options(datadist = "dd")
  f <- lrm(y ~ ., data = d_tr, x = TRUE, y = TRUE, penalty = 1)
  nom <- nomogram(f, fun = plogis, funlabel = "Disease Probability",
                  lp = FALSE, fun.at = c(0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95),
                  nint = 4)
  grDevices::pdf(file.path(OUT_MAIN, paste0("Fig4", panel_nom, "_nomogram_", tag, ".pdf")),
                 width = 120/25.4, height = 85/25.4, family = "Helvetica", useDingbats = FALSE)
  par(cex = 0.62)
  plot(nom)
  dev.off()
  cat(sprintf("%s nomogram (simplified, nint=4) saved.\n", tag))
}

build_nom(d76293, "ARDS", "A")
build_nom(d69528, "Sepsis", "B")
cat("Step 63 done.\n")
