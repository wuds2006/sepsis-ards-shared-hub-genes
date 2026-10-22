# ============================================================
# Step 24 —— Figure S2: 11 种机器学习分类器 ROC + 混淆矩阵
# 双队列 (ARDS GSE76293 / Sepsis GSE69528), 5 折分层 CV 的
# 样本外 (out-of-fold) 预测概率绘制 ROC, 阈值 0.5 出混淆矩阵
# 与 step2a_ML.R / step2a2_ML_sepsis.R 完全同参数 (seed/特征/算法)
# ============================================================
suppressPackageStartupMessages({
  library(glmnet); library(randomForest); library(e1071); library(xgboost)
  library(rpart); library(class); library(nnet)
  library(pROC); library(ggplot2); library(patchwork)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "05_ML_SHAP")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

PAL11 <- c("#6A5C95", "#D99AAE", "#E7A889", "#A9CBB8", "#CFE3E6",
           "#8FA2D6", "#5B608C", "#8C7AA9", "#E5C8D4", "#D9BFA8", "#7A9E8E")

createFolds_safe <- function(y, disease) {
  k <- 5
  idx_ctl <- which(y == "Control"); idx_dis <- which(y == disease)
  set.seed(1); idx_ctl <- sample(idx_ctl); idx_dis <- sample(idx_dis)
  f <- vector("list", k)
  for (i in 1:k) {
    f[[i]] <- c(idx_ctl[seq(i, length(idx_ctl), k)], idx_dis[seq(i, length(idx_dis), k)])
  }
  f
}

# 11 分类器 (与 step2a 完全一致, GBM 已按口径排除)
make_algos <- function(pos) list(
  LASSO = function(trX, trY, teX) {
    cv <- cv.glmnet(trX, as.numeric(trY==pos), family="binomial", alpha=1)
    as.numeric(predict(cv$glmnet.fit, teX, s=cv$lambda.min, type="response")) },
  Ridge = function(trX, trY, teX) {
    cv <- cv.glmnet(trX, as.numeric(trY==pos), family="binomial", alpha=0)
    as.numeric(predict(cv$glmnet.fit, teX, s=cv$lambda.min, type="response")) },
  ENet  = function(trX, trY, teX) {
    cv <- cv.glmnet(trX, as.numeric(trY==pos), family="binomial", alpha=0.5)
    as.numeric(predict(cv$glmnet.fit, teX, s=cv$lambda.min, type="response")) },
  RF    = function(trX, trY, teX) {
    m <- randomForest(trX, trY, ntree=500); as.numeric(predict(m, teX, type="prob")[,2]) },
  SVM   = function(trX, trY, teX) {
    m <- svm(trX, trY, probability=TRUE, kernel="radial")
    attr(predict(m, teX, probability=TRUE), "probabilities")[, pos] },
  KNN   = function(trX, trY, teX) {
    pr <- knn(trX, teX, trY, k=3, prob=TRUE)
    p <- attr(pr, "prob"); p[pr == levels(trY)[1]] <- 1 - p[pr == levels(trY)[1]]; p },
  NB    = function(trX, trY, teX) {
    m <- naiveBayes(trX, trY); predict(m, teX, type="raw")[, pos] },
  NN    = function(trX, trY, teX) {
    m <- nnet(y ~ ., data=data.frame(y=trY, trX), size=3, trace=FALSE, maxit=200, decay=0.1)
    as.numeric(predict(m, data.frame(teX), type="raw")) },
  LR    = function(trX, trY, teX) {
    d <- data.frame(y = trY, trX); m <- glm(y ~ ., data = d, family = binomial)
    predict(m, data.frame(teX), type = "response") },
  XGB   = function(trX, trY, teX) {
    dtrain <- xgb.DMatrix(trX, label = as.numeric(trY==pos))
    m <- xgb.train(params = list(objective="binary:logistic", eta=0.1, max_depth=3),
                   data = dtrain, nrounds = 100, verbose = 0)
    predict(m, xgb.DMatrix(teX)) },
  DT    = function(trX, trY, teX) {
    m <- rpart(y ~ ., data=data.frame(y=trY, trX), method="class")
    predict(m, data.frame(teX), type="prob")[, pos] }
)

run_cohort <- function(ds, disease) {
  cat("\n##########", ds, disease, "##########\n")
  expr <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_gene_expr.rds")))
  grp  <- readRDS(file.path(BASE, "01_DEG", paste0(ds, "_grp.rds")))
  cgs  <- readLines(file.path(BASE, "03_intersection", "Shared_CGs.txt"))
  cgs  <- intersect(cgs, rownames(expr))
  X <- scale(t(expr[cgs, , drop = FALSE]))
  y <- factor(grp$group, levels = c("Control", disease))
  cat(sprintf("特征: %d, 样本: %d (%s=%d, Control=%d)\n",
              ncol(X), nrow(X), disease, sum(y == disease), sum(y == "Control")))

  folds <- createFolds_safe(y, disease)
  algos <- make_algos(disease)

  oof <- matrix(NA_real_, nrow(X), length(algos),
                dimnames = list(rownames(X), names(algos)))
  for (nm in names(algos)) {
    for (i in seq_along(folds)) {
      tr <- unlist(folds[-i]); te <- folds[[i]]
      prob <- tryCatch(algos[[nm]](X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE]),
                       error = function(e) rep(NA_real_, length(te)))
      oof[te, nm] <- prob
    }
    cat(nm, "OOF AUC =",
        sprintf("%.3f", auc(roc(y, oof[, nm], levels = c("Control", disease),
                                direction = "<", quiet = TRUE))), "\n")
  }

  # ---- ROC 曲线 ----
  plist <- lapply(names(algos), function(nm) {
    ro <- roc(y, oof[, nm], levels = c("Control", disease), direction = "<", quiet = TRUE)
    data.frame(fpr = 1 - ro$specificities, tpr = ro$sensitivities,
               model = nm, auc = as.numeric(ro$auc))
  })
  df <- do.call(rbind, plist)
  df$model <- factor(df$model, levels = names(algos))
  aucs <- sapply(split(df, df$model), function(d) unique(d$auc))
  lab <- paste0(names(algos), " (AUC = ", sprintf("%.3f", aucs[names(algos)]), ")")

  p_roc <- ggplot(df, aes(x = fpr, y = tpr, color = model)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "#999999", linewidth = 0.4) +
    geom_path(linewidth = 0.55, key_glyph = "path") +
    scale_color_manual(values = setNames(PAL11, names(algos)),
                       labels = setNames(lab, names(algos))) +
    scale_x_continuous(limits = c(0, 1), expand = c(0.002, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0.002, 0)) +
    labs(x = "1 − Specificity", y = "Sensitivity",
         title = paste0(ds, " (", disease, "): ROC of 11 classifiers"),
         subtitle = "Pooled out-of-fold predictions, 5-fold stratified CV") +
    theme_pub() +
    theme(plot.title = element_text(size = 9, face = "bold"),
          plot.subtitle = element_text(size = 6.5, color = "#666666"),
          legend.position = c(0.62, 0.18),
          legend.title = element_blank(),
          legend.text = element_text(size = 5.8),
          legend.key.size = unit(0.35, "cm"),
          legend.background = element_rect(fill = "white", color = "#dddddd", linewidth = 0.3),
          axis.title = element_text(size = 7.5), axis.text = element_text(size = 6.5))

  # ---- 混淆矩阵 ----
  cm_list <- lapply(names(algos), function(nm) {
    pred <- factor(ifelse(oof[, nm] >= 0.5, disease, "Control"), levels = c(disease, "Control"))
    tb <- table(Reference = y, Predicted = pred)
    data.frame(model = nm,
               Reference = rep(c(disease, "Control"), each = 2),
               Predicted = rep(c(disease, "Control"), 2),
               n = as.vector(tb))
  })
  cm <- do.call(rbind, cm_list)
  cm$model <- factor(cm$model, levels = names(algos))
  cm$lab <- paste0(cm$Reference, "\n", cm$Predicted)

  p_cm <- ggplot(cm, aes(x = Predicted, y = Reference, fill = n)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = n), size = 2.7, fontface = "bold") +
    facet_wrap(~ model, nrow = 2) +
    scale_fill_gradient(low = "#CFE3E6", high = "#5B608C", name = "Count") +
    scale_x_discrete(position = "bottom") +
    labs(x = "Predicted", y = "Reference",
         title = paste0(ds, ": confusion matrices (threshold 0.5)")) +
    theme_pub() +
    theme(plot.title = element_text(size = 9, face = "bold"),
          strip.text = element_text(size = 6.5, face = "bold"),
          strip.background = element_rect(fill = "#F0F0F4", color = NA),
          axis.title = element_text(size = 7.5),
          axis.text = element_text(size = 6),
          panel.spacing = unit(0.12, "cm"),
          legend.position = "right",
          legend.key.height = unit(0.6, "cm"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6.5))

  write.csv(data.frame(sample = rownames(oof), oof, group = y),
            file.path(OUT, paste0(ds, "_OOF_predictions.csv")), row.names = FALSE)
  save_fig(p_roc, paste0("figS2_ROC_", ds), w = 5.2, h = 4.6, dir = OUT)
  save_fig(p_cm,  paste0("figS2_CM_", ds),  w = 9.2, h = 4.4, dir = OUT)
  list(roc = p_roc, cm = p_cm)
}

r1 <- run_cohort("GSE76293", "ARDS")
r2 <- run_cohort("GSE69528", "Sepsis")

# ---- composite: A/B ROC, C/D CM ----
comp <- (r1$roc + r2$roc) / (r1$cm + r2$cm) +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, face = "bold", color = "#2b2b2b"))

save_fig(comp, "FigureS2_ML_ROC_confusion", w = 14.6, h = 9.2, dir = OUT)
cat("\n===== Step 24 Figure S2 完成 =====\n")
cat("OK\n")
