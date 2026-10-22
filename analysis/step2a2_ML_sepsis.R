# ============================================================
# Step 2a2 —— Sepsis 队列: 11 种机器学习算法 + SHAP 解释
# 用共享 CGs(179) 预测 Sepsis 发生(GSE69528), 5折CV评估AUC, SHAP识别关键基因
# 与 step2a_ML.R (ARDS, GSE76293) 平行; 最后取两队列候选枢纽基因交集
# 算法: LR/LASSO/Ridge/ENet/RF/SVM/XGBoost/GBM/DT/KNN/NB/NN
# ============================================================
suppressPackageStartupMessages({
  library(glmnet); library(randomForest); library(e1071); library(xgboost)
  library(gbm); library(rpart); library(class); library(nnet)
  library(pROC); library(ggplot2); library(kernelshap); library(shapviz)
})
source("scripts/common.R")

BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "05_ML_SHAP", "sepsis")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

# ---- 数据准备 ----
expr <- readRDS(file.path(BASE, "01_DEG", "GSE69528_gene_expr.rds"))
grp  <- readRDS(file.path(BASE, "01_DEG", "GSE69528_grp.rds"))
cgs  <- readLines(file.path(BASE, "03_intersection", "Shared_CGs.txt"))
cgs  <- intersect(cgs, rownames(expr))
X <- t(expr[cgs, , drop = FALSE])          # 样本 x 基因
y <- factor(grp$group, levels = c("Control", "Sepsis"))
cat(sprintf("特征(CGs): %d, 样本: %d (Sepsis=%d, Control=%d)\n",
            length(cgs), nrow(X), sum(y == "Sepsis"), sum(y == "Control")))

X <- scale(X)

# ---- 11 算法 5折CV ----
createFolds_safe <- function(y) {
  k <- 5
  idx_ctl <- which(y == "Control"); idx_dis <- which(y == "Sepsis")
  set.seed(1); idx_ctl <- sample(idx_ctl); idx_dis <- sample(idx_dis)
  f <- vector("list", k)
  for (i in 1:k) {
    f[[i]] <- c(idx_ctl[seq(i, length(idx_ctl), k)], idx_dis[seq(i, length(idx_dis), k)])
  }
  f
}
folds <- createFolds_safe(y)
run_cv <- function(fit_predict, name) {
  aucs <- c()
  for (i in seq_along(folds)) {
    tr <- unlist(folds[-i]); te <- folds[[i]]
    prob <- tryCatch(fit_predict(X[tr, , drop = F], y[tr], X[te, , drop = F]),
                     error = function(e) rep(NA, length(te)))
    if (all(is.na(prob))) { aucs <- c(aucs, NA); next }
    ro <- tryCatch(roc(y[te], prob, levels = c("Control", "Sepsis"), direction = "<", quiet = TRUE)$auc,
                   error = function(e) NA)
    aucs <- c(aucs, ro)
  }
  c(model = name, mean_AUC = mean(aucs, na.rm = TRUE), sd_AUC = sd(aucs, na.rm = TRUE))
}

algos <- list(
  LR   = function(trX, trY, teX) {
    d <- data.frame(y = trY, trX); m <- glm(y ~ ., data = d, family = binomial)
    predict(m, data.frame(teX), type = "response") },
  LASSO = function(trX, trY, teX) {
    cv <- cv.glmnet(trX, as.numeric(trY == "Sepsis"), family = "binomial", alpha = 1)
    as.numeric(predict(cv$glmnet.fit, teX, s = cv$lambda.min, type = "response")) },
  Ridge = function(trX, trY, teX) {
    cv <- cv.glmnet(trX, as.numeric(trY == "Sepsis"), family = "binomial", alpha = 0)
    as.numeric(predict(cv$glmnet.fit, teX, s = cv$lambda.min, type = "response")) },
  ENet  = function(trX, trY, teX) {
    cv <- cv.glmnet(trX, as.numeric(trY == "Sepsis"), family = "binomial", alpha = 0.5)
    as.numeric(predict(cv$glmnet.fit, teX, s = cv$lambda.min, type = "response")) },
  RF    = function(trX, trY, teX) {
    m <- randomForest(trX, trY, ntree = 500); as.numeric(predict(m, teX, type = "prob")[, 2]) },
  SVM   = function(trX, trY, teX) {
    m <- svm(trX, trY, probability = TRUE, kernel = "radial")
    attr(predict(m, teX, probability = TRUE), "probabilities")[, "Sepsis"] },
  XGB   = function(trX, trY, teX) {
    dtrain <- xgb.DMatrix(trX, label = as.numeric(trY == "Sepsis"))
    m <- xgb.train(params = list(objective = "binary:logistic", eta = 0.1, max_depth = 3),
                   data = dtrain, nrounds = 100, verbose = 0)
    predict(m, xgb.DMatrix(teX)) },
  GBM   = function(trX, trY, teX) {
    d <- data.frame(y = as.numeric(trY == "Sepsis"), trX)
    m <- gbm(y ~ ., data = d, distribution = "bernoulli", n.trees = 100,
             interaction.depth = 2, shrinkage = 0.05, verbose = FALSE, n.minobsinnode = 3)
    predict(m, newdata = as.data.frame(teX), n.trees = 100, type = "response") },
  DT    = function(trX, trY, teX) {
    m <- rpart(y ~ ., data = data.frame(y = trY, trX), method = "class")
    predict(m, data.frame(teX), type = "prob")[, "Sepsis"] },
  KNN   = function(trX, trY, teX) {
    pr <- knn(trX, teX, trY, k = 3, prob = TRUE)
    p <- attr(pr, "prob"); p[pr == "Control"] <- 1 - p[pr == "Control"]; p },
  NB    = function(trX, trY, teX) {
    m <- naiveBayes(trX, trY); predict(m, teX, type = "raw")[, "Sepsis"] },
  NN    = function(trX, trY, teX) {
    m <- nnet(y ~ ., data = data.frame(y = trY, trX), size = 3, trace = FALSE, maxit = 200, decay = 0.1)
    as.numeric(predict(m, data.frame(teX), type = "raw")) }
)

res <- do.call(rbind, lapply(names(algos), function(nm) run_cv(algos[[nm]], nm)))
res <- as.data.frame(res, stringsAsFactors = FALSE)
res$mean_AUC <- as.numeric(res$mean_AUC); res$sd_AUC <- as.numeric(res$sd_AUC)
res <- res[order(res$mean_AUC, decreasing = TRUE), ]
print(res)
write.csv(res, file.path(OUT, "ML_AUC_results_sepsis.csv"), row.names = FALSE)

p_auc <- ggplot(res, aes(x = reorder(model, mean_AUC), y = mean_AUC, fill = mean_AUC)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean_AUC - sd_AUC, ymax = mean_AUC + sd_AUC), width = 0.2) +
  geom_text(aes(label = sprintf("%.3f", mean_AUC)), hjust = -0.2, size = 3) +
  coord_flip(ylim = c(0, 1.15)) +
  scale_fill_gradient(low = PAL["periwinkle"], high = PAL["purple"]) +
  labs(x = "Algorithm", y = "Mean AUC (5-fold CV)", title = "11 ML Algorithms: Sepsis Prediction") +
  theme_pub() + theme(legend.position = "none")
save_fig(p_auc, "ML_AUC_barplot_sepsis", w = 7, h = 5, dir = OUT)

# ---- 特征重要性(全样本拟合) + SHAP ----
rf_final <- randomForest(X, y, ntree = 1000, importance = TRUE)
imp_rf <- importance(rf_final)[, "MeanDecreaseGini"]
imp_rf <- sort(imp_rf, decreasing = TRUE)

dtrain_full <- xgb.DMatrix(X, label = as.numeric(y == "Sepsis"))
xgb_final <- xgb.train(params = list(objective = "binary:logistic", eta = 0.1, max_depth = 3),
                       data = dtrain_full, nrounds = 100, verbose = 0)
imp_xgb <- xgb.importance(model = xgb_final, feature_names = colnames(X))
xgb_gain <- setNames(imp_xgb$Gain, imp_xgb$Feature)
xgb_gain <- xgb_gain[colnames(X)]; xgb_gain[is.na(xgb_gain)] <- 0

cv_lasso <- cv.glmnet(X, as.numeric(y == "Sepsis"), family = "binomial", alpha = 1)
coef_lasso <- as.matrix(coef(cv_lasso, s = "lambda.min"))[-1, , drop = FALSE]
lasso_genes <- rownames(coef_lasso)[coef_lasso != 0]
cat("LASSO 选中基因数:", length(lasso_genes), "\n")

top_genes <- unique(c(lasso_genes, names(imp_rf)[1:30]))
top_genes <- head(top_genes, 30)
cat("SHAP 使用 top", length(top_genes), "基因\n")

X_top <- X[, top_genes, drop = FALSE]
rf_top <- randomForest(X_top, y, ntree = 1000, importance = TRUE)
rf_pred <- function(model, newdata) predict(model, newdata, type = "prob")[, 2]
bg <- X_top[1:min(20, nrow(X_top)), , drop = FALSE]
ks <- kernelshap(rf_top, X = X_top, bg_X = bg, pred_fun = rf_pred)
sv <- shapviz(ks)
p_shap <- sv_importance(sv, kind = "bee", max_display = 20) +
  scale_color_gradient2(low = PAL_DOWN, mid = "white", high = PAL_UP) +
  ggtitle("SHAP Feature Importance (Random Forest, Sepsis)") + theme_pub()
save_fig(p_shap, "SHAP_beeplot_sepsis", w = 8, h = 7, dir = OUT)

shap_imp <- data.frame(gene = top_genes,
                       mean_abs_shap = colMeans(abs(as.data.frame(ks$S))),
                       stringsAsFactors = FALSE)
shap_imp <- shap_imp[order(shap_imp$mean_abs_shap, decreasing = TRUE), ]
write.csv(shap_imp, file.path(OUT, "SHAP_importance_sepsis.csv"), row.names = FALSE)

# ---- Sepsis 候选枢纽基因: LASSO ∩ (SHAP top 或 RF top) ----
rf_genes <- names(imp_rf)[1:30]
shap_genes <- shap_imp$gene[1:20]
key_genes <- lasso_genes[lasso_genes %in% union(rf_genes, shap_genes)]
if (length(key_genes) < 5) key_genes <- shap_imp$gene[1:min(15, nrow(shap_imp))]
writeLines(key_genes, file.path(OUT, "ML_key_genes_sepsis.txt"))
cat("Sepsis 候选枢纽基因:", paste(key_genes, collapse = ", "), "\n")

write.csv(data.frame(gene = names(imp_rf), RF_gini = imp_rf), file.path(OUT, "RF_importance_sepsis.csv"), row.names = FALSE)
write.csv(data.frame(gene = names(xgb_gain), XGB_gain = xgb_gain), file.path(OUT, "XGB_importance_sepsis.csv"), row.names = FALSE)

# ---- 与 ARDS 候选枢纽基因取交集 → 核心基因 ----
ards_genes <- readLines(file.path(BASE, "05_ML_SHAP", "ML_key_genes.txt"))
shared <- intersect(ards_genes, key_genes)
writeLines(shared, file.path(OUT, "ML_shared_hub_genes.txt"))
cat("\nARDS 候选枢纽基因 (", length(ards_genes), "):", paste(ards_genes, collapse = ", "), "\n")
cat("交集核心基因 (", length(shared), "):", paste(shared, collapse = ", "), "\n")
cat("\n===== Step2a2 完成 =====\n")
