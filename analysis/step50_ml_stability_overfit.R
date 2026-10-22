# ============================================================
# Step 50 —— ML 过拟合量化 + 特征选择稳定性（审稿 W1/W2）
# B1: permutation test (500x) / apparent-vs-CV AUC / AUC=1.00 exact lower bounds
# B2: LASSO bootstrap stability / RF 50-seed inclusion / threshold perturbation Jaccard
# 镜像 step2a_ML.R 与 step2a2_ML_sepsis.R 的数据口径与算法定义
# ============================================================
suppressPackageStartupMessages({
  library(glmnet); library(randomForest); library(e1071); library(xgboost)
  library(gbm); library(rpart); library(class); library(nnet)
  library(pROC)
})
BASE <- "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT  <- file.path(BASE, "05_ML_SHAP", "stability")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cgs <- readLines(file.path(BASE, "03_intersection", "Shared_CGs.txt"))

make_algos <- function(pos) list(
  LR    = function(trX, trY, teX) {
    d <- data.frame(y = trY, trX); m <- glm(y ~ ., data = d, family = binomial)
    predict(m, data.frame(teX), type = "response") },
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
    attr(predict(m, teX, probability=TRUE), "probabilities")[,pos] },
  XGB   = function(trX, trY, teX) {
    dtrain <- xgb.DMatrix(trX, label = as.numeric(trY==pos))
    m <- xgb.train(params = list(objective="binary:logistic", eta=0.1, max_depth=3),
                   data = dtrain, nrounds = 100, verbose = 0)
    predict(m, xgb.DMatrix(teX)) },
  GBM   = function(trX, trY, teX) {
    d <- data.frame(y = as.numeric(trY==pos), trX)
    m <- gbm(y ~ ., data = d, distribution = "bernoulli", n.trees = 100,
             interaction.depth = 2, shrinkage = 0.05, verbose = FALSE, n.minobsinnode = 3)
    predict(m, newdata = as.data.frame(teX), n.trees = 100, type = "response") },
  DT    = function(trX, trY, teX) {
    m <- rpart(y ~ ., data=data.frame(y=trY, trX), method="class")
    predict(m, data.frame(teX), type="prob")[,pos] },
  KNN   = function(trX, trY, teX) {
    pr <- knn(trX, teX, trY, k=3, prob=TRUE)
    p <- attr(pr, "prob"); p[pr=="Control"] <- 1-p[pr=="Control"]; p },
  NB    = function(trX, trY, teX) {
    m <- naiveBayes(trX, trY); predict(m, teX, type="raw")[,pos] },
  NN    = function(trX, trY, teX) {
    m <- nnet(y ~ ., data=data.frame(y=trY, trX), size=3, trace=FALSE, maxit=200, decay=0.1)
    as.numeric(predict(m, data.frame(teX), type="raw")) }
)

createFolds_safe <- function(y, k = 5) {
  idx_ctl <- which(y == "Control"); idx_dis <- which(y != "Control")
  set.seed(1); idx_ctl <- sample(idx_ctl); idx_dis <- sample(idx_dis)
  f <- vector("list", k)
  for (i in 1:k) f[[i]] <- c(idx_ctl[seq(i, length(idx_ctl), k)], idx_dis[seq(i, length(idx_dis), k)])
  f
}

cv_auc <- function(algos, X, y, folds, lvl) {
  sapply(names(algos), function(nm) {
    aucs <- c()
    for (i in seq_along(folds)) {
      tr <- unlist(folds[-i]); te <- folds[[i]]
      prob <- tryCatch(algos[[nm]](X[tr,,drop=F], y[tr], X[te,,drop=F]), error=function(e) rep(NA, length(te)))
      ro  <- tryCatch(roc(y[te], prob, levels = c("Control", lvl), direction="<", quiet=TRUE)$auc, error=function(e) NA)
      aucs <- c(aucs, ro)
    }
    mean(aucs, na.rm = TRUE)
  })
}

run_cohort <- function(tag, gse, lvl) {
  expr <- readRDS(file.path(BASE, "01_DEG", paste0(gse, "_gene_expr.rds")))
  grp  <- readRDS(file.path(BASE, "01_DEG", paste0(gse, "_grp.rds")))
  ccgs <- intersect(cgs, rownames(expr))
  X <- scale(t(expr[ccgs, , drop = FALSE]))
  y <- factor(grp$group, levels = c("Control", lvl))
  algos <- make_algos(lvl)
  folds <- createFolds_safe(y)
  cat(sprintf("\n===== %s (%s): n=%d (%s=%d, Control=%d) =====\n",
              tag, gse, nrow(X), lvl, sum(y==lvl), sum(y=="Control")))

  # ---- B1(2) real CV AUC + apparent AUC ----
  real_cv <- cv_auc(algos, X, y, folds, lvl)
  apparent <- sapply(names(algos), function(nm) {
    prob <- tryCatch(algos[[nm]](X, y, X), error=function(e) rep(NA, nrow(X)))
    tryCatch(roc(y, prob, levels=c("Control", lvl), direction="<", quiet=TRUE)$auc, error=function(e) NA)
  })
  tab <- data.frame(model = names(real_cv), cv_AUC = round(real_cv, 3),
                    apparent_AUC = round(apparent, 3),
                    gap = round(apparent - real_cv, 3))
  write.csv(tab, file.path(OUT, sprintf("apparent_vs_cv_AUC_%s.csv", tag)), row.names = FALSE)
  print(tab)

  # ---- B1(1) permutation test 500x (仅 ARDS 小队列执行) ----
  if (tag == "ARDS") {
    B <- 500
    set.seed(2026)
    null_auc <- matrix(NA_real_, nrow = B, ncol = length(algos),
                       dimnames = list(NULL, names(algos)))
    t0 <- Sys.time()
    for (b in 1:B) {
      yb <- factor(sample(as.character(y)), levels = c("Control", lvl))
      # 标签打乱后分层折按原样本索引结构固定（与原 folds 一致）
      null_auc[b, ] <- tryCatch(cv_auc(algos, X, yb, folds, lvl),
                                error = function(e) rep(NA_real_, length(algos)))
      if (b %% 25 == 0) {
        el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
        cat(sprintf("  perm %d/%d (%.1f min elapsed, ETA %.1f min)\n",
                    b, B, el, el / b * (B - b)))
      }
    }
    perm <- data.frame(
      model = names(algos),
      real_cv_AUC = round(real_cv, 3),
      null_mean = round(colMeans(null_auc, na.rm=TRUE), 3),
      null_sd   = round(apply(null_auc, 2, sd, na.rm=TRUE), 3),
      null_max  = round(apply(null_auc, 2, max, na.rm=TRUE), 3))
    perm$perm_P <- sapply(names(algos), function(nm) {
      (1 + sum(null_auc[, nm] >= real_cv[[nm]], na.rm=TRUE)) / (1 + B)
    })
    write.csv(perm, file.path(OUT, "permutation_500x_results.csv"), row.names = FALSE)
    print(perm)
    saveRDS(null_auc, file.path(OUT, "permutation_null_dist.rds"))
  }

  # ---- B2(1) LASSO bootstrap 稳定性 (200x) ----
  set.seed(7)
  nboot <- 200
  sel_freq <- setNames(numeric(length(ccgs)), ccgs)
  for (b in 1:nboot) {
    idx <- sample(nrow(X), nrow(X), replace = TRUE)
    cv <- tryCatch(cv.glmnet(X[idx,,drop=F], as.numeric(y[idx]==lvl),
                             family="binomial", alpha=1), error=function(e) NULL)
    if (is.null(cv)) next
    cf <- as.matrix(coef(cv, s = "lambda.min"))[-1, ]
    g <- colnames(X)[cf != 0]
    sel_freq[g] <- sel_freq[g] + 1
  }
  lb <- data.frame(gene = names(sel_freq), boot_selection_freq = sel_freq / nboot)
  lb <- lb[order(-lb$boot_selection_freq), ]
  write.csv(lb, file.path(OUT, sprintf("lasso_bootstrap_stability_%s.csv", tag)), row.names = FALSE)
  cat(sprintf("[%s LASSO] bootstrap 200x: %d 基因入选率>=0.5, %d >=0.9, 最高 %s (%.2f)\n",
              tag, sum(lb$boot_selection_freq >= 0.5), sum(lb$boot_selection_freq >= 0.9),
              lb$gene[1], lb$boot_selection_freq[1]))

  # ---- 原始 LASSO / RF / SHAP 排名（基线规则复现） ----
  orig_dir <- if (tag == "ARDS") file.path(BASE, "05_ML_SHAP") else file.path(BASE, "05_ML_SHAP", "sepsis")
  rf_imp <- read.csv(file.path(orig_dir, if (tag=="ARDS") "RF_importance.csv" else "RF_importance_sepsis.csv"))
  shp_imp <- read.csv(file.path(orig_dir, if (tag=="ARDS") "SHAP_importance.csv" else "SHAP_importance_sepsis.csv"))
  rf_rank <- rf_imp$gene[order(-rf_imp$RF_gini)]
  shap_rank <- shp_imp$gene[order(-shp_imp$mean_abs_shap)]
  cvf <- cv.glmnet(X, as.numeric(y==lvl), family="binomial", alpha=1)
  cf <- as.matrix(coef(cvf, s="lambda.min"))[-1, ]
  lasso_genes <- colnames(X)[cf != 0]
  cat(sprintf("[%s] LASSO 原始选择 %d 基因\n", tag, length(lasso_genes)))

  key_set <- function(rf_k, shap_m) {
    ks <- lasso_genes[lasso_genes %in% union(rf_rank[1:rf_k], shap_rank[1:shap_m])]
    if (length(ks) < 5) ks <- shap_rank[1:min(15, length(shap_rank))]
    ks
  }
  base_set <- key_set(30, 20)

  # ---- B2(2) 阈值扰动 Jaccard ----
  jac <- function(a, b) length(intersect(a,b)) / length(union(a,b))
  grid <- expand.grid(rf_k = c(20,30,40), shap_m = c(10,20,30))
  grid$key_set <- sapply(seq_len(nrow(grid)), function(i) paste(key_set(grid$rf_k[i], grid$shap_m[i]), collapse=","))
  grid$jaccard_vs_base <- sapply(seq_len(nrow(grid)), function(i) jac(key_set(grid$rf_k[i], grid$shap_m[i]), base_set))
  write.csv(grid, file.path(OUT, sprintf("threshold_perturbation_jaccard_%s.csv", tag)), row.names = FALSE)
  cat(sprintf("[%s] 阈值扰动 Jaccard: 中位 %.2f, 范围 %.2f-%.2f\n",
              tag, median(grid$jaccard_vs_base), min(grid$jaccard_vs_base), max(grid$jaccard_vs_base)))
  # 六个 hub 基因在各扰动下是否保留
  hubs <- c("CDK5RAP2","CD177","GRB10","UPP1","EXOSC4","VNN1")
  keep <- sapply(seq_len(nrow(grid)), function(i) {
    ks <- key_set(grid$rf_k[i], grid$shap_m[i])
    mean(hubs %in% ks)
  })
  cat(sprintf("[%s] 六 hub 基因在各扰动下的保留率: %s\n", tag,
              paste(sprintf("%.0f%%", 100*keep), collapse=", ")))

  # ---- B2(3) RF 50 种子入选频率 ----
  seeds <- 101:150
  top30_freq <- setNames(numeric(length(ccgs)), ccgs)
  final_freq <- setNames(numeric(length(ccgs)), ccgs)
  for (s in seeds) {
    m <- randomForest(X, y, ntree = 1000, importance = TRUE, seed = s)
    gini <- importance(m)[, "MeanDecreaseGini"]
    rk <- names(sort(gini, decreasing = TRUE))
    top30_freq[rk[1:30]] <- top30_freq[rk[1:30]] + 1
    ks <- lasso_genes[lasso_genes %in% union(rk[1:30], shap_rank[1:20])]
    if (length(ks) < 5) ks <- shap_rank[1:15]
    final_freq[ks] <- final_freq[ks] + 1
  }
  rs <- data.frame(gene = names(top30_freq),
                   rf_top30_freq = top30_freq / length(seeds),
                   final_candidate_freq = final_freq / length(seeds))
  rs <- rs[order(-rs$final_candidate_freq, -rs$rf_top30_freq), ]
  write.csv(rs, file.path(OUT, sprintf("rf_seed_stability_%s.csv", tag)), row.names = FALSE)
  cat(sprintf("[%s RF-50seed] final 候选入选率>=0.8 的基因: %s\n", tag,
              paste(rs$gene[rs$final_candidate_freq >= 0.8], collapse=", ")))
  invisible(list(base_set = base_set, lasso_boot = lb, rf_seed = rs))
}

res_a <- run_cohort("ARDS",   "GSE76293", "ARDS")
res_s <- run_cohort("Sepsis", "GSE69528", "Sepsis")

# ---- B1(3) AUC = 1.000 精确单侧置信下限 ----
lb_one <- function(n1, n0, alpha) alpha^(1/(n1*n0))
bounds <- data.frame(
  cohort = c("GSE76293 (ARDS training)", "GSE243066 (ARDS in-cohort refit)",
             "GSE100159 (Sepsis in-cohort refit)"),
  n_cases = c(12, 34, 33), n_controls = c(12, 15, 12))
bounds$one_sided_97.5_LB <- mapply(lb_one, bounds$n_cases, bounds$n_controls, 0.025)
bounds$one_sided_95_LB   <- mapply(lb_one, bounds$n_cases, bounds$n_controls, 0.05)
bounds[] <- lapply(bounds, function(x) if (is.numeric(x)) round(x, 3) else x)
write.csv(bounds, file.path(OUT, "auc_exact_lower_bounds.csv"), row.names = FALSE)
print(bounds)

cat("\n=== Step 50 全部完成 ===\n")
