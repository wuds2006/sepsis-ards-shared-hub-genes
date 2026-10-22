#!/usr/bin/env python3
"""CellOracle 虚拟敲除（UPP1/EXOSC4）—— 基于官方 hg38 base GRN"""
import os, sys
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT = f"{BASE}/08_singlecell/celloracle"
os.makedirs(OUT, exist_ok=True)

import celloracle as co
print("celloracle:", co.__version__)

# ============ 1. 数据准备：从 raw counts 构建 h5ad ============
print("\n[1] 构建 raw counts AnnData ...")
counts = scipy.io.mmread(f"{BASE}/08_singlecell/cellchat_input/counts.mtx").tocsr()
genes = pd.read_csv(f"{BASE}/08_singlecell/cellchat_input/genes.tsv", header=None)[0].tolist()
barcodes = pd.read_csv(f"{BASE}/08_singlecell/cellchat_input/barcodes.tsv", header=None)[0].tolist()
meta = pd.read_csv(f"{BASE}/08_singlecell/cellchat_input/meta.csv", index_col=0)
print(f"counts: {counts.shape} (基因 x 细胞)")

# AnnData 是 细胞 x 基因
adata = sc.AnnData(X=counts.T.tocsr())
adata.var_names = genes
adata.obs_names = barcodes
adata.obs["celltype"] = meta.loc[barcodes, "celltype"].values
adata.obs["group"] = meta.loc[barcodes, "group"].values
print(f"AnnData: {adata.shape}, celltype: {adata.obs['celltype'].value_counts().to_dict()}")

# 加 X_umap（从已有 h5ad 复制）
h5ad = sc.read(f"{BASE}/08_singlecell/GSE151263_scanpy.h5ad")
# 对齐 barcode
common = list(set(barcodes) & set(h5ad.obs_names))
print(f"对齐 barcode: {len(common)}")
adata = adata[common].copy()
h5ad_sub = h5ad[common].copy()
adata.obsm["X_umap"] = h5ad_sub.obsm["X_umap"]

# 下采样到 3000 细胞（内存限制，避免 knn_imputation OOM）
if adata.n_obs > 3000:
    np.random.seed(42)
    keep = np.random.choice(adata.n_obs, 3000, replace=False)
    adata = adata[keep].copy()
    print(f"下采样到 {adata.n_obs} 细胞")

# 基因筛选：按表达量(总和)排序选 top 2000 高表达基因 + 强制纳入敲除目标
print("筛选高表达基因 ...")
gene_total = np.array(adata.X.sum(axis=0)).flatten()
top_idx = np.argsort(gene_total)[-2000:]
hv = adata.var_names[top_idx].tolist()
for g in ["UPP1", "EXOSC4"]:
    if g in adata.var_names and g not in hv:
        hv.append(g)
adata = adata[:, hv].copy()
print(f"筛选后基因数: {adata.n_vars}, 包含 UPP1: {'UPP1' in adata.var_names}, EXOSC4: {'EXOSC4' in adata.var_names}")

# ============ 2. 加载 base GRN ============
print("\n[2] 加载 hg38 base GRN ...")
base_grn = pd.read_parquet(f"{OUT}/hg38_base_GRN.parquet")
print(f"base GRN: {base_grn.shape}, 列: {list(base_grn.columns)}")

# ============ 3. 初始化 Oracle ============
print("\n[3] 初始化 Oracle ...")
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(
    adata=adata,
    cluster_column_name="celltype",
    embedding_name="X_umap")
oracle.import_TF_data(TF_info_matrix=base_grn)
print("Oracle 初始化完成, 基因数:", oracle.adata.n_vars)

# ============ 4. PCA + KNN imputation ============
print("\n[4] PCA + KNN imputation ...")
oracle.perform_PCA()
n_comps = 30
oracle.knn_imputation(n_pca_dims=n_comps, k=20, balanced=True,
                      b_sight=1500, b_maxl=1500, n_jobs=1)
print("KNN imputation 完成")

# ============ 5. GRN 拟合 ============
print("\n[5] fit_GRN_for_simulation ...")
oracle.fit_GRN_for_simulation(alpha=10, verbose_level=1)
# 用 coef_matrix_per_cluster（dict: cluster -> DataFrame[TF x 靶基因]），合并为整体（按 |coef| 平均）
per_cluster = oracle.coef_matrix_per_cluster
coef = pd.concat(list(per_cluster.values())).groupby(level=0).mean()
coef.to_csv(f"{OUT}/CellOracle_coef_matrix.csv")
print(f"coef_matrix: {coef.shape} (TF x 靶基因, {len(per_cluster)} cluster 平均)")

# ============ 6. 找 Hub 基因上游 TF + 敲除 ============
hub_targets = ["UPP1", "EXOSC4"]
# 找每个 Hub 基因的上游调控 TF（|coef| 最大的 TF）
upstream_TF = {}
for g in hub_targets:
    if g in coef.columns:
        col = coef[g].abs().sort_values(ascending=False)
        top_tfs = col[col > 0].head(10).index.tolist()
        upstream_TF[g] = top_tfs
        print(f"\n[6] {g} 上游调控 TF Top10: {', '.join(top_tfs[:6])}")
    else:
        print(f"\n[6] {g} 不在 coef_matrix 列中")

# 敲除上游 TF，观察 Hub 基因表达变化
print("\n[7] 敲除上游 TF，观察 Hub 基因变化 ...")
results = {}
all_ko_tfs = list(dict.fromkeys([tf for tfs in upstream_TF.values() for tf in tfs[:3]]))  # 取每基因 top3 TF，去重
print(f"敲除 TF 列表: {all_ko_tfs}")
for tf in all_ko_tfs:
    if tf not in oracle.adata.var_names and tf not in coef.index:
        print(f"  {tf} 不在 GRN, 跳过"); continue
    print(f"\n  敲除 TF: {tf} ...")
    try:
        oracle.simulate_shift(perturb_condition={tf: 0.0}, n_propagation=3)
        # 从 simulated_count 提取 Hub 基因的变化
        sim = oracle.adata.layers["simulated_count"]
        sim_df = pd.DataFrame(sim.toarray() if sp.issparse(sim) else sim,
                              index=oracle.adata.obs_names, columns=oracle.adata.var_names)
        for g in hub_targets:
            if g in sim_df.columns:
                delta = sim_df[g].mean() - np.array(oracle.adata.layers["raw_count"][:, list(oracle.adata.var_names).index(g)].todense()).flatten().mean() if sp.issparse(oracle.adata.layers["raw_count"]) else 0
                results[f"{tf}->{g}"] = {"KO_TF": tf, "target": g, "simulated_mean": sim_df[g].mean(), "delta": delta}
                print(f"    {tf} -> {g}: 模拟均值 {sim_df[g].mean():.4f}")
    except Exception as e:
        print(f"  {tf} 敲除失败: {e}")

pd.DataFrame(results.values()).to_csv(f"{OUT}/CellOracle_TF_KO_hub_effect.csv", index=False)
print("\n===== CellOracle TF 敲除对 Hub 基因影响完成 =====")

