#!/usr/bin/env python3
"""CellOracle 虚拟敲除 —— GSE167363 (Sepsis CD14+ Mono)
方法: 基于官方 hg38 base GRN, fit_GRN_for_simulation 拟合 GRN,
敲除 hub 基因(UPP1/CDK5RAP2/EXOSC4/VNN1)的上游调控 TF, 观察 hub 基因表达变化。
(UPP1 等 hub 基因为靶基因非 TF, CellOracle 只能敲 TF -> 敲上游 TF 间接验证)
"""
import os, sys
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT = f"{BASE}/08_singlecell/GSE167363/celloracle"
INP = f"{BASE}/08_singlecell/GSE167363/input"
GRN_PARQUET = f"{BASE}/08_singlecell/celloracle/hg38_base_GRN.parquet"
os.makedirs(OUT, exist_ok=True)

import celloracle as co
print("celloracle:", co.__version__)

HUB = ["CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1"]
# 表达量足够、预期有 GRN 关联的敲除目标
KNOCKABLE = ["UPP1", "CDK5RAP2", "EXOSC4", "VNN1"]

# ============ 1. 数据准备 ============
print("\n[1] 构建 raw counts AnnData ...")
counts = scipy.io.mmread(f"{INP}/counts.mtx").tocsr()
genes = pd.read_csv(f"{INP}/genes.tsv", header=None)[0].tolist()
meta = pd.read_csv(f"{INP}/meta.csv", index_col=0)
umap = pd.read_csv(f"{INP}/umap.csv", index_col=0)
print(f"counts: {counts.shape} (基因 x 细胞)")

adata = sc.AnnData(X=counts.T.tocsr())
adata.var_names = genes
adata.obs_names = meta.index.tolist()
adata.obs["celltype"] = meta["celltype"].values
adata.obs["group"] = meta["group"].values
adata.obsm["X_umap"] = umap.loc[adata.obs_names].values
print(f"AnnData: {adata.shape}, celltype: {adata.obs['celltype'].value_counts().to_dict()}")

# 下采样 3000 细胞
if adata.n_obs > 3000:
    np.random.seed(42)
    keep = np.random.choice(adata.n_obs, 3000, replace=False)
    adata = adata[keep].copy()
    print(f"下采样到 {adata.n_obs} 细胞")

# 基因筛选: top 2000 高表达 + 强制纳入敲除目标
gene_total = np.array(adata.X.sum(axis=0)).flatten()
top_idx = np.argsort(gene_total)[-2000:]
hv = adata.var_names[top_idx].tolist()
for g in KNOCKABLE:
    if g in adata.var_names and g not in hv:
        hv.append(g)
adata = adata[:, hv].copy()
print(f"筛选后基因数: {adata.n_vars}; hub 在列: {[g for g in KNOCKABLE if g in adata.var_names]}")

# ============ 2. base GRN ============
print("\n[2] 加载 hg38 base GRN ...")
base_grn = pd.read_parquet(GRN_PARQUET)
print(f"base GRN: {base_grn.shape}")

# ============ 3-5. Oracle ============
print("\n[3] 初始化 Oracle ...")
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name="celltype",
                                   embedding_name="X_umap")
oracle.import_TF_data(TF_info_matrix=base_grn)

print("\n[4] PCA + KNN imputation ...")
oracle.perform_PCA()
oracle.knn_imputation(n_pca_dims=30, k=20, balanced=True,
                      b_sight=1500, b_maxl=1500, n_jobs=1)

print("\n[5] fit_GRN_for_simulation ...")
oracle.fit_GRN_for_simulation(alpha=10, verbose_level=1)
per_cluster = oracle.coef_matrix_per_cluster
coef = pd.concat(list(per_cluster.values())).groupby(level=0).mean()
coef.to_csv(f"{OUT}/CellOracle_coef_matrix.csv")
print(f"coef_matrix: {coef.shape} (TF x 靶基因, {len(per_cluster)} cluster 平均)")

# ============ 6. 找 hub 基因上游 TF ============
hub_targets = [g for g in KNOCKABLE if g in coef.columns]
print(f"\n[6] GRN 中含 hub 靶基因: {hub_targets}")
upstream_TF = {}
for g in hub_targets:
    col = coef[g].abs().sort_values(ascending=False)
    top_tfs = col[col > 0].head(10).index.tolist()
    upstream_TF[g] = top_tfs
    print(f"  {g} 上游 TF Top6: {', '.join(top_tfs[:6])}")

# ============ 7. 敲除上游 TF ============
print("\n[7] 敲除上游 TF，观察 hub 基因变化 ...")
results = []
all_ko_tfs = list(dict.fromkeys([tf for tfs in upstream_TF.values() for tf in tfs[:3]]))
print(f"敲除 TF 列表: {all_ko_tfs}")
for tf in all_ko_tfs:
    print(f"\n  敲除 TF: {tf} ...")
    try:
        oracle.simulate_shift(perturb_condition={tf: 0.0}, n_propagation=3)
        sim = oracle.adata.layers["simulated_count"]
        sim_df = pd.DataFrame(sim.toarray() if sp.issparse(sim) else sim,
                              index=oracle.adata.obs_names, columns=oracle.adata.var_names)
        raw = oracle.adata.layers["raw_count"]
        raw_df = pd.DataFrame(raw.toarray() if sp.issparse(raw) else raw,
                              index=oracle.adata.obs_names, columns=oracle.adata.var_names)
        for g in hub_targets:
            if g in sim_df.columns:
                delta = sim_df[g].mean() - raw_df[g].mean()
                results.append({"KO_TF": tf, "target": g,
                                "raw_mean": raw_df[g].mean(),
                                "simulated_mean": sim_df[g].mean(),
                                "delta": delta})
                print(f"    {tf} -> {g}: delta={delta:+.4f}")
    except Exception as e:
        print(f"  {tf} 敲除失败: {e}")

df_res = pd.DataFrame(results)
df_res.to_csv(f"{OUT}/CellOracle_TF_KO_hub_effect.csv", index=False)

# ============ 8. 汇总图 ============
if len(df_res) > 0:
    import matplotlib.pyplot as plt
    # 每 hub 基因: TF 敲除效应排序条形图
    hub_pal = ["#6A5C95", "#D99AAE", "#E7A889", "#A9CBB8", "#CFE3E6", "#8FA2D6", "#5B608C"]
    for g in hub_targets:
        sub = df_res[df_res.target == g].sort_values("delta")
        fig, ax = plt.subplots(figsize=(6, 4.5))
        colors = ["#8FA2D6" if v < 0 else "#D99AAE" for v in sub.delta]
        ax.barh(sub.KO_TF, sub.delta, color=colors)
        ax.axvline(0, color="#5B608C", lw=0.8)
        ax.set_xlabel("Expression delta after TF KO")
        ax.set_title(f"CellOracle: upstream TF KO effect on {g}")
        fig.tight_layout()
        for ext, dpi in (("pdf", None), ("png", 300), ("tiff", 500)):
            if dpi:
                fig.savefig(f"{OUT}/CellOracle_KO_effect_{g}.{ext}", dpi=dpi)
            else:
                fig.savefig(f"{OUT}/CellOracle_KO_effect_{g}.{ext}")
        plt.close(fig)

    # 热图: TF x hub delta
    pv = df_res.pivot_table(index="KO_TF", columns="target", values="delta")
    pv.to_csv(f"{OUT}/CellOracle_TF_KO_matrix.csv")
    fig, ax = plt.subplots(figsize=(6, max(4, 0.32 * len(pv))))
    vmax = np.nanmax(np.abs(pv.values))
    im = ax.imshow(pv.values, cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(pv.shape[1]), pv.columns, rotation=45, ha="right")
    ax.set_yticks(range(pv.shape[0]), pv.index)
    for i in range(pv.shape[0]):
        for j in range(pv.shape[1]):
            ax.text(j, i, f"{pv.values[i,j]:+.3f}", ha="center", va="center", fontsize=6)
    fig.colorbar(im, label="Expression delta")
    ax.set_title("CellOracle TF KO -> hub gene effect (GSE167363)")
    fig.tight_layout()
    for ext, dpi in (("pdf", None), ("png", 300), ("tiff", 500)):
        if dpi:
            fig.savefig(f"{OUT}/CellOracle_KO_heatmap.{ext}", dpi=dpi)
        else:
            fig.savefig(f"{OUT}/CellOracle_KO_heatmap.{ext}")
    plt.close(fig)

print("\n===== CellOracle (GSE167363) 完成 =====")
