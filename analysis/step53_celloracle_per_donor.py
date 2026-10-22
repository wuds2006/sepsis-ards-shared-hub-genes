#!/usr/bin/env python3
"""Step53 —— 供体级 CellOracle 敲除效应（审稿 W4/问题4）。

GSE167363 (CD14+ monocytes) 每个供体单独建 GRN：
  HC1, HC2, S2_T0, S2_T6（同一患者两时点）, S3
  对每个供体运行与 step3n 完全相同的 Oracle 流程，
  敲除 CEBPB（+该供体 UPP1 排名前 3 的上游 TF），记录 UPP1 等 hub 基因 delta。
输出: 08_singlecell/celloracle_per_donor/
  CellOracle_CEBPB_KO_per_donor.csv   每供体 CEBPB->hub delta（核心交付）
  CellOracle_TF_KO_hub_effect_<donor>.csv  全量
"""
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp

warnings.filterwarnings("ignore")

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT = f"{BASE}/08_singlecell/celloracle_per_donor"
os.makedirs(OUT, exist_ok=True)
GRN_PARQUET = f"{BASE}/08_singlecell/celloracle/hg38_base_GRN.parquet"

HUB = ["CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1"]
KNOCKABLE = ["CDK5RAP2", "UPP1", "EXOSC4", "VNN1"]   # CD177/GRB10 <1% 不适敲除（同主文口径）

_lock = os.path.expanduser("~/Library/Caches/genomepy/0.16.4/cache.lock")
try:
    os.utime(_lock, None)
except OSError:
    pass
import genomepy  # noqa: F401
import celloracle as co
print("celloracle:", co.__version__)

cd = f"{BASE}/08_singlecell/GSE167363/input"
counts = scipy.io.mmread(f"{cd}/counts.mtx").tocsr()   # genes x cells
genes = pd.read_csv(f"{cd}/genes.tsv", header=None)[0].tolist()
genes = [g.split("\t")[0] for g in genes]
meta = pd.read_csv(f"{cd}/meta.csv", index_col=0)
umap = pd.read_csv(f"{cd}/umap.csv", index_col=0)
pos = {c: i for i, c in enumerate(meta.index)}

def load_donor(sample):
    sel = meta.index[meta["sample"] == sample]
    idx = [pos[c] for c in sel]
    adata = sc.AnnData(X=counts[:, idx].T.tocsr())
    adata.var_names = genes
    adata.obs_names = list(sel)
    adata.obs["celltype"] = meta.loc[sel, "celltype"].values
    adata.obsm["X_umap"] = umap.loc[sel].values
    gene_total = np.array(adata.X.sum(axis=0)).flatten()
    top_idx = np.argsort(gene_total)[-2000:]
    hv = adata.var_names[top_idx].tolist()
    for g in KNOCKABLE:
        if g in adata.var_names and g not in hv:
            hv.append(g)
    return adata[:, hv].copy()

samples = ["GSM5511355_S3_T0"]   # 补跑：仅 S3（其余 4 个已完成）
all_cebpb = []

for smp in samples:
    don = smp.split("_", 1)[1] if "_" in smp else smp
    tag = smp
    print(f"\n========== 供体 {tag} ==========")
    adata = load_donor(smp)
    print(f"{tag}: {adata.n_obs} 细胞, {adata.n_vars} 基因")
    ctypes = adata.obs["celltype"].values.copy()
    raw_dense = np.asarray(adata.X.todense())

    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name="celltype",
                                       embedding_name="X_umap")
    oracle.import_TF_data(TF_info_matrix=pd.read_parquet(GRN_PARQUET))
    oracle.perform_PCA()
    oracle.knn_imputation(n_pca_dims=30, k=20, balanced=True,
                          b_sight=1500, b_maxl=1500, n_jobs=1)
    oracle.fit_GRN_for_simulation(alpha=10, verbose_level=0)

    per_cluster = oracle.coef_matrix_per_cluster
    coef = pd.concat(list(per_cluster.values())).groupby(level=0).mean()

    # 该供体 hub 基因的上游 TF（按平均 |coef| 排序）
    upstream = {}
    for g in HUB:
        if g in coef.columns:
            col = coef[g].abs().sort_values(ascending=False)
            upstream[g] = col[col > 0].head(10).index.tolist()
    ko_tfs = list(dict.fromkeys([tf for tfs in upstream.values() for tf in tfs[:3]]))[:8]
    if "CEBPB" not in ko_tfs:
        ko_tfs = ["CEBPB"] + ko_tfs      # 确保敲 CEBPB
    else:
        ko_tfs = ["CEBPB"] + [t for t in ko_tfs if t != "CEBPB"]
    print(f"{tag} 上游 TF (KO 列表): {ko_tfs}")

    hub_in = [g for g in HUB if g in coef.columns]
    rows = []
    for tf in ko_tfs:
        try:
            oracle.simulate_shift(perturb_condition={tf: 0.0}, n_propagation=3)
        except Exception as e:
            print(f"  {tf} simulate 失败: {e}")
            continue
        sim = oracle.adata.layers["simulated_count"]
        sim_d = np.asarray(sim.todense()) if sp.issparse(sim) else np.asarray(sim)
        varn = list(oracle.adata.var_names)
        for g in hub_in:
            if g not in varn:
                continue
            j = varn.index(g)
            raw_v = raw_dense[:, list(adata.var_names).index(g)]
            sim_v = sim_d[:, j]
            r = dict(donor=tag, KO_TF=tf, target=g,
                     raw_mean=float(raw_v.mean()),
                     simulated_mean=float(sim_v.mean()),
                     delta=float(sim_v.mean() - raw_v.mean()))
            for ct in sorted(set(ctypes)):
                m = ctypes == ct
                r[f"delta_{ct}"] = float(sim_v[m].mean() - raw_v[m].mean())
            rows.append(r)
        print(f"  {tf} 完成")

    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/CellOracle_TF_KO_hub_effect_{tag}.csv", index=False)
    cb = df[(df["KO_TF"] == "CEBPB") & (df["target"] == "UPP1")]
    if len(cb):
        all_cebpb.append(cb.iloc[0].to_dict())
        print(f"  >> {tag} CEBPB->UPP1 delta = {cb.iloc[0]['delta']:.3f}")

import glob
prev = [pd.read_csv(f) for f in glob.glob(f"{OUT}/CellOracle_TF_KO_hub_effect_GSM*.csv")]
allrows = pd.concat(prev, ignore_index=True)
cb = allrows[(allrows["KO_TF"] == "CEBPB") & (allrows["target"] == "UPP1")]
res = cb[["donor", "raw_mean", "simulated_mean", "delta"]].copy()
res.to_csv(f"{OUT}/CellOracle_CEBPB_KO_per_donor.csv", index=False)
print("\n===== 供体级 CEBPB->UPP1 效应 =====")
print(res[["donor", "raw_mean", "simulated_mean", "delta"]].to_string(index=False))
print("\n=== Step 53 完成 ===")
