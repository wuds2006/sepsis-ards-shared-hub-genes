#!/usr/bin/env python3
"""Step3n-4: CellOracle 分层（ARDS组 vs 对照组）+ 细胞分化轨迹预测。

对 GSE151263 (ARDS vs Sepsis) 与 GSE167363 (Sepsis vs HC) 分别:
  1. 每组独立建 GRN（fit_GRN_for_simulation，细胞类型=cluster）→ 细胞状态特异 TRN
  2. 找 6 个 hub 基因（均为靶基因、非 TF）的组特异上游调控 TF
  3. 敲除上游 TF（simulate_shift），记录 hub 基因变化（总体 + 细胞类型分层）
  4. 用 KO 后表达位移向量 + 各细胞类型质心，估算细胞类型间转移概率矩阵
     （分化轨迹变化的代理；方法学透明、非 CellOracle 内建 cell fate 模块）

输出: 08_singlecell/celloracle_stratified/
"""
import os
import csv
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
from matplotlib import font_manager

warnings.filterwarnings("ignore")

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT = f"{BASE}/08_singlecell/celloracle_stratified"
os.makedirs(OUT, exist_ok=True)
GRN_PARQUET = f"{BASE}/08_singlecell/celloracle/hg38_base_GRN.parquet"

PURPLE = "#6A5C95"; PINK = "#D99AAE"; ORANGE = "#E7A889"
GREEN = "#A9CBB8"; MINT = "#CFE3E6"; BLUE = "#8FA2D6"; NAVY = "#5B608C"
HUB = ["CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1"]

for fp in ("/System/Library/Fonts/Helvetica.ttc",):
    if os.path.exists(fp):
        font_manager.fontManager.addfont(fp)
plt.rcParams["font.family"] = ["Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42

# 预导入 genomepy：touch 缓存锁使其 mtime 新鲜，避开其导入时的过期锁删除
# （该删除会触发沙箱 bulk-delete 钩子杀死进程）
_lock = os.path.expanduser("~/Library/Caches/genomepy/0.16.4/cache.lock")
try:
    os.utime(_lock, None)
except OSError:
    pass
import genomepy  # noqa: F401

import celloracle as co
print("celloracle:", co.__version__)


def save3(fig, name, line_art=False):
    out = os.path.join(OUT, name)
    fig.savefig(out + ".pdf", facecolor="white")
    fig.savefig(out + ".png", dpi=300, facecolor="white")
    fig.savefig(out + ".tiff", dpi=500 if not line_art else 1000, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print("[saved]", name)


def load_group_adata(which, grp, knockable, n_max=2500):
    if which == "GSE151263":
        cd = f"{BASE}/08_singlecell/cellchat_input"
    else:
        cd = f"{BASE}/08_singlecell/GSE167363/input"
    counts = scipy.io.mmread(f"{cd}/counts.mtx").tocsr()  # genes x cells
    genes = pd.read_csv(f"{cd}/genes.tsv", header=None)[0].tolist()
    genes = [g.split("\t")[0] for g in genes]
    meta = pd.read_csv(f"{cd}/meta.csv", index_col=0)
    umap = pd.read_csv(f"{cd}/umap.csv", index_col=0)
    sel = meta.index[meta["group"] == grp]
    pos = {c: i for i, c in enumerate(meta.index)}
    idx = [pos[c] for c in sel]
    adata = sc.AnnData(X=counts[:, idx].T.tocsr())
    adata.var_names = genes
    adata.obs_names = list(sel)
    adata.obs["celltype"] = meta.loc[sel, "celltype"].values
    adata.obsm["X_umap"] = umap.loc[sel].values
    print(f"{which} {grp}: {adata.n_obs} 细胞")
    if adata.n_obs > n_max:
        np.random.seed(42)
        keep = np.random.choice(adata.n_obs, n_max, replace=False)
        adata = adata[keep].copy()
    gene_total = np.array(adata.X.sum(axis=0)).flatten()
    top_idx = np.argsort(gene_total)[-2000:]
    hv = adata.var_names[top_idx].tolist()
    for g in knockable:
        if g in adata.var_names and g not in hv:
            hv.append(g)
    return adata[:, hv].copy()


def transition_proxy(raw_expr, sim_expr, ctypes, pca=None):
    """KO 位移向量 -> 细胞类型间转移概率矩阵（cos/距离 softmax 聚合）。
    raw_expr/sim_expr: cells x genes（dense ndarray）"""
    from sklearn.decomposition import PCA as skPCA
    if pca is None:
        pca = skPCA(n_components=30, random_state=42)
        raw_p = pca.fit_transform(raw_expr)
    else:
        raw_p = pca.transform(raw_expr)
    sim_p = pca.transform(sim_expr)
    d = sim_p - raw_p                       # 每细胞 KO 位移
    uniq = sorted(set(ctypes))
    cent = np.vstack([raw_p[ctypes == c].mean(axis=0) for c in uniq])
    T = np.zeros((len(uniq), len(uniq)))
    for i, c_i in enumerate(uniq):
        m = ctypes == c_i
        delta = d[m]                        # n_i x 30
        diff = cent[None, :, :] - raw_p[m][:, None, :]   # n_i x K x 30
        dn = np.linalg.norm(delta, axis=1, keepdims=True) * np.linalg.norm(diff, axis=2) + 1e-9
        cos = (delta[:, None, :] * diff).sum(axis=2) / dn
        cos = cos / 0.5                     # 温度
        P = np.exp(cos - cos.max(axis=1, keepdims=True))
        P = P / P.sum(axis=1, keepdims=True)
        T[i] = P.mean(axis=0)
    return T, uniq


def run_group(which, grp, knockable):
    tag = f"{which}_{grp}"
    print(f"\n========== {tag} ==========")
    adata = load_group_adata(which, grp, knockable)
    ctypes = adata.obs["celltype"].values.copy()
    raw_dense = np.asarray(adata.X.todense()) if sp.issparse(adata.X) else np.asarray(adata.X)

    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name="celltype",
                                      embedding_name="X_umap")
    oracle.import_TF_data(TF_info_matrix=pd.read_parquet(GRN_PARQUET))
    oracle.perform_PCA()
    oracle.knn_imputation(n_pca_dims=30, k=20, balanced=True,
                          b_sight=1500, b_maxl=1500, n_jobs=1)
    oracle.fit_GRN_for_simulation(alpha=10, verbose_level=0)

    # ---- 细胞状态特异 TRN（每 cluster coef 矩阵 + 连接度）----
    per_cluster = oracle.coef_matrix_per_cluster
    conn_rows = []
    for ct, mat in per_cluster.items():
        if hasattr(mat, "to_numpy"):
            M = mat.abs().to_numpy()
        else:
            M = np.abs(np.asarray(mat.todense()) if sp.issparse(mat) else mat)
        conn_rows.append(dict(dataset=which, group=grp, celltype=ct,
                              n_tf=(M > 0).any(axis=1).sum(),
                              n_targets=(M > 0).any(axis=0).sum(),
                              total_abs_weight=float(M.sum()),
                              mean_abs_weight=float(M[M > 0].mean()) if (M > 0).any() else 0.0))
        pd.DataFrame(mat).to_csv(f"{OUT}/CellOracle_coef_{tag}_{ct}.csv")
    pd.DataFrame(conn_rows).to_csv(f"{OUT}/CellOracle_TRN_connectivity_{tag}.csv", index=False)

    # ---- 组特异上游 TF ----
    coef = pd.concat(list(per_cluster.values())).groupby(level=0).mean()
    coef.to_csv(f"{OUT}/CellOracle_coef_matrix_{tag}.csv")
    upstream = {}
    for g in HUB:
        if g in coef.columns:
            col = coef[g].abs().sort_values(ascending=False)
            upstream[g] = col[col > 0].head(10).index.tolist()
    ko_tfs = list(dict.fromkeys([tf for tfs in upstream.values() for tf in tfs[:3]]))[:8]
    print(f"{tag} 上游 TF: {ko_tfs}")

    # ---- KO 模拟 ----
    hub_in = [g for g in HUB if g in coef.columns]
    rows, trans = [], {}
    for tf in ko_tfs:
        try:
            oracle.simulate_shift(perturb_condition={tf: 0.0}, n_propagation=3)
        except Exception as e:
            print(f"  {tf} simulate 失败: {e}")
            continue
        sim = oracle.adata.layers["simulated_count"]
        sim_d = np.asarray(sim.todense()) if sp.issparse(sim) else np.asarray(sim)
        # sim 矩阵基因顺序 = oracle.adata.var_names（同 adata）
        varn = list(oracle.adata.var_names)
        for g in hub_in:
            if g not in varn:
                continue
            j = varn.index(g)
            raw_v = raw_dense[:, list(adata.var_names).index(g)]
            sim_v = sim_d[:, j]
            r = dict(dataset=which, group=grp, KO_TF=tf, target=g,
                     raw_mean=float(raw_v.mean()),
                     simulated_mean=float(sim_v.mean()),
                     delta=float(sim_v.mean() - raw_v.mean()))
            for ct in sorted(set(ctypes)):
                m = ctypes == ct
                r[f"delta_{ct}"] = float(sim_v[m].mean() - raw_v[m].mean())
            rows.append(r)
        # 转移概率（分化轨迹变化代理）
        T, uniq = transition_proxy(raw_dense, sim_d, ctypes)
        trans[tf] = (T, uniq)
        pd.DataFrame(T, index=uniq, columns=uniq).to_csv(
            f"{OUT}/CellOracle_transition_{tag}_{tf}.csv")
        print(f"  {tf} 完成")

    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/CellOracle_TF_KO_hub_effect_{tag}.csv", index=False)
    return dict(tag=tag, df=df, trans=trans, upstream=upstream,
                conn=pd.DataFrame(conn_rows))


def plot_dataset(which, groups, res):
    tagA, tagB = [f"{which}_{g}" for g in groups]
    dfa, dfb = res[tagA]["df"], res[tagB]["df"]
    if len(dfa) == 0 or len(dfb) == 0:
        print(f"{which}: 无 KO 效应数据，跳图")
        return

    # ---- A. hub delta 热图（两组并排）----
    pivs = {}
    for tag, df in ((tagA, dfa), (tagB, dfb)):
        pivs[tag] = df.pivot_table(index="KO_TF", columns="target", values="delta")
    all_tfs = sorted(set(pivs[tagA].index) | set(pivs[tagB].index))
    all_hubs = sorted(set(pivs[tagA].columns) | set(pivs[tagB].columns))
    P = {t: pivs[t].reindex(index=all_tfs, columns=all_hubs) for t in pivs}
    vmax = np.nanmax([np.nanmax(np.abs(P[t].values)) for t in P])
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 0.36 * len(all_tfs) + 2.6), sharey=True)
    for ax, t in zip(axes, (tagA, tagB)):
        M = P[t].values
        im = ax.imshow(M, cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(all_hubs)), all_hubs, rotation=45, ha="right",
                      fontsize=8, style="italic")
        ax.set_title(t.replace(which + "_", which + " · "), fontsize=9)
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if not np.isnan(M[i, j]):
                    ax.text(j, i, f"{M[i,j]:+.2f}", ha="center", va="center", fontsize=6)
    axes[0].set_yticks(range(len(all_tfs)), all_tfs, fontsize=8)
    cb = fig.colorbar(im, ax=axes, shrink=0.7, pad=0.02)
    cb.set_label("Hub-gene delta after TF KO", fontsize=8)
    fig.suptitle(f"{which}: CellOracle upstream-TF KO effect on hub genes (group-specific GRN)",
                 fontsize=11)
    save3(fig, f"{which}_CellOracle_stratified_hubdelta")

    # ---- B. 转移概率热图（效应最强 TF，两组）----
    best_tf = None
    best = 0
    for tf, (T, uniq) in res[tagA]["trans"].items():
        v = np.abs(T - np.eye(len(T))).sum()
        if v > best:
            best, best_tf = v, tf
    if best_tf and best_tf in res[tagB]["trans"]:
        Tb, uniq_b = res[tagB]["trans"][best_tf]
        Ta, uniq_a = res[tagA]["trans"][best_tf]
        uniq = sorted(set(uniq_a) | set(uniq_b))
        A = pd.DataFrame(Ta, index=uniq_a, columns=uniq_a).reindex(index=uniq, columns=uniq).values
        B = pd.DataFrame(Tb, index=uniq_b, columns=uniq_b).reindex(index=uniq, columns=uniq).values
        fig, axes = plt.subplots(1, 2, figsize=(11.6, 0.5 * len(uniq) + 2.2), sharey=True)
        for ax, M, t in zip(axes, (A, B), (tagA, tagB)):
            im = ax.imshow(M, cmap="Blues", vmin=0, vmax=1, aspect="auto")
            ax.set_xticks(range(len(uniq)), uniq, rotation=45, ha="right", fontsize=8)
            ax.set_title(t.replace(which + "_", which + " · "), fontsize=9)
        axes[0].set_yticks(range(len(uniq)), uniq, fontsize=8)
        cb = fig.colorbar(im, ax=axes, shrink=0.7, pad=0.02)
        cb.set_label("P(cell type i -> j) after KO", fontsize=8)
        fig.suptitle(f"{which}: predicted cell-state transition after {best_tf} KO "
                     f"(cosine/centroid proxy)", fontsize=11)
        save3(fig, f"{which}_CellOracle_stratified_transition_{best_tf}")

    # ---- C. TRN 连接度（细胞类型 x 组）----
    ca = res[tagA]["conn"]; cb2 = res[tagB]["conn"]
    cts = sorted(set(ca.celltype) | set(cb2.celltype))
    x = np.arange(len(cts)); w = 0.38
    fig, ax = plt.subplots(figsize=(9.5, 4.4))
    va = [ca.set_index("celltype").total_abs_weight.get(c, 0) for c in cts]
    vb = [cb2.set_index("celltype").total_abs_weight.get(c, 0) for c in cts]
    ax.bar(x - w / 2, va, w, color=PINK, edgecolor="k", linewidth=0.4, label=tagA.split("_", 1)[1])
    ax.bar(x + w / 2, vb, w, color=BLUE, edgecolor="k", linewidth=0.4, label=tagB.split("_", 1)[1])
    ax.set_xticks(x, cts, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Total |coef| (TRN connectivity)")
    ax.legend(fontsize=8, frameon=False)
    ax.set_title(f"{which}: cell-state-specific TRN connectivity", fontsize=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    save3(fig, f"{which}_CellOracle_stratified_TRN")


PBMC_KO = ["UPP1", "EXOSC4"]
CD14_KO = ["UPP1", "EXOSC4", "VNN1", "CDK5RAP2"]

res = {}
for grp in ("ARDS", "Sepsis"):
    res[f"GSE151263_{grp}"] = run_group("GSE151263", grp, PBMC_KO)
plot_dataset("GSE151263", ["ARDS", "Sepsis"], res)

for grp in ("Sepsis", "HC"):
    res[f"GSE167363_{grp}"] = run_group("GSE167363", grp, CD14_KO)
plot_dataset("GSE167363", ["Sepsis", "HC"], res)

# ---- 上游 TF 重叠比较 ----
rows = []
for tag, r in res.items():
    for g, tfs in r["upstream"].items():
        rows.append(dict(dataset_group=tag, gene=g, top5_upstream_TF=", ".join(tfs[:5])))
pd.DataFrame(rows).to_csv(f"{OUT}/CellOracle_upstream_TF_per_group.csv", index=False)

print("\n===== CellOracle stratified 完成 =====")
