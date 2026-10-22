#!/usr/bin/env python
"""Step3n-2: 代谢通路分析 — GSE151263 (PBMC) & GSE167363 (CD14+ Mono)。

使用 MSigDB Hallmark 代谢基因集（h.all.v2024.1.Hs.symbols.gmt）：
  GLYCOLYSIS / OXIDATIVE_PHOSPHORYLATION / FATTY_ACID_METABOLISM /
  CHOLESTEROL_HOMEOSTASIS / HEME_METABOLISM / BILE_ACID_METABOLISM / XENOBIOTIC_METABOLISM

分析内容：
  1. 每细胞代谢通路活性评分（scanpy score_genes，log-normalized counts）
  2. 细胞类型 × 通路 热图（按分组分层）
  3. 分组间通路差异（disease vs control，Wilcoxon，按数据集）
  4. 6 个 hub 基因表达与代谢通路的单细胞 Spearman 相关
  5. PBMC UMAP 上糖酵解/OXPHOS 活性叠加

输出（08_singlecell/metabolism/）:
  GSE151263_metab_pathway_heatmap.{pdf,png,tiff}
  GSE151263_metab_group_difference.{pdf,png,tiff}
  GSE151263_metab_hub_correlation.{pdf,png,tiff}
  GSE151263_metab_umap_overlay.{pdf,png,tiff}
  GSE167363 同名文件（无 UMAP 面板时仅前 3 张）
  metabolism_celltype_scores_{ds}.csv / metabolism_group_stats_{ds}.csv /
  metabolism_hub_correlation_{ds}.csv
"""
import os
import csv
from collections import OrderedDict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import scipy.sparse as sp
from matplotlib import font_manager
from scipy.stats import spearmanr, mannwhitneyu

import scanpy as sc
import anndata as ad

ROOT = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
SC = os.path.join(ROOT, "08_singlecell")
OUT = os.path.join(SC, "metabolism")
os.makedirs(OUT, exist_ok=True)

PURPLE = "#6A5C95"; PINK = "#D99AAE"; ORANGE = "#E7A889"
GREEN = "#A9CBB8"; MINT = "#CFE3E6"; BLUE = "#8FA2D6"; NAVY = "#5B608C"
PAL7 = [PURPLE, PINK, ORANGE, GREEN, MINT, BLUE, NAVY]
HUB = ["CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1"]

PATHWAYS = OrderedDict([
    ("Glycolysis", "HALLMARK_GLYCOLYSIS"),
    ("Oxidative phosphorylation", "HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
    ("Fatty acid metabolism", "HALLMARK_FATTY_ACID_METABOLISM"),
    ("Cholesterol homeostasis", "HALLMARK_CHOLESTEROL_HOMEOSTASIS"),
    ("Heme metabolism", "HALLMARK_HEME_METABOLISM"),
    ("Bile acid metabolism", "HALLMARK_BILE_ACID_METABOLISM"),
    ("Xenobiotic metabolism", "HALLMARK_XENOBIOTIC_METABOLISM"),
])
PLABELS = list(PATHWAYS.keys())

for fp in ("/System/Library/Fonts/Helvetica.ttc",):
    if os.path.exists(fp):
        font_manager.fontManager.addfont(fp)
plt.rcParams["font.family"] = ["Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42


def save3(fig, name, line_art=False):
    out = os.path.join(OUT, name)
    fig.savefig(out + ".pdf", facecolor="white")
    fig.savefig(out + ".png", dpi=300, facecolor="white")
    fig.savefig(out + ".tiff", dpi=500 if not line_art else 1000, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print("[saved]", name)


def load_gmt():
    gsets = {}
    for line in open(os.path.join(ROOT, "00_rawdata/hallmark/h.all.v2024.1.Hs.symbols.gmt")):
        f = line.rstrip("\n").split("\t")
        gsets[f[0]] = [g for g in f[2:] if g]
    return gsets


GMT = load_gmt()


def load_dataset(which):
    if which == "GSE151263":
        cd = os.path.join(SC, "cellchat_input")
    else:
        cd = os.path.join(SC, "GSE167363/input")
    cnt = sio.mmread(os.path.join(cd, "counts.mtx")).T.tocsr()  # cells x genes
    genes = [l.strip().split("\t")[0] for l in open(os.path.join(cd, "genes.tsv"))]
    rows = list(csv.DictReader(open(os.path.join(cd, "meta.csv"))))
    cell_key = "cell" if rows and "cell" in rows[0] else None
    import pandas as pd
    obs = pd.DataFrame(rows)
    if cell_key is None:
        # 第一列为无名细胞 barcode（GSE167363 meta.csv 首列为空表头）
        cell_key = obs.columns[0]
    obs = obs.set_index(cell_key)
    A = ad.AnnData(X=sp.csr_matrix(cnt, dtype=np.float32),
                   obs=obs, var=pd.DataFrame(index=genes))
    return A


def analyze(which, group_case, group_ctrl, umap_xy=None):
    print("\n=====", which, "=====")
    A = load_dataset(which)
    # log-normalize（CPM -> log1p）
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)

    # 通路评分（存于 obs）
    for label, gs in PATHWAYS.items():
        genes = [g for g in GMT[gs] if g in A.var_names]
        sc.tl.score_genes(A, gene_list=genes, score_name=f"metab_{gs}", random_state=42)
    S = np.column_stack([A.obs[f"metab_{PATHWAYS[p]}"].to_numpy() for p in PLABELS])

    ctypes = np.array(A.obs["celltype"])
    groups = np.array(A.obs["group"])
    uniq_ct = sorted(set(ctypes))
    case_col, ctrl_col = PINK, BLUE

    # ---------- 1. 细胞类型 x 通路热图（分组分层） ----------
    mat = {g: np.zeros((len(uniq_ct), len(PLABELS))) for g in (group_case, group_ctrl)}
    for gi, g in enumerate((group_case, group_ctrl)):
        m = groups == g
        for i, ct in enumerate(uniq_ct):
            mm = m & (ctypes == ct)
            if mm.sum() == 0:
                mat[g][i] = np.nan
                continue
            mat[g][i] = S[mm].mean(axis=0)

    with open(os.path.join(OUT, f"metabolism_celltype_scores_{which}.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["group", "celltype"] + PLABELS)
        for g in (group_case, group_ctrl):
            for i, ct in enumerate(uniq_ct):
                w.writerow([g, ct] + [f"{v:.4f}" for v in mat[g][i]])

    allvals = np.concatenate([mat[g].ravel() for g in (group_case, group_ctrl)])
    allvals = allvals[~np.isnan(allvals)]
    vmin, vmax = np.percentile(allvals, 2), np.percentile(allvals, 98)
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 0.5 * len(uniq_ct) + 2.4), sharey=True)
    for ax, g in zip(axes, (group_case, group_ctrl)):
        M = np.ma.masked_invalid(mat[g])
        im = ax.imshow(M, cmap="RdBu_r", vmin=vmin, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(PLABELS)), ["Glyco", "OXPHOS", "FA", "Chol", "Heme", "Bile", "Xeno"],
                      rotation=40, ha="right", fontsize=8)
        ax.set_title(g, fontsize=10)
    axes[0].set_yticks(range(len(uniq_ct)), uniq_ct, fontsize=8)
    cb = fig.colorbar(im, ax=axes, shrink=0.75, pad=0.02)
    cb.set_label("Mean module score", fontsize=8)
    fig.suptitle(f"{which}: metabolic pathway activity by cell type", fontsize=11)
    save3(fig, f"{which}_metab_pathway_heatmap")

    # ---------- 2. 分组间通路差异 ----------
    stats_rows = []
    for j, p in enumerate(PLABELS):
        x = S[groups == group_case, j]
        y = S[groups == group_ctrl, j]
        try:
            u, pv = mannwhitneyu(x, y)
        except ValueError:
            pv = 1.0
        stats_rows.append(dict(pathway=p, mean_case=float(x.mean()), mean_ctrl=float(y.mean()),
                               delta=float(x.mean() - y.mean()), p=float(pv),
                               p_log10=-np.log10(max(pv, 1e-300))))
    with open(os.path.join(OUT, f"metabolism_group_stats_{which}.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(stats_rows[0].keys())); w.writeheader()
        for r in stats_rows: w.writerow(r)

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.3), gridspec_kw={"width_ratios": [1, 1]})
    ax = axes[0]
    x = np.arange(len(PLABELS)); wdt = 0.38
    ax.bar(x - wdt / 2, [r["mean_case"] for r in stats_rows], wdt, color=case_col,
           edgecolor="k", linewidth=0.4, label=group_case)
    ax.bar(x + wdt / 2, [r["mean_ctrl"] for r in stats_rows], wdt, color=ctrl_col,
           edgecolor="k", linewidth=0.4, label=group_ctrl)
    ax.set_xticks(x, ["Glyco", "OXPHOS", "FA", "Chol", "Heme", "Bile", "Xeno"], fontsize=8)
    ax.set_ylabel("Mean pathway score"); ax.legend(fontsize=8, frameon=False)
    ax.set_title("Mean activity", fontsize=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    ax = axes[1]
    cols = [PURPLE if r["delta"] > 0 else ORANGE for r in stats_rows]
    ax.barh(np.arange(len(PLABELS)), [r["delta"] for r in stats_rows], color=cols,
            edgecolor="k", linewidth=0.4)
    for i, r in enumerate(stats_rows):
        star = "***" if r["p"] < 0.001 else ("**" if r["p"] < 0.01 else
                                             ("*" if r["p"] < 0.05 else "ns"))
        ax.text(r["delta"] + (0.004 if r["delta"] >= 0 else -0.004), i, star, va="center",
                ha="left" if r["delta"] >= 0 else "right", fontsize=7.5)
    ax.set_yticks(np.arange(len(PLABELS)), PLABELS, fontsize=8)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel(f"Delta ({group_case} - {group_ctrl})")
    ax.set_title("Disease vs control (Mann-Whitney)", fontsize=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.suptitle(f"{which}: metabolic reprogramming ({group_case} vs {group_ctrl})", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save3(fig, f"{which}_metab_group_difference")

    # ---------- 3. hub 基因与代谢通路相关 ----------
    cor = np.full((len(HUB), len(PLABELS)), np.nan)
    pv = np.full_like(cor, np.nan)
    Xexpr = A.to_df()  # 基因表达矩阵（log-normalized）
    for i, g in enumerate(HUB):
        if g not in Xexpr.columns:
            continue
        gx = Xexpr[g].to_numpy()
        for j, p in enumerate(PLABELS):
            sx = S[:, j]
            keep = (gx > 0)
            if keep.sum() < 30:
                keep = np.ones_like(gx, dtype=bool)
            r_, p_ = spearmanr(gx[keep], sx[keep])
            cor[i, j] = r_; pv[i, j] = p_
    with open(os.path.join(OUT, f"metabolism_hub_correlation_{which}.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["gene"] + PLABELS)
        for i, g in enumerate(HUB):
            w.writerow([g] + [f"{cor[i,j]:.3f} (p={pv[i,j]:.1e})" if not np.isnan(cor[i, j]) else "NA"
                              for j in range(len(PLABELS))])

    fig, ax = plt.subplots(figsize=(8.8, 4.6))
    im = ax.imshow(cor, cmap="RdYlBu_r", vmin=-0.35, vmax=0.35)
    ax.set_xticks(range(len(PLABELS)), PLABELS, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(HUB)), HUB, fontsize=8, style="italic")
    for i in range(len(HUB)):
        for j in range(len(PLABELS)):
            if np.isnan(cor[i, j]):
                continue
            star = "***" if pv[i, j] < 0.001 else ("**" if pv[i, j] < 0.01 else
                                                  ("*" if pv[i, j] < 0.05 else ""))
            ax.text(j, i, star, ha="center", va="center", fontsize=7)
    cb = fig.colorbar(im, ax=ax, shrink=0.85)
    cb.set_label("Spearman rho", fontsize=8)
    ax.set_title(f"{which}: hub-gene expression vs metabolic pathway activity\n(* p<0.05, ** p<0.01, *** p<0.001)",
                 fontsize=10)
    save3(fig, f"{which}_metab_hub_correlation")

    # ---------- 4. UMAP 叠加（仅 PBMC 有现成 UMAP） ----------
    if which == "GSE151263":
        hb = ad.read_h5ad(os.path.join(SC, "GSE151263_scanpy.h5ad"))
        hb_obs_index = list(hb.obs_names)
        common = [c for c in hb_obs_index if c in A.obs_names]
        pos = {c: i for i, c in enumerate(A.obs_names)}
        idx = np.array([pos[c] for c in common])
        gly = S[idx, 0]
        oxp = S[idx, 1]
        xy = hb[common].obsm["X_umap"]
        fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6))
        for ax, val, ttl, cmap in ((axes[0], gly, "Glycolysis", "magma"),
                                   (axes[1], oxp, "Oxidative phosphorylation", "viridis")):
            o = np.argsort(val)
            sc_ = ax.scatter(xy[o, 0], xy[o, 1], c=val[o], cmap=cmap, s=1.5, lw=0)
            cb = fig.colorbar(sc_, ax=ax, shrink=0.8)
            cb.set_label("Module score", fontsize=8)
            ax.set_title(ttl, fontsize=10); ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle(f"{which}: metabolic pathway activity on UMAP", fontsize=11)
        save3(fig, f"{which}_metab_umap_overlay")


analyze("GSE151263", "ARDS", "Sepsis")
analyze("GSE167363", "Sepsis", "HC")
print("\nall done.")
