#!/usr/bin/env python3
"""
Step 3m —— GSE151263 PBMC 细胞注释证据图（与 GSE167363 marker_dotplot 对称）
输出:
  GSE151263_marker_dotplot.{pdf,png,tiff}   经典 marker × 7 类细胞 dotplot
  GSE151263_cluster_marker_scores.csv       各细胞类型 marker 打分汇总
口径: 与 hub dotplot 一致（h5ad z-score X, pct>0 = 表达率）
"""
import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import csv
import subprocess
import os

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
SC = os.path.join(BASE, "08_singlecell")

NATURE = ["#6A5C95", "#D99AAE", "#E7A889", "#A9CBB8", "#CFE3E6", "#8FA2D6", "#5B608C"]

# 经典 PBMC marker（与注释打分所用一致，去重保序）
MARKERS = [
    ("CD4_T", ["CD3D", "CD3E", "CD4", "IL7R"]),
    ("CD8_T", ["CD8A", "CD8B"]),
    ("NK", ["NKG7", "GNLY", "KLRD1"]),
    ("B", ["MS4A1", "CD79A", "CD19"]),
    ("Mono", ["CD14", "LYZ", "FCN1", "S100A8"]),
    ("DC", ["FCER1A", "CST3"]),
    ("Platelet", ["PPBP", "PF4"]),
]
CT_ORDER = ["CD4_T", "CD8_T", "NK", "B", "Mono", "DC", "Platelet"]

f = h5py.File(os.path.join(SC, "GSE151263_scanpy.h5ad"), "r")
obs = f["obs"]
var_index = np.array([v.decode() if isinstance(v, bytes) else v for v in f["var"]["_index"][:]])
g = obs["celltype"]
cats = [x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]]
celltype = np.array(cats)[g["codes"][:]]
X = f["X"]

# 汇总每个 marker 基因 × 每个细胞类型: mean z-score, pct>0
rows = []
plot_data = {}  # (gene, ct) -> (mean, pct)
gene_list = []
for ct_marker, genes in MARKERS:
    for gene in genes:
        if gene in gene_list:
            continue
        gene_list.append(gene)
        if gene not in var_index:
            for ct in CT_ORDER:
                plot_data[(gene, ct)] = (np.nan, 0.0)
            continue
        gi = int(np.where(var_index == gene)[0][0])
        col = X[:, gi]
        for ct in CT_ORDER:
            m = celltype == ct
            sub = col[m]
            mean = float(sub.mean())
            pct = float((sub > 0).mean() * 100)
            plot_data[(gene, ct)] = (mean, pct)
            rows.append({
                "gene": gene, "celltype": ct, "marker_of": ct_marker,
                "mean_zscore": round(mean, 4), "pct_expr": round(pct, 2),
            })

with open(os.path.join(SC, "GSE151263_cluster_marker_scores.csv"), "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=["gene", "celltype", "marker_of", "mean_zscore", "pct_expr"])
    w.writeheader()
    w.writerows(rows)
print("CSV 写出:", len(rows), "行")

# ---- dotplot ----
fig, ax = plt.subplots(figsize=(9.2, 5.6))
ng, nc = len(gene_list), len(CT_ORDER)
x = np.repeat(np.arange(ng), nc)
y = np.tile(np.arange(nc), ng)
means = np.array([plot_data[(ge, ct)][0] for ge in gene_list for ct in CT_ORDER])
pcts = np.array([plot_data[(ge, ct)][1] for ge in gene_list for ct in CT_ORDER])
means_n = np.nan_to_num(means, nan=-2.0)
sizes = np.clip(pcts, 0.5, 100) * 2.6
sc = ax.scatter(x, y, s=sizes, c=means_n, cmap="RdPu", vmin=-0.5, vmax=2.0,
                edgecolors="#555555", linewidths=0.25, zorder=3)
ax.set_xticks(np.arange(ng))
ax.set_xticklabels(gene_list, rotation=45, ha="right", fontsize=8.5)
ax.set_yticks(np.arange(nc))
ax.set_yticklabels(CT_ORDER, fontsize=9.5)
ax.set_xlim(-0.7, ng - 0.3)
ax.set_ylim(nc - 0.5, -0.6)
# marker 归属分隔线
bounds = []
acc = 0
for _, genes in MARKERS:
    acc += len([g2 for g2 in genes if g2 in gene_list[:acc + len(genes)]])  # 简化: 顺序唯一
    bounds.append(acc - 0.5)
for b in bounds[:-1]:
    ax.axvline(b, color="#bbbbbb", lw=0.6, zorder=1)
# 顶部归属标注
starts = [0] + [int(b + 0.5) for b in bounds[:-1]]
ends = [int(b + 0.5) for b in bounds]
for (ct_marker, _), s0, e0 in zip(MARKERS, starts, ends):
    ax.text((s0 + e0 - 1) / 2, -1.0, ct_marker, ha="center", va="bottom",
            fontsize=8, color="#444444", fontstyle="italic")
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
ax.grid(color="#f0f0f0", lw=0.5, zorder=0)
ax.set_axisbelow(True)
cb = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.10)
cb.set_label("Mean z-score", fontsize=8.5)
cb.ax.tick_params(labelsize=7.5)
handles = [Line2D([], [], marker="o", ls="", markersize=np.sqrt(p * 2.6),
                  markerfacecolor="#cccccc", markeredgecolor="#555555",
                  label=f"{p}%") for p in (20, 50, 80)]
leg = ax.legend(handles=handles, title="Pct. expressed", loc="center left",
                bbox_to_anchor=(1.14, 0.32), frameon=False, fontsize=7.5, title_fontsize=8)
ax.set_title("Canonical marker expression by cell type (GSE151263 PBMC)",
             fontsize=10.5, fontweight="bold", loc="left", pad=26)
fig.subplots_adjust(left=0.10, right=0.82, top=0.86, bottom=0.20)

stem = os.path.join(SC, "GSE151263_marker_dotplot")
fig.savefig(stem + ".pdf")
fig.savefig(stem + ".png", dpi=300)
plt.close(fig)
# tiff 500dpi（重渲染保证 dpi 元数据）
fig = plt.figure(figsize=(9.2, 5.6))
ax = fig.add_subplot(111)
sc = ax.scatter(x, y, s=sizes, c=means_n, cmap="RdPu", vmin=-0.5, vmax=2.0,
                edgecolors="#555555", linewidths=0.25, zorder=3)
ax.set_xticks(np.arange(ng)); ax.set_xticklabels(gene_list, rotation=45, ha="right", fontsize=8.5)
ax.set_yticks(np.arange(nc)); ax.set_yticklabels(CT_ORDER, fontsize=9.5)
ax.set_xlim(-0.7, ng - 0.3); ax.set_ylim(nc - 0.5, -0.6)
for b in bounds[:-1]:
    ax.axvline(b, color="#bbbbbb", lw=0.6, zorder=1)
for (ct_marker, _), s0, e0 in zip(MARKERS, starts, ends):
    ax.text((s0 + e0 - 1) / 2, -1.0, ct_marker, ha="center", va="bottom",
            fontsize=8, color="#444444", fontstyle="italic")
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
ax.grid(color="#f0f0f0", lw=0.5, zorder=0); ax.set_axisbelow(True)
cb = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.10)
cb.set_label("Mean z-score", fontsize=8.5); cb.ax.tick_params(labelsize=7.5)
ax.legend(handles=handles, title="Pct. expressed", loc="center left",
          bbox_to_anchor=(1.14, 0.32), frameon=False, fontsize=7.5, title_fontsize=8)
ax.set_title("Canonical marker expression by cell type (GSE151263 PBMC)",
             fontsize=10.5, fontweight="bold", loc="left", pad=26)
fig.subplots_adjust(left=0.10, right=0.82, top=0.86, bottom=0.20)
fig.savefig(stem + ".tiff", dpi=500, format="tiff", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)
print("三格式完成:", stem)
