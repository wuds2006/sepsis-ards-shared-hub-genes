#!/usr/bin/env python
"""Step3k-1: GSE151263 (PBMC) hub 图件 + 表达比例（基于 scanpy h5ad）。

背景: GSE151263_seurat_annotated.rds 文件损坏(gzip 截断)，但 GSE151263_scanpy.h5ad
完好且为各模块分析的实际数据源（13,996 细胞 × 21,210 基因, X=z-score, obsm 含 UMAP）。
X>0 占比与既有 KO 可行性表口径一致(UPP1 38.0%/EXOSC4 4.8%/VNN1 1.0%)。

输出（08_singlecell/）:
  GSE151263_hub_expression_pct.csv / GSE151263_hub_celltype_pct.csv
  GSE151263_hub_dotplot.{pdf,png,tiff}      点图: 基因×细胞型, 色=均值 韦=表达比例
  GSE151263_hub_featureplot.{pdf,png,tiff}  6 基因 UMAP 特征图
  GSE151263_hub_violin_group.{pdf,png,tiff} 6 基因 ARDS vs Sepsis 小提琴(+Mann-Whitney P)
  GSE151263_celltype_prop.{pdf,png,tiff}    细胞型组成条图
  GSE151263_UMAP.tiff                        (gs 500dpi 从矢量 PDF 渲染)
规范: PDF 矢量(fonttype 42) + PNG 300dpi + TIFF 500dpi LZW RGB，Nature 7 色。
"""
import os
import subprocess

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from scipy.stats import mannwhitneyu

ROOT = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
SC = os.path.join(ROOT, "08_singlecell")

PURPLE = "#6A5C95"; PINK = "#D99AAE"; ORANGE = "#E7A889"
GREEN = "#A9CBB8"; MINT = "#CFE3E6"; BLUE = "#8FA2D6"; NAVY = "#5B608C"
HUB = ["CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1"]
HUB_COLORS = dict(zip(HUB, [PURPLE, PINK, ORANGE, GREEN, BLUE, NAVY]))

for fp in ("/System/Library/Fonts/Helvetica.ttc",):
    if os.path.exists(fp):
        font_manager.fontManager.addfont(fp)
plt.rcParams["font.family"] = ["Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["axes.unicode_minus"] = True


def save3(fig, stem, w=None, h=None, dpi_tiff=500):
    fig.savefig(os.path.join(SC, stem + ".pdf"), facecolor="white")
    fig.savefig(os.path.join(SC, stem + ".png"), dpi=300, facecolor="white")
    fig.savefig(os.path.join(SC, stem + ".tiff"), dpi=dpi_tiff,
                facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print("saved:", stem)


def despine(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# ---------------------------------------------------------------- 读数据
print("[1] 读取 h5ad ...")
f = h5py.File(os.path.join(SC, "GSE151263_scanpy.h5ad"), "r")
obs = f["obs"]
var = f["var"]["_index"][:]
umap = np.array(f["obsm"]["X_umap"])
ct_codes = obs["celltype"]["codes"][:]
ct_cats = [c.decode() for c in obs["celltype"]["categories"][:]]
celltype = np.array([ct_cats[c] for c in ct_codes])
gr_codes = obs["group"]["codes"][:]
gr_cats = [c.decode() for c in obs["group"]["categories"][:]]
group = np.array([gr_cats[c] for c in gr_codes])

Xcols = {}
for g in HUB:
    i = int(np.where(var == g.encode())[0][0])
    Xcols[g] = np.array(f["X"][:, i])
f.close()

# ---------------------------------------------------------------- 2) 表达比例
print("[2] hub 表达比例 ...")
pct_df = [(g, float((Xcols[g] > 0).mean())) for g in HUB]
with open(os.path.join(SC, "GSE151263_hub_expression_pct.csv"), "w") as fh:
    fh.write('"gene","pct_expressing"\n')
    for g, p in pct_df:
        fh.write(f'"{g}",{round(p, 4)}\n')
print("   ", " | ".join(f"{g} {p:.1%}" for g, p in pct_df))

cts = sorted(set(celltype))
lines = ["celltype," + ",".join(HUB)]
for ct in cts:
    m = celltype == ct
    lines.append(ct + "," + ",".join(f"{float((Xcols[g][m] > 0).mean()):.4f}" for g in HUB))
with open(os.path.join(SC, "GSE151263_hub_celltype_pct.csv"), "w") as fh:
    fh.write("\n".join(lines) + "\n")

# ---------------------------------------------------------------- 3) DotPlot
print("[3] hub DotPlot ...")
mean_mat = np.zeros((len(cts), len(HUB)))
pct_mat = np.zeros((len(cts), len(HUB)))
for i, ct in enumerate(cts):
    m = celltype == ct
    for j, g in enumerate(HUB):
        mean_mat[i, j] = Xcols[g][m].mean()
        pct_mat[i, j] = (Xcols[g][m] > 0).mean()

fig, ax = plt.subplots(figsize=(7.2, 4.8))
vmax = np.abs(mean_mat).max()
norm = matplotlib.colors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    "rdw", ["#2166AC", "#FFFFFF", "#B2182B"])
for i in range(len(cts)):
    for j in range(len(HUB)):
        ax.scatter(j, i, s=30 + pct_mat[i, j] * 600, color=cmap(norm(mean_mat[i, j])),
                   edgecolors="#666666", linewidths=0.5, zorder=3)
ax.set_xticks(range(len(HUB)))
ax.set_xticklabels(HUB, fontsize=9, style="italic")
ax.set_yticks(range(len(cts)))
ax.set_yticklabels(cts, fontsize=9)
ax.set_xlim(-0.6, len(HUB) - 0.4)
ax.set_ylim(len(cts) - 0.4, -0.6)
despine(ax)
ax.spines["left"].set_position(("outward", 0))
ax.grid(color="#efefef", lw=0.5, zorder=1)
ax.set_axisbelow(True)
cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.75, pad=0.02)
cb.set_label("Mean expression (scaled)", fontsize=8)
cb.ax.tick_params(labelsize=7)
# size legend
for k, p in enumerate([0.1, 0.25, 0.5]):
    ax.scatter([], [], s=30 + p * 600, color="#DDDDDD", edgecolors="#666666",
               linewidths=0.5, label=f"{p:.0%}")
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.02), fontsize=7, frameon=False,
          title="Fraction expressing", title_fontsize=7)
ax.set_title("Hub gene expression by cell type (GSE151263 PBMC)",
             fontsize=9.5, fontweight="bold", loc="left", color="#2b2b2b", pad=10)
fig.subplots_adjust(left=0.16, right=0.84, top=0.90, bottom=0.12)
save3(fig, "GSE151263_hub_dotplot")

# ---------------------------------------------------------------- 4) FeaturePlot
print("[4] hub FeaturePlot ...")
fig, axes = plt.subplots(2, 3, figsize=(11, 7))
order = np.lexsort((umap[:, 1], umap[:, 0]))
for k, (ax, g) in enumerate(zip(axes.flat, HUB)):
    v = Xcols[g]
    vmax_g = np.percentile(v, 99)
    sc = ax.scatter(umap[:, 0], umap[:, 1], c=v, cmap="viridis", s=1.2, alpha=0.8,
                    vmin=0, vmax=max(vmax_g, 0.1), rasterized=True)
    ax.set_title(g, fontsize=10, style="italic")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color("#cccccc")
    cb = fig.colorbar(sc, ax=ax, shrink=0.7)
    cb.ax.tick_params(labelsize=6)
    cb.set_label("scaled expr.", fontsize=6)
fig.suptitle("Hub gene expression on UMAP (GSE151263 PBMC)", fontsize=10.5,
             fontweight="bold", x=0.02, ha="left")
fig.subplots_adjust(left=0.04, right=0.98, top=0.92, bottom=0.03, wspace=0.14, hspace=0.18)
save3(fig, "GSE151263_hub_featureplot")

# ---------------------------------------------------------------- 5) Violin by group
print("[5] hub Violin ...")
fig, axes = plt.subplots(2, 3, figsize=(11, 6.6))
for ax, g in zip(axes.flat, HUB):
    parts = ax.violinplot([Xcols[g][group == "ARDS"], Xcols[g][group == "Sepsis"]],
                          positions=[1, 2], showextrema=False, widths=0.8)
    for pc, c in zip(parts["bodies"], [PINK, BLUE]):
        pc.set_facecolor(c); pc.set_edgecolor("black"); pc.set_linewidth(0.6); pc.set_alpha(0.85)
    for pos, c in ((1, PINK), (2, BLUE)):
        d = Xcols[g][group == ("ARDS" if c == PINK else "Sepsis")]
        q1, med, q3 = np.percentile(d, [25, 50, 75])
        ax.vlines(pos, q1, q3, color="black", lw=1.1)
        ax.hlines(med, pos - 0.13, pos + 0.13, color="black", lw=1.6)
    u = mannwhitneyu(Xcols[g][group == "ARDS"], Xcols[g][group == "Sepsis"])
    p = u.pvalue
    ptxt = f"P = {p:.2e}" if p < 0.001 else f"P = {p:.3f}"
    ymax = max(np.percentile(Xcols[g], 99.5), 0.5)
    ax.set_ylim(min(np.percentile(Xcols[g], 0.5) - 0.2, -0.3), ymax * 1.18)
    ax.text(1.5, ymax * 1.10, ptxt, ha="center", fontsize=7.5, color="#2b2b2b")
    ax.set_xticks([1, 2]); ax.set_xticklabels(["ARDS", "Sepsis"], fontsize=8.5)
    ax.set_title(g, fontsize=9.5, style="italic")
    ax.set_ylabel("Scaled expression", fontsize=7.5)
    ax.tick_params(labelsize=7)
    despine(ax)
fig.suptitle("Hub gene expression: ARDS vs Sepsis (GSE151263 PBMC)",
             fontsize=10.5, fontweight="bold", x=0.02, ha="left")
fig.subplots_adjust(left=0.06, right=0.98, top=0.90, bottom=0.06, wspace=0.28, hspace=0.42)
save3(fig, "GSE151263_hub_violin_group")

# ---------------------------------------------------------------- 6) 细胞型比例图
print("[6] celltype proportion ...")
prop = {}
for gr in ("ARDS", "Sepsis"):
    m = group == gr
    vc = np.unique(celltype[m], return_counts=True)
    prop[gr] = {c: n / m.sum() for c, n in zip(*vc)}
fig, ax = plt.subplots(figsize=(7.6, 4.6))
x = np.arange(len(cts))
for k, (gr, c) in enumerate((("ARDS", PINK), ("Sepsis", BLUE))):
    vals = [prop[gr].get(ct, 0) * 100 for ct in cts]
    bars = ax.bar(x + (k - 0.5) * 0.4, vals, width=0.38, color=c, edgecolor="black",
                  linewidth=0.4, label=gr, zorder=3)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.4, f"{v:.1f}", ha="center",
                fontsize=6.6, color="#2b2b2b")
ax.set_xticks(x)
ax.set_xticklabels(cts, fontsize=9)
ax.set_ylabel("Cell fraction (%)", fontsize=8.5)
ax.legend(frameon=False, fontsize=8.5)
despine(ax)
ax.grid(axis="y", color="#efefef", lw=0.5, zorder=1)
ax.set_axisbelow(True)
ax.set_title("Cell type composition (GSE151263 PBMC)", fontsize=9.5,
             fontweight="bold", loc="left", color="#2b2b2b", pad=10)
fig.subplots_adjust(left=0.11, right=0.97, top=0.90, bottom=0.12)
save3(fig, "GSE151263_celltype_prop")

# ---------------------------------------------------------------- 7) UMAP tiff
print("[7] GSE151263_UMAP.tiff (gs 500dpi from vector PDF) ...")
subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=tiff24nc",
                "-r500", "-sOutputFile=" + os.path.join(SC, "GSE151263_UMAP.tiff"),
                os.path.join(SC, "GSE151263_UMAP.pdf")], check=True)
print("    done")
print("all done.")
