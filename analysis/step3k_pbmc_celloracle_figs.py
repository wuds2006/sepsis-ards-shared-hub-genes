#!/usr/bin/env python
"""Step3k-2: GSE151263 CellOracle 图件补齐 + UMAP TIFF。

1) GSE151263_UMAP.tiff — 从现有矢量 PDF 用 gs 500dpi 渲染（彩色图规范）
2) CellOracle PBMC 图件（此前只有 CSV，无图）:
   - CellOracle_KO_effect_UPP1 / _EXOSC4 — TF 敲除 delta 条形图
   - CellOracle_KO_heatmap — 上游 TF coef 热图（top TF x hub 基因）
规范: PDF 矢量(fonttype 42) + PNG 300dpi + TIFF 500dpi LZW RGB，Nature 7 色。
"""
import csv
import os
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

ROOT = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
SC = os.path.join(ROOT, "08_singlecell")
CO = os.path.join(SC, "celloracle")

PURPLE = "#6A5C95"; PINK = "#D99AAE"; ORANGE = "#E7A889"
GREEN = "#A9CBB8"; MINT = "#CFE3E6"; BLUE = "#8FA2D6"; NAVY = "#5B608C"
HUB_COLORS = {"CDK5RAP2": PURPLE, "CD177": PINK, "GRB10": ORANGE,
              "UPP1": GREEN, "EXOSC4": BLUE, "VNN1": NAVY}

for fp in ("/System/Library/Fonts/Helvetica.ttc",):
    if os.path.exists(fp):
        font_manager.fontManager.addfont(fp)
plt.rcParams["font.family"] = ["Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["axes.unicode_minus"] = True


def save3(fig, stem, outdir, w, h, dpi_tiff=500):
    fig.savefig(os.path.join(outdir, stem + ".pdf"), facecolor="white")
    fig.savefig(os.path.join(outdir, stem + ".png"), dpi=300, facecolor="white")
    fig.savefig(os.path.join(outdir, stem + ".tiff"), dpi=dpi_tiff,
                facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print("saved:", stem)


# ---------------------------------------------------------------- 1) UMAP tiff
print("[1] GSE151263_UMAP.tiff (from vector PDF, 500dpi color) ...")
subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=tiff24nc",
                "-r500", "-sOutputFile=" + os.path.join(SC, "GSE151263_UMAP.tiff"),
                os.path.join(SC, "GSE151263_UMAP.pdf")], check=True)
print("    done")

# ---------------------------------------------------------------- 2) CellOracle 图件
print("[2] CellOracle PBMC figures ...")

# 2a. TF KO effect (delta)
rows = list(csv.DictReader(open(os.path.join(CO, "CellOracle_TF_KO_hub_effect.csv"))))
eff = {}
for r in rows:
    eff.setdefault(r["target"], []).append((r["KO_TF"], float(r["delta"])))

for target, items in eff.items():
    items = sorted(items, key=lambda x: x[1])          # 负向效应最大在前
    tfs = [t for t, _ in items]
    deltas = [d for _, d in items]
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    cols = [HUB_COLORS.get(target, PURPLE)] * len(tfs)
    bars = ax.bar(range(len(tfs)), deltas, color=cols, width=0.62,
                  edgecolor="white", linewidth=0.8, zorder=3)
    ax.axhline(0, color="#555555", lw=0.9, zorder=2)
    for i, (b, d) in enumerate(zip(bars, deltas)):
        va = "bottom" if d < 0 else "top"
        off = 0.012 * (max(deltas) - min(deltas))
        ax.text(b.get_x() + b.get_width() / 2, d - off if d < 0 else d + off,
                f"{d:.2f}", ha="center", va=va, fontsize=7.0, color="#2b2b2b")
    ax.set_xticks(range(len(tfs)))
    ax.set_xticklabels(tfs, fontsize=8.0, rotation=30, ha="right")
    ax.set_ylabel(f"Simulated Δ{target} (a.u.)", fontsize=8.5)
    ax.set_title(f"CellOracle: TF knockdown effect on {target}\n(GSE151263 PBMC)",
                 fontsize=9.5, fontweight="bold", loc="left", color="#2b2b2b", pad=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color("#cccccc"); ax.spines["bottom"].set_color("#cccccc")
    ax.grid(axis="y", color="#e8e8e8", lw=0.6, zorder=1)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=7.5)
    fig.subplots_adjust(left=0.13, right=0.96, top=0.84, bottom=0.26)
    save3(fig, f"CellOracle_KO_effect_{target}", CO, 7.0, 4.2)

# 2b. 上游 TF coef 热图（top TF x hub 基因）
cm = list(csv.reader(open(os.path.join(CO, "CellOracle_coef_matrix.csv"))))
header = cm[0][1:]
mat_rows = cm[1:]
tfs_all = [r[0] for r in mat_rows]
mat = np.array([[float(v) for v in r[1:]] for r in mat_rows])

targets = [t for t in eff.keys()]                     # UPP1 / EXOSC4
tcols = [header.index(t) for t in targets]
sub = mat[:, tcols]                                   # TF x target
# 每个 target 取 |coef| top10 的 TF 并集
sel = set()
for j, t in enumerate(targets):
    order = np.argsort(-np.abs(sub[:, j]))[:10]
    sel.update(tfs_all[i] for i in order)
sel_tfs = sorted(sel)
rows_idx = [tfs_all.index(t) for t in sel_tfs]
hm = sub[rows_idx, :]

order = np.argsort(-np.abs(hm).max(axis=1))           # 按最大 |coef| 排序
hm = hm[order]
sel_tfs = [sel_tfs[i] for i in order]

fig, ax = plt.subplots(figsize=(4.6, 0.34 * len(sel_tfs) + 1.6))
vmax = np.abs(hm).max()
im = ax.imshow(hm, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_xticks(range(len(targets)))
ax.set_xticklabels(targets, fontsize=9, style="italic")
ax.set_yticks(range(len(sel_tfs)))
ax.set_yticklabels(sel_tfs, fontsize=7.5)
for i in range(hm.shape[0]):
    for j in range(hm.shape[1]):
        ax.text(j, i, f"{hm[i, j]:.3f}", ha="center", va="center",
                fontsize=6.6, color="white" if abs(hm[i, j]) > vmax * 0.55 else "#2b2b2b")
cb = fig.colorbar(im, ax=ax, shrink=0.7)
cb.set_label("TF -> target coef", fontsize=8)
cb.ax.tick_params(labelsize=7)
ax.set_title("CellOracle: upstream TF coefficients\n(GSE151263 PBMC)",
             fontsize=9.5, fontweight="bold", loc="left", color="#2b2b2b", pad=10)
fig.subplots_adjust(left=0.22, right=0.93, top=0.88, bottom=0.06)
save3(fig, "CellOracle_KO_heatmap", CO, 4.6, 0.34 * len(sel_tfs) + 1.6)

print("done.")
