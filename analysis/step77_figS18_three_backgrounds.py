#!/usr/bin/env python3
"""
Step 77 -- 重绘 Figure S18: 加入背景 C (DEG-matched) 后的三背景阴性对照图
A: best glucocorticoid -log10(P) 分布 (背景 A 候选基因 n=87 命中 + 背景 C DEG-matched n=87 命中)
B: top-10 糖皮质激素命中率四组对比 (全转录组 0% / 候选 17% / DEG-matched 18% / 真实 100%)
输出三格式: pdf + png@300 + tiff@1000 (PIL LZW)
"""
import csv, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
FIGDIR = f"{BASE}/AJRCCM_submission/05_supplementary/Supplementary_Figures_Tables/Figures"
OUT = f"{FIGDIR}/FigureS18_CMap_null_distribution"
REAL = 3.66e-4
CAND, CBLUE, REALC, DEGC = "#8FA2D6", "#56B4E9", "#D95F4E", "#E7A889"

def neglog(p):
    try:
        return -math.log10(float(p))
    except Exception:
        return None

# 背景 A: 候选基因 (Old_CMAP_up)
cand_vals, cand_top10, cand_n = [], 0, 0
with open(f"{BASE}/09_docking_drug/cmap_negative_control_results.csv") as fh:
    for r in csv.DictReader(fh):
        if r["library"] == "Old_CMAP_up" and r["background"] == "A_nonhub_candidates":
            cand_n += 1
            if r["gluco_in_top10"] == "True":
                cand_top10 += 1
            if r["best_gluco_p"]:
                v = neglog(r["best_gluco_p"])
                if v is not None:
                    cand_vals.append(v)

# 背景 C: DEG-matched
deg_vals, deg_top10, deg_n, deg_hits = [], 0, 0, 0
with open(f"{BASE}/09_docking_drug/cmap_deg_matched_null.csv") as fh:
    for r in csv.DictReader(fh):
        deg_n += 1
        if r["gluco_in_top10"] == "True":
            deg_top10 += 1
        if r["best_gluco_p"]:
            deg_hits += 1
            v = neglog(r["best_gluco_p"])
            if v is not None:
                deg_vals.append(v)

real_x = -math.log10(REAL)
print(f"A: 候选 n_hits={len(cand_vals)}/{cand_n}, top10={cand_top10} | C: DEG n_hits={len(deg_vals)}/{deg_n}, top10={deg_top10} | real -log10P={real_x:.2f}")

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), dpi=300)
bins = [x * 0.1 for x in range(10, 60)]
ax = axes[0]
ax.hist(cand_vals, bins=bins, alpha=0.85, color=CAND, label=f"non-hub candidates (n={len(cand_vals)} with hits)")
ax.hist(deg_vals, bins=bins, alpha=0.55, color=DEGC, label=f"DEG-matched (n={len(deg_vals)} with hits)")
ax.axvline(real_x, color=REALC, lw=1.6)
ymax = ax.get_ylim()[1]
ax.text(real_x + 0.08, ymax * 0.92, "six-hub-gene\nsignature", color=REALC, fontsize=7.5, ha="left", va="top")
ax.set_xlabel("best glucocorticoid  −log$_{10}$(P)", fontsize=8.5)
ax.set_ylabel("random gene sets", fontsize=8.5)
ax.set_title("A  Null distributions (Old_CMAP_up)", fontsize=9, loc="left")
ax.legend(fontsize=6.5, frameon=False)
ax.tick_params(labelsize=7.5)
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)

ax = axes[1]
labels = ["random\nmeasured genes", "random non-hub\ncandidates", "DEG-matched\nrandom sets", "six-hub-gene\nsignature"]
vals = [0, cand_top10, deg_top10, 100]
colors = [CAND, CBLUE, DEGC, REALC]
bars = ax.bar(labels, vals, color=colors, width=0.62)
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width() / 2, v + 2, f"{v}%", ha="center", fontsize=7.5)
ax.set_ylabel("glucocorticoid in top 10 (%)", fontsize=8.5)
ax.set_ylim(0, 115)
ax.set_title("B  Top-10 hit frequency", fontsize=9, loc="left")
ax.tick_params(axis="x", labelsize=7, rotation=0)
ax.tick_params(axis="y", labelsize=7.5)
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)

fig.tight_layout()
fig.savefig(f"{OUT}.pdf")
fig.savefig(f"{OUT}.png", dpi=300)
# TIFF @1000 dpi 经 PNG 中转 (PIL LZW)
img = Image.open(f"{OUT}.png")
w, h = img.size
scale = 1000 / 300
img.resize((int(w * scale), int(h * scale)), Image.LANCZOS).save(
    f"{OUT}.tiff", compression="tiff_lzw")
print("三格式导出完成:", f"{OUT}.pdf/.png/.tiff")
