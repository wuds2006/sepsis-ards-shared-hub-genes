#!/usr/bin/env python
"""Export Figure S3 panels C, D, E as standalone figures.

  C  Sensitivity-methods forest (IVW / MR-Egger / weighted median / weighted mode)
  D  SMR-HEIDI results for the five hub genes (5 genes x 2 outcomes)
  E  LDSC SNP heritability + genetic correlation

Each panel saved as vector PDF + PNG/TIFF at 1000 dpi (line-art standard,
also satisfies colour >= 300 dpi). Nature 7-color palette.
"""
import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D

ROOT = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
MR = os.path.join(ROOT, "04_MR_SMR")
OUT = os.path.join(MR, "figS3_parts")
os.makedirs(OUT, exist_ok=True)

PURPLE = "#6A5C95"; PINK = "#D99AAE"; ORANGE = "#E7A889"
GREEN = "#A9CBB8"; MINT = "#CFE3E6"; BLUE = "#8FA2D6"; NAVY = "#5B608C"

for fp in ("/System/Library/Fonts/Helvetica.ttc",):
    if os.path.exists(fp):
        font_manager.fontManager.addfont(fp)
plt.rcParams["font.family"] = ["Helvetica", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = True
plt.rcParams["pdf.fonttype"] = 42  # editable text in vector PDF


def save_all(fig, stem):
    """Save one panel as PDF (vector) + PNG/TIFF at 1000 dpi."""
    fig.savefig(os.path.join(OUT, stem + ".pdf"), facecolor="white")
    fig.savefig(os.path.join(OUT, stem + ".png"), dpi=1000, facecolor="white")
    fig.savefig(os.path.join(OUT, stem + ".tiff"), dpi=1000, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print("saved:", stem)


# ------------------------------------------------------------------ panel C
rows = list(csv.DictReader(open(os.path.join(MR, "MR_sepsis4980_ARDS_results.csv"))))
methods = [
    ("Inverse variance weighted", "IVW"),
    ("MR Egger", "MR-Egger"),
    ("Weighted median", "Weighted median"),
    ("Weighted mode", "Weighted mode"),
]
lab, b, lo, hi = [], [], [], []
for m, short in methods:
    r = next(x for x in rows if x["method"] == m)
    beta = float(r["b"]); se = float(r["se"]); p = float(r["pval"])
    lab.append(f"{short}\nβ = {beta:.3f}\nP = {p:.3f}")
    b.append(beta); lo.append(beta - 1.96 * se); hi.append(beta + 1.96 * se)

het = list(csv.DictReader(open(os.path.join(MR, "MR_sepsis4980_heterogeneity.csv"))))
q_ivw = next(x for x in het if x["method"] == "Inverse variance weighted")
ple = list(csv.DictReader(open(os.path.join(MR, "MR_sepsis4980_pleiotropy.csv"))))[0]

fig, ax = plt.subplots(figsize=(7.0, 4.6))
ys = range(len(lab))
cols = [NAVY, PURPLE, PINK, ORANGE]
for i, y in enumerate(ys):
    ax.plot([lo[i], hi[i]], [y, y], color=cols[i], lw=2.2, zorder=3,
            solid_capstyle="round")
    ax.plot([lo[i], lo[i]], [y - 0.12, y + 0.12], color=cols[i], lw=2.2, zorder=3)
    ax.plot([hi[i], hi[i]], [y - 0.12, y + 0.12], color=cols[i], lw=2.2, zorder=3)
    ax.scatter([b[i]], [y], s=90, color=cols[i], zorder=4, edgecolors="white",
               linewidths=1.2)
ax.axvline(0, color="#b00020", lw=1.0, ls="--", zorder=2)
ax.text(0.03, -0.52, "null (β = 0)", color="#b00020", fontsize=7.0)
ax.set_yticks(list(ys))
ax.set_yticklabels(lab, fontsize=7.6)
ax.set_ylim(-0.6, len(lab) - 0.4)
ax.set_xlabel("Causal estimate β (log odds of ARDS per log odds of sepsis), 95% CI",
              fontsize=7.8)
ax.tick_params(axis="x", labelsize=7.2)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.spines["left"].set_color("#cccccc"); ax.spines["bottom"].set_color("#cccccc")
ax.grid(axis="x", color="#e8e8e8", lw=0.6, zorder=1)
ax.set_axisbelow(True)
ax.set_title("C  Sensitivity analyses (all null)", fontsize=9.0, fontweight="bold",
             loc="left", color="#2b2b2b", pad=10)
fig.text(0.26, 0.055,
        f"Heterogeneity (IVW Cochran Q): Q = {float(q_ivw['Q']):.2f}, "
        f"df = {q_ivw['Q_df']}, P = {float(q_ivw['Q_pval']):.3f}    |    "
        f"MR-Egger intercept: {float(ple['egger_intercept']):.4f} "
        f"(P = {float(ple['pval']):.3f})\n"
        "Steiger filtering: correct causal direction "
        "(P = 2.65 × 10$^{-34}$)",
        ha="left", va="center", fontsize=6.6, color="#555555", linespacing=1.6)
fig.subplots_adjust(left=0.26, right=0.97, top=0.90, bottom=0.24)
save_all(fig, "FigureS3C_sensitivity_forest")

# ------------------------------------------------------------------ panel D
gene_ensg = {
    "CDK5RAP2": "ENSG00000136861", "CD177": "ENSG00000204936",
    "GRB10": "ENSG00000106070", "UPP1": "ENSG00000183696",
    "VNN1": "ENSG00000112299",
}
smr = []          # (gene, outcome, p_smr, p_heidi)
for g, ensg in gene_ensg.items():
    for out in ("sepsis", "ards"):
        fn = os.path.join(MR, f"smr_{g}_{out}_heidi.smr")
        for r in csv.DictReader(open(fn), delimiter="\t"):
            if r["probeID"] == ensg:
                smr.append((g, out, float(r["p_SMR"]), float(r["p_HEIDI"])))

fig, ax = plt.subplots(figsize=(7.0, 4.6))
genes = list(gene_ensg.keys())
ypos = {g: i for i, g in enumerate(genes)}
thr = -math.log10(0.05)
ax.axvline(thr, color="#b00020", lw=1.0, ls="--", zorder=2)
ax.text(thr + 0.03, 4.42, "P = 0.05", color="#b00020", fontsize=7.0,
        rotation=90, va="top")
for g, out, ps, ph in smr:
    y = ypos[g] + (0.18 if out == "sepsis" else -0.18)
    c = PINK if out == "sepsis" else BLUE
    ax.scatter([-math.log10(ps)], [y], marker="o", s=72, color=c, zorder=4,
               edgecolors="white", linewidths=1.0)
    ax.scatter([-math.log10(ph)], [y], marker="s", s=62, facecolors="none",
               edgecolors=c, linewidths=1.8, zorder=4)
ax.set_yticks(range(len(genes)))
ax.set_yticklabels(genes, fontsize=8.4, style="italic")
ax.set_xlabel("−log10(P value)", fontsize=7.8)
ax.set_xlim(0, 1.5)
ax.set_ylim(-0.6, 4.6)
ax.tick_params(axis="x", labelsize=7.2)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.spines["bottom"].set_color("#cccccc")
ax.grid(axis="x", color="#e8e8e8", lw=0.6, zorder=1)
ax.set_axisbelow(True)
ax.set_title("D  SMR-HEIDI: no causality, no pleiotropy (10 tests)",
             fontsize=9.0, fontweight="bold", loc="left", color="#2b2b2b",
             pad=10)
handles = [
    Line2D([0], [0], marker="o", color="none", markerfacecolor=PINK,
           markeredgecolor="white", markersize=8.5,
           label="p_SMR (sepsis)"),
    Line2D([0], [0], marker="s", color="none", markerfacecolor="none",
           markeredgecolor=PINK, markeredgewidth=1.8, markersize=8,
           label="p_HEIDI (sepsis)"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=BLUE,
           markeredgecolor="white", markersize=8.5,
           label="p_SMR (ARDS)"),
    Line2D([0], [0], marker="s", color="none", markerfacecolor="none",
           markeredgecolor=BLUE, markeredgewidth=1.8, markersize=8,
           label="p_HEIDI (ARDS)"),
]
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.16),
          fontsize=6.8, frameon=False, ncol=2, handletextpad=0.2,
          columnspacing=1.2)
fig.subplots_adjust(left=0.15, right=0.97, top=0.90, bottom=0.24)
save_all(fig, "FigureS3D_SMR_HEIDI")

# ------------------------------------------------------------------ panel E
fig, axes = plt.subplots(1, 2, figsize=(7.0, 4.6),
                         gridspec_kw={"wspace": 0.45})
# left: SNP heritability
ax = axes[0]
traits = ["Sepsis\n(UKB)", "ARDS\n(FinnGen)"]
h2 = [0.0022, 0.0019]; h2se = [0.0007, 0.0006]
xs = [0, 1]
ax.bar(xs, h2, width=0.55, color=[PINK, BLUE], zorder=3)
ax.errorbar(xs, h2, yerr=h2se, fmt="none", ecolor="#555555", elinewidth=1.4,
            capsize=5, zorder=4)
for x, v, se in zip(xs, h2, h2se):
    ax.text(x, v + se + 0.00012, f"h² = {v:.4f}\n(SE {se:.4f})", ha="center",
            fontsize=7.0, color="#2b2b2b")
ax.set_xticks(xs); ax.set_xticklabels(traits, fontsize=7.8)
ax.set_ylabel("Observed-scale SNP heritability (h²)", fontsize=7.8)
ax.set_ylim(0, 0.0035)
ax.tick_params(labelsize=7.2)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.spines["left"].set_color("#cccccc"); ax.spines["bottom"].set_color("#cccccc")
ax.set_title("E  LDSC: weak germline signal", fontsize=9.0, fontweight="bold",
             loc="left", color="#2b2b2b", pad=10)
# right: genetic correlation
ax = axes[1]
rg, rgse = -0.0515, 0.245
ax.errorbar([0], [rg], yerr=[[rgse], [rgse]], fmt="o", color=NAVY, markersize=9,
            elinewidth=1.6, capsize=6, zorder=4)
ax.axhline(0, color="#b00020", lw=1.0, ls="--", zorder=2)
ax.text(0.28, 0.05, "rg = 0", color="#b00020", fontsize=7.2)
ax.text(0.02, rg - rgse - 0.09,
        "rg = −0.05\n(SE 0.25)\nP = 0.83", fontsize=7.6, color="#2b2b2b")
ax.set_xlim(-0.55, 0.55); ax.set_ylim(-0.42, 0.42)
ax.set_xticks([])
ax.set_ylabel("Genetic correlation (rg)", fontsize=7.8)
ax.tick_params(labelsize=7.2)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.spines["left"].set_color("#cccccc"); ax.spines["bottom"].set_visible(False)
ax.set_title("   LDSC: no shared genetic basis", fontsize=9.0,
             fontweight="bold", loc="left", color="#2b2b2b", pad=10)
fig.subplots_adjust(left=0.11, right=0.97, top=0.88, bottom=0.10)
save_all(fig, "FigureS3E_LDSC")
