#!/usr/bin/env python
"""Figure S3 -- Mendelian randomization shows no independent causal effect.

Panels:
  A  MR scatter plot (sepsis UKB ieu-b-4980 -> ARDS FinnGen J10), existing
  B  MR forest plot (per-SNP Wald ratios + IVW), existing
  C  Sensitivity-methods forest (IVW / MR-Egger / weighted median / weighted mode),
     with heterogeneity and pleiotropy annotations   [new]
  D  SMR-HEIDI results for the five hub genes with cis-eQTL (5 genes x 2 outcomes)
     [new, values parsed from smr_*_heidi.smr]
  E  LDSC genetic correlation and SNP heritability  [new, values from ldsc logs]

Composite exported as FigureS3_MR_validation.{pdf,png,tiff}.
Nature 7-color palette. Composite is color -> PNG/TIFF >= 300 dpi (500 used).
"""
import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from PIL import Image, ImageDraw, ImageFont

ROOT = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
MR = os.path.join(ROOT, "04_MR_SMR")
PARTS = os.path.join(MR, "figS3_parts")
os.makedirs(PARTS, exist_ok=True)

PURPLE = "#6A5C95"; PINK = "#D99AAE"; ORANGE = "#E7A889"
GREEN = "#A9CBB8"; MINT = "#CFE3E6"; BLUE = "#8FA2D6"; NAVY = "#5B608C"

for fp in ("/System/Library/Fonts/Helvetica.ttc",):
    if os.path.exists(fp):
        font_manager.fontManager.addfont(fp)
plt.rcParams["font.family"] = ["Helvetica", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = True

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

fig, ax = plt.subplots(figsize=(7.0, 4.6), dpi=500)
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
ax.set_title("Sensitivity analyses (all null)", fontsize=9.0, fontweight="bold",
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
fig.savefig(os.path.join(PARTS, "C_sensitivity_forest.png"), dpi=500,
            facecolor="white")
plt.close(fig)

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

fig, ax = plt.subplots(figsize=(7.0, 4.6), dpi=500)
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
ax.set_title("SMR-HEIDI: no causality, no pleiotropy (10 tests)", fontsize=9.0,
             fontweight="bold", loc="left", color="#2b2b2b", pad=10)
from matplotlib.lines import Line2D
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
fig.savefig(os.path.join(PARTS, "D_smr_heidi.png"), dpi=500, facecolor="white")
plt.close(fig)

# ------------------------------------------------------------------ panel E
fig, axes = plt.subplots(1, 2, figsize=(7.0, 4.6), dpi=500,
                         gridspec_kw={"wspace": 0.45})
# left: SNP heritability
ax = axes[0]
traits = ["Sepsis\n(UKB)", "ARDS\n(FinnGen)"]
h2 = [0.0022, 0.0019]; h2se = [0.0007, 0.0006]
xs = [0, 1]
bars = ax.bar(xs, h2, width=0.55, color=[PINK, BLUE], zorder=3)
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
ax.set_title("LDSC: weak germline signal", fontsize=9.0, fontweight="bold",
             loc="left", color="#2b2b2b", pad=10)
# right: genetic correlation
ax = axes[1]
rg, rgse = -0.0515, 0.245
ax.errorbar([0], [rg], yerr=[[rgse], [rgse]], fmt="o", color=NAVY, markersize=9,
            elinewidth=1.6, capsize=6, zorder=4)
ax.axhline(0, color="#b00020", lw=1.0, ls="--", zorder=2)
ax.text(0.28, 0.05, "rg = 0", color="#b00020", fontsize=7.2)
ax.text(0.02, rg - rgse - 0.09,
        f"rg = −0.05\n(SE 0.25)\nP = 0.83", fontsize=7.6, color="#2b2b2b")
ax.set_xlim(-0.55, 0.55); ax.set_ylim(-0.42, 0.42)
ax.set_xticks([])
ax.set_ylabel("Genetic correlation (rg)", fontsize=7.8)
ax.tick_params(labelsize=7.2)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.spines["left"].set_color("#cccccc"); ax.spines["bottom"].set_visible(False)
ax.set_title("LDSC: no shared genetic basis", fontsize=9.0, fontweight="bold",
             loc="left", color="#2b2b2b", pad=10)
fig.subplots_adjust(left=0.11, right=0.97, top=0.88, bottom=0.10)
fig.savefig(os.path.join(PARTS, "E_ldsc.png"), dpi=500, facecolor="white")
plt.close(fig)

# ------------------------------------------------------------------ composite
def load_rgb(path, target_w):
    im = Image.open(path)
    if im.mode != "RGB":
        bg = Image.new("RGB", im.size, (255, 255, 255))
        if im.mode == "RGBA":
            bg.paste(im, mask=im.split()[3])
        else:
            bg.paste(im.convert("RGB"))
        im = bg
    w, h = im.size
    return im.resize((target_w, round(h * target_w / w)), Image.LANCZOS)

A = load_rgb(os.path.join(MR, "MR_sepsis4980_scatter.png"), 1550)
B = load_rgb(os.path.join(MR, "MR_sepsis4980_forest.png"), 1550)
C = load_rgb(os.path.join(PARTS, "C_sensitivity_forest.png"), 1010)
D = load_rgb(os.path.join(PARTS, "D_smr_heidi.png"), 1010)
E = load_rgb(os.path.join(PARTS, "E_ldsc.png"), 1010)

PAD, GAP, LBL = 170, 30, 150
row1_h = max(A.size[1], B.size[1])
row2_h = max(C.size[1], D.size[1], E.size[1])
W = PAD * 2 + A.size[0] + GAP + B.size[0]
H = LBL + PAD * 3 + row1_h + GAP + row2_h
canvas = Image.new("RGB", (W, H), "white")
draw = ImageDraw.Draw(canvas)
HF = "/System/Library/Fonts/Helvetica.ttc"
fB = ImageFont.truetype(HF, 96)

def letter(x, y, ch):
    draw.rectangle([x, y, x + 120, y + 120], fill="#2b2b2b")
    bb = draw.textbbox((0, 0), ch, font=fB)
    draw.text((x + (120 - (bb[2] - bb[0])) // 2 - bb[0],
               y + (120 - (bb[3] - bb[1])) // 2 - bb[1]), ch, fill="white", font=fB)

# row 1
y1 = LBL + PAD
canvas.paste(A, (PAD, y1)); canvas.paste(B, (PAD + A.size[0] + GAP, y1))
letter(PAD - 130, y1 + 40, "A")
letter(PAD + A.size[0] + GAP - 130, y1 + 40, "B")
# row 2
y2 = y1 + row1_h + GAP + PAD
w2 = C.size[0] + D.size[0] + E.size[0] + GAP * 2
x0 = (W - w2) // 2
canvas.paste(C, (x0, y2))
canvas.paste(D, (x0 + C.size[0] + GAP, y2))
canvas.paste(E, (x0 + C.size[0] + D.size[0] + GAP * 2, y2))
letter(x0 - 130, y2 + 40, "C")
letter(x0 + C.size[0] + GAP - 130, y2 + 40, "D")
letter(x0 + C.size[0] + D.size[0] + GAP * 2 - 130, y2 + 40, "E")

out = os.path.join(MR, "FigureS3_MR_validation")
canvas.save(out + ".png", dpi=(500, 500))
canvas.save(out + ".tiff", dpi=(500, 500), compression="tiff_lzw")
canvas.save(out + ".pdf", resolution=500)
print("saved:", out + ".{png,tiff,pdf}", canvas.size)
