#!/usr/bin/env python
"""Standalone export of Figure 7 panel C
"Drug repositioning: glucocorticoids prioritized via CMAP / GEO perturbations".

Reproduces the PIL panel (scripts/make_missing_panels.py, Figure 7C) as a
vector PDF plus 1000-dpi PNG/TIFF (flat-color bar chart = line art, so the
>=1000 dpi rule applies). Data identical to the composite panel.

Outputs (AJRCCM_submission/07_figure_exports/Figure7/parts/):
  C_cmap_glucocorticoids.pdf  (vector)
  C_cmap_glucocorticoids.png  (1000 dpi)
  C_cmap_glucocorticoids.tiff (1000 dpi)
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUTDIR = os.path.join(ROOT, "AJRCCM_submission/07_figure_exports/Figure7/parts")
os.makedirs(OUTDIR, exist_ok=True)

# Nature 7-color palette (project convention)
PURPLE = "#6A5C95"; PINK = "#D99AAE"; ORANGE = "#E7A889"

# Helvetica (matches PIL panel); fall back to DejaVu Sans if unavailable
for fp in ("/System/Library/Fonts/Helvetica.ttc",):
    if os.path.exists(fp):
        font_manager.fontManager.addfont(fp)
plt.rcParams["font.family"] = ["Helvetica", "DejaVu Sans"]

hits = [
    ("Prednisolone",  4.433148837263066e-05, "GEO perturbations (rank 1)", "VNN1, EXOSC4, GRB10, UPP1"),
    ("Dexamethasone", 3.664321680484968e-04, "CMAP up (rank 1)",           "GRB10, CDK5RAP2"),
    ("Hydrocortisone", 4.417388600957901e-03, "GEO perturbations (rank 4)", "VNN1, GRB10, UPP1"),
]
colors = [PURPLE, PINK, ORANGE]

neglog = [-math.log10(p) for _, p, _, _ in hits]
maxv = max(neglog) * 1.18            # same headroom as PIL panel
thr = -math.log10(0.05)

# PIL panel: 1800 x 1620 px @300 dpi -> 6.0 x 5.4 in
fig, ax = plt.subplots(figsize=(6.0, 5.4), dpi=1000)

ys = [2, 1, 0]                       # top-to-bottom order
bar_h = 0.42
for (name, p, src, genes), v, c, y in zip(hits, neglog, colors, ys):
    ax.barh(y, v, height=bar_h, color=c, zorder=3)
    ax.text(v + maxv * 0.018, y, f"P = {p:.1e}",
            va="center", ha="left", fontsize=8.2, fontweight="bold",
            color="#2b2b2b", zorder=4)
    ax.text(-maxv * 0.015, y + 0.155, name,
            va="center", ha="right", fontsize=8.4, fontweight="bold",
            color="#2b2b2b")
    ax.text(-maxv * 0.015, y - 0.02, src,
            va="center", ha="right", fontsize=6.6, color="#888888")
    ax.text(-maxv * 0.015, y - 0.195, f"overlap: {genes}",
            va="center", ha="right", fontsize=6.0, color="#888888")

# significance threshold (P = 0.05)
ax.axvline(thr, color="#b00020", lw=1.0, zorder=2)
ax.text(thr, 2.52, "P = 0.05", color="#b00020", fontsize=7.0,
        fontweight="bold", ha="center", va="bottom")

# gridlines + x axis
ax.set_xlim(0, maxv)
ax.set_ylim(-0.55, 2.75)
ax.set_xticks([1, 2, 3, 4, 5])
ax.grid(axis="x", color="#e0e0e0", lw=0.7, zorder=1)
ax.set_axisbelow(True)
ax.set_xlabel("-log10(P value)", fontsize=8.4, color="#2b2b2b")
ax.tick_params(axis="x", labelsize=7.2, colors="#666666", length=0)
ax.set_yticks([])
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.spines["bottom"].set_color("#cccccc")

fig.subplots_adjust(left=0.30, right=0.965, top=0.86, bottom=0.12)
fig.text(0.02, 0.955,
         "Drug repositioning: glucocorticoids prioritized via CMAP / GEO perturbations",
         fontsize=8.8, fontweight="bold", color="#2b2b2b", ha="left", va="top")

base = os.path.join(OUTDIR, "C_cmap_glucocorticoids")
fig.savefig(base + ".pdf", facecolor="white")    # vector
fig.savefig(base + ".png", dpi=1000, facecolor="white")
fig.savefig(base + ".tiff", dpi=1000, facecolor="white",
            pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)

# flatten PNG/TIFF to RGB on white (project convention: no alpha channel)
from PIL import Image
for ext, kw in ((".png", {"dpi": (1000, 1000)}),
                (".tiff", {"dpi": (1000, 1000), "compression": "tiff_lzw"})):
    fp = base + ext
    im = Image.open(fp)
    if im.mode == "RGBA":
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(im, mask=im.split()[3])
        bg.save(fp, **kw)
        im.close(); bg.close()
print("saved:", base + ".pdf/.png/.tiff")
