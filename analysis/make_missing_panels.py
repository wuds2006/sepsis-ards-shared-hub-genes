#!/usr/bin/env python
"""Generate the two missing main-figure panels:
- Figure1/parts/A_workflow_schematic.png : study workflow schematic (Nature 7 palette)
- Figure7/parts/C_cmap_glucocorticoids.png : CMAP/GEO glucocorticoid bar chart
"""
import os
import csv
from PIL import Image, ImageDraw, ImageFont

ROOT = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
F1 = os.path.join(ROOT, "AJRCCM_submission/07_figure_exports/Figure1")
F7 = os.path.join(ROOT, "AJRCCM_submission/07_figure_exports/Figure7")
os.makedirs(os.path.join(F1, "parts"), exist_ok=True)
os.makedirs(os.path.join(F7, "parts"), exist_ok=True)

# Nature 7-color palette (project convention)
PURPLE = "#6A5C95"; PINK = "#D99AAE"; ORANGE = "#E7A889"
GREEN = "#A9CBB8"; MINT = "#CFE3E6"; BLUE = "#8FA2D6"; NAVY = "#5B608C"

HF = "/System/Library/Fonts/Helvetica.ttc"
def font(sz, idx=0):
    return ImageFont.truetype(HF, sz, index=idx)

# ---------------------------------------------------------------- Figure 1A
W, H = 3000, 1230
img = Image.new("RGB", (W, H), "white")
d = ImageDraw.Draw(img)

def box(x0, y0, x1, y1, fill, outline, lines, txt_color="#2b2b2b", fs=44, oblique_idx=None):
    d.rounded_rectangle([x0, y0, x1, y1], radius=22, fill=fill, outline=outline, width=5)
    f = font(fs, oblique_idx if oblique_idx is not None else 0)
    total_h = len(lines) * (fs + 12) - 12
    ty = y0 + (y1 - y0 - total_h) // 2
    for ln in lines:
        bb = d.textbbox((0, 0), ln, font=f)
        tw = bb[2] - bb[0]
        d.text((x0 + (x1 - x0 - tw) // 2, ty), ln, fill=txt_color, font=f)
        ty += fs + 12

def arrow(x0, y0, x1, y1, color="#8a8a8a", w=8):
    d.line([x0, y0, x1, y1], fill=color, width=w)
    import math
    ang = math.atan2(y1 - y0, x1 - x0)
    L = 26
    for sgn in (1, -1):
        a = ang + sgn * 0.5
        d.line([x1, y1, x1 - L * math.cos(a - 1.5707 * 0) , y1 - L * math.sin(a)], fill=color, width=w)

# column x-centres
CX = [430, 1200, 1970, 2670]
BW = 560  # box width

# Stage 1: discovery cohorts
box(CX[0]-BW//2, 120, CX[0]+BW//2, 300, MINT, BLUE,
    ["GSE76293: ARDS", "12 vs 12 (PMN)"])
box(CX[0]-BW//2, 380, CX[0]+BW//2, 560, MINT, BLUE,
    ["GSE69528: Sepsis", "83 vs 28 (whole blood)"])
d.text((CX[0]-BW//2, 40), "Discovery cohorts", fill="#2b2b2b", font=font(48, 1))

# Stage 2: DEG + WGCNA
box(CX[1]-BW//2, 120, CX[1]+BW//2, 300, "#ffffff", PURPLE,
    ["limma DEGs", "557 (ARDS) / 1,458 (sepsis)"])
box(CX[1]-BW//2, 380, CX[1]+BW//2, 560, "#ffffff", PURPLE,
    ["WGCNA modules", "13 (ARDS) / 11 (sepsis)"])
d.text((CX[1]-BW//2, 40), "DEG + WGCNA", fill="#2b2b2b", font=font(48, 1))

# Stage 3: intersection -> enrichment
box(CX[2]-BW//2, 120, CX[2]+BW//2, 300, PINK, "#b06d85",
    ["179 shared", "candidate genes"])
box(CX[2]-BW//2, 380, CX[2]+BW//2, 560, "#ffffff", PINK,
    ["GO / KEGG", "enrichment (clusterProfiler)"])
d.text((CX[2]-BW//2, 40), "Intersection", fill="#2b2b2b", font=font(48, 1))

# Stage 4: ML + hub genes
box(CX[3]-BW//2, 120, CX[3]+BW//2, 300, "#ffffff", ORANGE,
    ["11 ML classifiers", "+ SHAP"])
box(CX[3]-BW//2, 380, CX[3]+BW//2, 560, ORANGE, "#c07a55",
    ["6 hub genes"])
d.text((CX[3]-BW//2, 40), "Prioritization", fill="#2b2b2b", font=font(48, 1))

# horizontal arrows between stages
for i in range(3):
    arrow(CX[i]+BW//2+18, 210, CX[i+1]-BW//2-18, 210)
    arrow(CX[i]+BW//2+18, 470, CX[i+1]-BW//2-18, 470)

# hub gene names band (italic, Helvetica-Oblique index 2)
fy = 640
d.rounded_rectangle([150, fy, W-150, fy+110], radius=22, fill="#f4f1f8", outline=PURPLE, width=4)
genes = "CDK5RAP2   CD177   GRB10   UPP1   EXOSC4   VNN1"
fg = font(52, 2)
bb = d.textbbox((0, 0), genes, font=fg)
d.text(((W - (bb[2]-bb[0])) // 2, fy + 24), genes, fill=PURPLE, font=fg)

# downstream applications row
apps = [
    ("External validation", "GSE243066 / GSE100159", BLUE),
    ("Immune microenv.", "M2 / complement (GSVA)", GREEN),
    ("Single-cell 2 cohorts", "CellChat / slingshot / KO", PINK),
    ("MR / SMR-HEIDI", "causal assessment", ORANGE),
    ("Docking + repositioning", "Vina / CMAP", NAVY),
    ("Nomogram (rms/lrm)", "ROC / calibration / DCA", PURPLE),
]
ay0, ah = 830, 210
bw2 = (W - 300 - 5 * 40) // 6
for i, (t1, t2, c) in enumerate(apps):
    x0 = 150 + i * (bw2 + 40)
    box(x0, ay0, x0 + bw2, ay0 + ah, "#ffffff", c, [t1, t2], fs=36)
# arrow from hub band down to applications
arrow(W // 2, fy + 118, W // 2, ay0 - 10)

# caption strip
d.text((150, ay0 + ah + 50), "Multi-omics integration pipeline: discovery -> prioritization -> multi-level validation and translation",
       fill="#666666", font=font(36))

img.save(os.path.join(F1, "parts/A_workflow_schematic.png"), dpi=(300, 300))
print("saved Figure1/parts/A_workflow_schematic.png", img.size)

# ---------------------------------------------------------------- Figure 7C
hits = [
    ("Prednisolone", 4.433148837263066e-05, "GEO perturbations (rank 1)", "VNN1, EXOSC4, GRB10, UPP1"),
    ("Dexamethasone", 3.664321680484968e-04, "CMAP up (rank 1)", "GRB10, CDK5RAP2"),
    ("Hydrocortisone", 4.417388600957901e-03, "GEO perturbations (rank 4)", "VNN1, GRB10, UPP1"),
]
W2, H2 = 1800, 1620
im2 = Image.new("RGB", (W2, H2), "white")
d2 = ImageDraw.Draw(im2)
TOP_PAD = 120  # room for the composite panel-letter box
d2.text((60, TOP_PAD - 20), "Drug repositioning: glucocorticoids prioritized via CMAP / GEO perturbations",
        fill="#2b2b2b", font=font(44, 1))

import math
neglog = [-math.log10(p) for _, p, _, _ in hits]
maxv = max(neglog) * 1.18
x0, y0, bw3 = 470, TOP_PAD + 120, 380
gap = 170
plot_w = W2 - x0 - 140
bar_h = 110
colors = [PURPLE, PINK, ORANGE]
thr = -math.log10(0.05)
# axis gridlines
for v in [1, 2, 3, 4, 5]:
    gx = x0 + plot_w * v / maxv
    d2.line([gx, y0 - 30, gx, y0 + 3 * (bar_h + gap)], fill="#e0e0e0", width=3)
    d2.text((gx - 14, y0 + 3 * (bar_h + gap) + 14), str(v), fill="#666666", font=font(36))
for i, ((name, p, src, genes2), v, c) in enumerate(zip(hits, neglog, colors)):
    y = y0 + i * (bar_h + gap)
    w = plot_w * v / maxv
    d2.rounded_rectangle([x0, y, x0 + w, y + bar_h], radius=10, fill=c)
    d2.text((x0 + w + 24, y + bar_h // 2 - 26), f"P = {p:.1e}", fill="#2b2b2b", font=font(40, 1))
    d2.text((60, y + 6), name, fill="#2b2b2b", font=font(42, 1))
    d2.text((60, y + 60), src, fill="#888888", font=font(32))
    d2.text((60, y + bar_h + 12), f"overlap: {genes2}", fill="#888888", font=font(28))
# significance threshold line (P = 0.05)
tx = x0 + plot_w * thr / maxv
d2.line([tx, y0 - 30, tx, y0 + 3 * (bar_h + gap)], fill="#b00020", width=4)
d2.text((tx - 130, y0 - 80), "P = 0.05", fill="#b00020", font=font(34, 1))
d2.text((x0 + plot_w // 2 - 200, y0 + 3 * (bar_h + gap) + 64), "-log10(P value)", fill="#2b2b2b", font=font(40))

im2.save(os.path.join(F7, "parts/C_cmap_glucocorticoids.png"), dpi=(300, 300))
print("saved Figure7/parts/C_cmap_glucocorticoids.png", im2.size)
