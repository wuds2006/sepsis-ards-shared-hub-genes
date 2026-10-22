#!/usr/bin/env python
"""Assemble Figure 6 (A-H) composite, 4-row layout:
Row 1: A  PBMC UMAP strip (full width)
Row 2: B | C | D  PBMC CellChat | slingshot | CellOracle
Row 3: E  CD14+ monocyte UMAP strip (full width)
Row 4: F | G | H  CD14+ CellChat | slingshot+AUCell | CellOracle
"""
from PIL import Image, ImageDraw, ImageFont
import os

ROOT = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
F6 = os.path.join(ROOT, "AJRCCM_submission/07_figure_exports/Figure6")
SC = os.path.join(ROOT, "08_singlecell/GSE167363")

CELL_W = 1050
GAP = 50
MARGIN = 60
LETTER_SZ = 80
FULL_W = CELL_W * 3 + GAP * 2

A = os.path.join(F6, "A_UMAP.png")
B = os.path.join(F6, "B_cellchat_diff.png")
C = os.path.join(F6, "C_slingshot.png")
D = os.path.join(F6, "parts/D_celloracle_CEBPB.png")
E = os.path.join(SC, "GSE167363_UMAP.png")
F = os.path.join(SC, "cellchat/diff_heatmap.png")
G = [os.path.join(SC, "slingshot/GSE167363_slingshot_gene_trend.png"),
     os.path.join(SC, "aucell/GSE167363_AUCell_violin.png")]
H = os.path.join(SC, "celloracle/CellOracle_KO_heatmap.png")

font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", LETTER_SZ)

def scaled(img, target_w):
    w, h = img.size
    return img.resize((target_w, int(round(h * target_w / w))), Image.LANCZOS)

def load_cell(src, target_w):
    if isinstance(src, list):
        sub = [scaled(Image.open(p).convert("RGB"), target_w) for p in src]
        h = sum(i.size[1] for i in sub) + 30
        canvas = Image.new("RGB", (target_w, h), "white")
        y = 0
        for i in sub:
            canvas.paste(i, (0, y)); y += i.size[1] + 30
        return canvas
    return scaled(Image.open(src).convert("RGB"), target_w)

rows = []  # list of (cells_with_letters, height)
row1 = [(load_cell(A, FULL_W), "A")]
row2 = [(load_cell(s, CELL_W), l) for s, l in zip([B, C, D], "BCD")]
row3 = [(load_cell(E, FULL_W), "E")]
row4 = [(load_cell(s, CELL_W), l) for s, l in zip([F, G, H], "FGH")]
rows = [(row1, row1[0][0].size[1]), (row2, max(c.size[1] for c, _ in row2)),
        (row3, row3[0][0].size[1]), (row4, max(c.size[1] for c, _ in row4))]

# pad cells within each row to common height
for idx in (1, 3):
    row = rows[idx][0]
    h = max(c.size[1] for c, _ in row)
    padded = []
    for c, l in row:
        if c.size[1] < h:
            pad = Image.new("RGB", (c.size[0], h), "white")
            pad.paste(c, (0, (h - c.size[1]) // 2))
            c = pad
        padded.append((c, l))
    rows[idx] = (padded, h)

W = MARGIN * 2 + FULL_W
Hh = MARGIN * 2 + sum(h for _, h in rows) + GAP * (len(rows) - 1)
canvas = Image.new("RGB", (W, Hh), "white")
draw = ImageDraw.Draw(canvas)

y = MARGIN
for cells, h in rows:
    x = MARGIN
    for c, l in cells:
        canvas.paste(c, (x, y))
        draw.rectangle([x, y, x + 105, y + 105], fill="white")
        draw.text((x + 12, y + 6), l, fill="black", font=font)
        x += c.size[0] + GAP
    y += h + GAP

for ext, kw in [("png", dict(dpi=(300, 300))),
                ("tiff", dict(compression="tiff_lzw", dpi=(300, 300))),
                ("pdf", dict(resolution=300))]:
    out = os.path.join(F6, f"Figure6_composite.{ext}")
    canvas.save(out, **kw)
    print("saved", out, canvas.size)
