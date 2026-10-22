#!/usr/bin/env python
"""Assemble Figure 8 composite:
Row 1: A ARDS nomogram | B Sepsis nomogram
Row 2: C ARDS ROC + calibration + DCA strip
Row 3: D Sepsis ROC + calibration + DCA strip
"""
from PIL import Image, ImageChops, ImageDraw, ImageFont
import os

ROOT = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
F8 = os.path.join(ROOT, "AJRCCM_submission/07_figure_exports/Figure8")

CELL_W = 1500          # nomogram cell width
GAP = 50
MARGIN = 60
LETTER_SZ = 80
FULL_W = CELL_W * 2 + GAP

font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", LETTER_SZ)

def trim(im, pad=12):
    diff = ImageChops.difference(im, Image.new("RGB", im.size, (255, 255, 255)))
    bbox = diff.getbbox()
    if bbox:
        im = im.crop((max(0, bbox[0] - pad), max(0, bbox[1] - pad),
                      min(im.size[0], bbox[2] + pad), min(im.size[1], bbox[3] + pad)))
    return im

def scaled(img, target_w):
    w, h = img.size
    return img.resize((target_w, int(round(h * target_w / w))), Image.LANCZOS)

def load(p, w, do_trim=False):
    im = Image.open(p).convert("RGB")
    if do_trim:
        im = trim(im)
    return scaled(im, w)

row1 = [(load(os.path.join(F8, "ARDS_nomogram.png"), CELL_W, do_trim=True), "A"),
        (load(os.path.join(F8, "Sepsis_nomogram.png"), CELL_W, do_trim=True), "B")]
row2 = [(load(os.path.join(F8, "ARDS_ROC_cal_DCA.png"), FULL_W), "C")]
row3 = [(load(os.path.join(F8, "Sepsis_ROC_cal_DCA.png"), FULL_W), "D")]

rows = []
for row in (row1, row2, row3):
    h = max(c.size[1] for c, _ in row)
    padded = []
    for c, l in row:
        if c.size[1] < h:
            pad = Image.new("RGB", (c.size[0], h), "white")
            pad.paste(c, (0, (h - c.size[1]) // 2))
            c = pad
        padded.append((c, l))
    rows.append((padded, h))

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
    out = os.path.join(F8, f"Figure8_composite.{ext}")
    canvas.save(out, **kw)
    print("saved", out, canvas.size)
