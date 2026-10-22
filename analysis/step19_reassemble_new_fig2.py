#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""为用户新图序的「Figure 2」(富集) 目录拼装 composite: 2x2
A 179CGs GO sankey | B 179CGs KEGG chord
C 6hub GO sankey   | D 6hub KEGG chord
输出: Figure 2/Figure2_composite.{pdf,png,tiff} 300dpi
"""
import os
from PIL import Image, ImageDraw, ImageFont

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
F2 = os.path.join(BASE, "AJRCCM_submission", "07_figure_exports", "Figure 2")
FONT = "/System/Library/Fonts/Helvetica.ttc"
OUT_W = 3170
GAP = 40


def load(p, tw):
    im = Image.open(p)
    if im.mode in ("RGBA", "LA", "P"):
        im = im.convert("RGBA")
        bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
        im = Image.alpha_composite(bg, im).convert("RGB")
    else:
        im = im.convert("RGB")
    w, h = im.size
    return im.resize((tw, int(round(h * tw / w))), Image.LANCZOS)


panels = [
    ("179CGs_GO_sankey_bubble.png", "A"),
    ("179CGs_KEGG_chord.png", "B"),
    ("6hub_GO_sankey_bubble.png", "C"),
    ("6hub_KEGG_chord.png", "D"),
]
cw = (OUT_W - GAP) // 2
imgs = [(load(os.path.join(F2, p), cw), l) for p, l in panels]
r1h = max(im.size[1] for im, _ in imgs[:2])
r2h = max(im.size[1] for im, _ in imgs[2:])
H = r1h + r2h + GAP
canvas = Image.new("RGB", (OUT_W, H), "white")
d = ImageDraw.Draw(canvas)
f = ImageFont.truetype(FONT, 80)
positions = [(0, 0), (cw + GAP, 0), (0, r1h + GAP), (cw + GAP, r1h + GAP)]
for (im, letter), (x, y) in zip(imgs, positions):
    canvas.paste(im, (x, y))
    d.rectangle([x, y, x + 105, y + 105], fill="white")
    d.text((x + 12, y + 6), letter, fill="black", font=f)

stem = os.path.join(F2, "Figure2_composite")
canvas.save(stem + ".png", dpi=(300, 300))
canvas.save(stem + ".tiff", dpi=(300, 300), compression="tiff_lzw")
canvas.save(stem + ".pdf")
print("saved", stem + ".*", canvas.size)
