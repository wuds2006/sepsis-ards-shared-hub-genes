#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""仅重拼 Figure 2 composite (F 面板换为双队列 ML+SHAP 2x2)。
与 assemble_figures.py 中 assemble() 同逻辑; 不触碰用户在 'Figure 1/2/3' 中的并行重组。
"""
import os
from PIL import Image, ImageDraw, ImageFont

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
FX = os.path.join(BASE, "AJRCCM_submission", "07_figure_exports")
FONT = "/System/Library/Fonts/Helvetica.ttc"
W_FULL = 3170
PAD = 20
LETTER_BOX = 110


def scaled(img, target_w):
    w, h = img.size
    return img.resize((target_w, int(round(h * target_w / w))), Image.LANCZOS)


def load(src, target_w):
    return scaled(Image.open(src).convert("RGB"), target_w)


def assemble(rows, out_base, out_w=W_FULL):
    colw = {1: out_w, 2: (out_w - PAD) // 2, 3: (out_w - 2 * PAD) // 3,
            4: (out_w - 3 * PAD) // 4}
    row_imgs = []
    for row in rows:
        cells = []
        for imgs, letter, cw in row:
            tw = out_w if cw == FULL else (out_w - (cw - 1) * PAD) // cw
            cells.append((load(imgs, tw) if isinstance(imgs, str)
                          else load_hstack(imgs, tw), letter))
        rh = max(im.size[1] for im, _ in cells)
        rm = Image.new("RGB", (out_w, rh), "white")
        x = 0
        for im, _ in cells:
            rm.paste(im, (x, 0))
            x += im.size[0] + PAD
        row_imgs.append((rm, cells))
    H = sum(r.size[1] for r, _ in row_imgs) + PAD * (len(row_imgs) - 1)
    canvas = Image.new("RGB", (out_w, H), "white")
    d = ImageDraw.Draw(canvas)
    f = ImageFont.truetype(FONT, 80)
    y = 0
    for rm, cells in row_imgs:
        canvas.paste(rm, (0, y))
        x = 0
        for im, letter in cells:
            d.rectangle([x, y, x + LETTER_BOX, y + LETTER_BOX], fill="white")
            d.text((x + 8, y + 4), letter, fill="black", font=f)
            x += im.size[0] + PAD
        y += rm.size[1] + PAD
    canvas.save(out_base + ".png", dpi=(300, 300))
    canvas.save(out_base + ".tiff", dpi=(300, 300), compression="tiff_lzw")
    canvas.save(out_base + ".pdf")
    print("saved", out_base + ".*", canvas.size)


def load_hstack(paths, target_w):
    imgs = [Image.open(p).convert("RGB") for p in paths]
    tw = (target_w - PAD * (len(imgs) - 1)) // len(imgs)
    imgs = [scaled(im, tw) for im in imgs]
    h = max(im.size[1] for im in imgs)
    c = Image.new("RGB", (tw * len(imgs) + PAD * (len(imgs) - 1), h), "white")
    x = 0
    for im in imgs:
        c.paste(im, (x, 0))
        x += im.size[0] + PAD
    return c


FULL = 1
CELL4 = 4
F2E = os.path.join(FX, "Figure2_expression_ROC")
F2 = os.path.join(FX, "Figure2")
assemble([
    [(os.path.join(F2E, "Figure2_boxplot_GSE76293.png"), "A", CELL4),
     (os.path.join(F2E, "Figure2_boxplot_GSE243066.png"), "B", CELL4),
     (os.path.join(F2E, "Figure2_boxplot_GSE69528.png"), "C", CELL4),
     (os.path.join(F2E, "Figure2_boxplot_GSE100159.png"), "D", CELL4)],
    [(os.path.join(F2E, "Figure2_ROC_AUC_CI.png"), "E", FULL)],
    [(os.path.join(F2, "parts", "F_ML_SHAP_dual.png"), "F", FULL)],
], os.path.join(F2, "Figure2_composite"))
