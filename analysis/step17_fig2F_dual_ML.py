#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Figure 2F 双队列 ML+SHAP 面板 (2x2):
上行: ARDS (GSE76293) AUC 条形图 | Sepsis (GSE69528) AUC 条形图
下行: ARDS SHAP 蜂群图     | Sepsis SHAP 蜂群图
输出: Figure2/parts/F_ML_SHAP_dual.png (300dpi 等效)
"""
from PIL import Image, ImageDraw, ImageFont
import os

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
SRC = os.path.join(BASE, "05_ML_SHAP")
OUT = os.path.join(BASE, "AJRCCM_submission", "07_figure_exports", "Figure2", "parts")
os.makedirs(OUT, exist_ok=True)

FONT = "/System/Library/Fonts/Helvetica.ttc"

def load(p):
    return Image.open(p).convert("RGB")

# 行 1: 两个 AUC 条形图 (2100x1500); 行 2: 两个 SHAP 蜂群图 (2400x2100)
r1 = [load(os.path.join(SRC, "ML_AUC_barplot.png")),
      load(os.path.join(SRC, "sepsis", "ML_AUC_barplot_sepsis.png"))]
r2 = [load(os.path.join(SRC, "SHAP_beeplot.png")),
      load(os.path.join(SRC, "sepsis", "SHAP_beeplot_sepsis.png"))]

# 统一缩放: 每行内两图等高; 行宽合计对齐
W1 = 2100  # 每列目标宽 (AUC 列)
# 行 1: AUC 图 2100x1500 -> 目标 1500 宽
r1 = [im.resize((1500, int(im.size[1] * 1500 / im.size[0])), Image.LANCZOS) for im in r1]
# 行 2: SHAP 图 2400x2100 -> 目标 1500 宽
r2 = [im.resize((1500, int(im.size[1] * 1500 / im.size[0])), Image.LANCZOS) for im in r2]

GAP = 40
row1_w = sum(im.size[0] for im in r1) + GAP
row2_w = sum(im.size[0] for im in r2) + GAP
W = max(row1_w, row2_w)
h1 = max(im.size[1] for im in r1)
h2 = max(im.size[1] for im in r2)
LABEL_H = 90  # 队列标签条
H = h1 + h2 + LABEL_H * 2 + GAP * 3

canvas = Image.new("RGB", (W, H), "white")
d = ImageDraw.Draw(canvas)
f_label = ImageFont.truetype(FONT, 56)

def paste_row(imgs, y0, label_left, label_right):
    # 左标签 (起始 x=150 避开 composite 字母框 ~110px)
    d.text((150, y0), label_left, fill="#2b2b2b", font=f_label)
    d.text((W - 40 - d.textlength(label_right, font=f_label), y0), label_right,
           fill="#2b2b2b", font=f_label)
    y = y0 + LABEL_H
    x = 0
    for im in imgs:
        canvas.paste(im, (x, y))
        x += im.size[0] + GAP
    return y + max(im.size[1] for im in imgs) + GAP

y = 40
y = paste_row(r1, y, "ARDS training cohort (GSE76293)", "Sepsis training cohort (GSE69528)")
paste_row(r2, y, "ARDS training cohort (GSE76293)", "Sepsis training cohort (GSE69528)")

out_png = os.path.join(OUT, "F_ML_SHAP_dual.png")
canvas.save(out_png, dpi=(300, 300))
print("saved", out_png, canvas.size)
