#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 80 —— 重建 Figure3_composite（表达箱线图 + ROC + ML AUC/SHAP）
修复 step61 的坐标单位混用 bug（mm/pt 混算导致 B/C/E 面板被放到页面外）。
布局（与图注一致）：
  Row1: A Fig2A_boxplot_GSE76293 | B Fig2B_boxplot_GSE243066 | C Fig2C_boxplot_GSE69528
  Row2: D Fig2D_boxplot_GSE100159 | E Fig2E1_ROC_GSE76293 | Fig2E2_ROC_GSE243066
  Row3:                           |   Fig2E3_ROC_GSE69528   | Fig2E4_ROC_GSE100159
  Row4: F Fig2F1_ML_AUC_barplot   | Fig2F2_SHAP_barplot     (2 栏)
产出：Figures_final/Figure3_composite.{pdf,png(300dpi),tiff(1000dpi LZW)}
"""
import os
import pymupdf
from PIL import Image

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission"
P = f"{BASE}/03_figures/Figures_ggplot2_NPG/main_panels"
FIN = f"{BASE}/03_figures/Figures_final"

MM = 72 / 25.4                      # points per mm
PAGE_W = 178 * MM                   # double-column width
MARGIN = 8 * MM
COL_GAP = 4 * MM
ROW_GAP = 6 * MM
TOP = 14 * MM
THIRD = (PAGE_W - 2 * MARGIN - 2 * COL_GAP) / 3
HALF = (PAGE_W - 2 * MARGIN - COL_GAP) / 2
XA = MARGIN
XB = XA + THIRD + COL_GAP
XC = XB + THIRD + COL_GAP

def panel_h(name, w_pt):
    d = pymupdf.open(f"{P}/{name}.pdf")
    r = d[0].rect
    d.close()
    return w_pt * r.height / r.width

h_box = panel_h("Fig2A_boxplot_GSE76293", THIRD)
h_roc = panel_h("Fig2E1_ROC_GSE76293", THIRD)
h_f1 = panel_h("Fig2F1_ML_AUC_barplot", HALF)
h_f2 = panel_h("Fig2F2_SHAP_barplot", HALF)
h_f = max(h_f1, h_f2)

PAGE_H = TOP + h_box + ROW_GAP + max(h_box, h_roc) + ROW_GAP + h_roc + ROW_GAP + h_f + 4 * MM

doc = pymupdf.open()
page = doc.new_page(width=PAGE_W, height=PAGE_H)
cache = {}

def place(name, x, y, w, letter=None):
    if name not in cache:
        cache[name] = pymupdf.open(f"{P}/{name}.pdf")
    d = cache[name]
    r = d[0].rect
    h = w * r.height / r.width
    page.show_pdf_page(pymupdf.Rect(x, y, x + w, y + h), d, 0)
    if letter:
        page.insert_text((x - 5.5 * MM, y + 3.4), letter,
                         fontsize=9, fontname="hebo", color=(0, 0, 0))
    return h

y = TOP
place("Fig2A_boxplot_GSE76293", XA, y, THIRD, "A")
place("Fig2B_boxplot_GSE243066", XB, y, THIRD, "B")
place("Fig2C_boxplot_GSE69528", XC, y, THIRD, "C")
y += h_box + ROW_GAP
place("Fig2D_boxplot_GSE100159", XA, y, THIRD, "D")
place("Fig2E1_ROC_GSE76293", XB, y, THIRD, "E")
place("Fig2E2_ROC_GSE243066", XC, y, THIRD)
y += max(h_box, h_roc) + ROW_GAP
place("Fig2E3_ROC_GSE69528", XB, y, THIRD)
place("Fig2E4_ROC_GSE100159", XC, y, THIRD)
y += h_roc + ROW_GAP
place("Fig2F1_ML_AUC_barplot", MARGIN, y, HALF, "F")
place("Fig2F2_SHAP_barplot", MARGIN + HALF + COL_GAP, y, HALF)

out_pdf = f"{FIN}/Figure3_composite.pdf"
doc.save(out_pdf, deflate=True)
doc.close()
print("saved:", out_pdf, f"({PAGE_W/72:.2f} x {PAGE_H/72:.2f} in)")

# raster outputs
d = pymupdf.open(out_pdf)
pg = d[0]
pix = pg.get_pixmap(matrix=pymupdf.Matrix(300 / 72, 300 / 72), alpha=False)
img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
img.save(f"{FIN}/Figure3_composite.png", dpi=(300, 300))
pix = pg.get_pixmap(matrix=pymupdf.Matrix(1000 / 72, 1000 / 72), alpha=False)
img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
img.save(f"{FIN}/Figure3_composite.tiff", format="TIFF",
         compression="tiff_lzw", dpi=(1000, 1000))
d.close()
for ext in ["pdf", "png", "tiff"]:
    print(ext, round(os.path.getsize(f"{FIN}/Figure3_composite.{ext}") / 1e6, 2), "MB")
