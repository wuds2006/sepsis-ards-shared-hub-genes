#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 61 —— 仅重建"表达箱线图 + ROC + ML AUC/SHAP"composite
（当前编号 = Figure 3；step39 内部旧命名 Figure 2，不可全量重跑以免覆盖新编号文件）
布局与 step39 Figure 2 节完全一致。产出 composites/Figure2_composite.pdf
（保持旧命名惯例）→ 同步 Figures_final/Figure3_composite.pdf
"""
import pymupdf, shutil, os

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission"
P = f"{BASE}/03_figures/Figures_ggplot2_NPG/main_panels"
OUT = f"{BASE}/03_figures/Figures_ggplot2_NPG/composites"
FIN = f"{BASE}/03_figures/Figures_final"

MM = 72 / 25.4
PAGE_W = 178 * MM
ROW_GAP = 6.0
COL_GAP = 4.0
MARGIN_X = 8.0
TOP = 14.0
X1 = 8.0
X2 = 8.0 + (178 - 2 * MARGIN_X) / 2 + COL_GAP / 2
THIRD = (PAGE_W - 2 * MARGIN_X - 2 * COL_GAP) / 3
XA = 8.0
XB = XA + THIRD + COL_GAP
XC = XB + THIRD + COL_GAP

def mm(v): return v * MM

def pdf_h(path, w_mm):
    d = pymupdf.open(path); r = d[0].rect; d.close()
    return w_mm * r.height / r.width

class Sheet:
    def __init__(self, height_mm):
        self.doc = pymupdf.open()
        self.page = self.doc.new_page(width=PAGE_W, height=mm(height_mm))
        self.cache = {}
    def _pdf(self, path):
        if path not in self.cache:
            self.cache[path] = pymupdf.open(path)
        return self.cache[path]
    def letter(self, x_mm, y_mm, text):
        self.page.insert_text((mm(x_mm), mm(y_mm) + 3.4), text,
                              fontsize=9, fontname="hebo", color=(0, 0, 0))
    def pdf_panel(self, path, x_mm, y_mm, w_mm, letter=None):
        d = self._pdf(path)
        r = d[0].rect
        h = w_mm * r.height / r.width
        self.page.show_pdf_page(pymupdf.Rect(mm(x_mm), mm(y_mm),
                                             mm(x_mm) + mm(w_mm), mm(y_mm) + h), d, 0)
        if letter:
            self.letter(x_mm - 5.5, y_mm, letter)
        return h
    def save(self, name):
        out = os.path.join(OUT, name)
        self.doc.save(out)
        self.doc.close()
        print("saved:", out)

w3 = THIRD
hb = pdf_h(f"{P}/Fig2A_boxplot_GSE76293.pdf", w3)
hr = pdf_h(f"{P}/Fig2E1_ROC_GSE76293.pdf", w3)
wF = (178 - 2 * MARGIN_X - COL_GAP) / 2 - 6
hF = pdf_h(f"{P}/Fig2F1_ML_AUC_barplot.pdf", wF)
H = TOP + hb + ROW_GAP + max(hb, hr) + ROW_GAP + hr + ROW_GAP + hF + 6
s = Sheet(H)
y = TOP
s.pdf_panel(f"{P}/Fig2A_boxplot_GSE76293.pdf", XA, y, w3, "A")
s.pdf_panel(f"{P}/Fig2B_boxplot_GSE243066.pdf", XB, y, w3, "B")
s.pdf_panel(f"{P}/Fig2C_boxplot_GSE69528.pdf", XC, y, w3, "C")
y += hb + ROW_GAP
s.pdf_panel(f"{P}/Fig2D_boxplot_GSE100159.pdf", XA, y, w3, "D")
s.pdf_panel(f"{P}/Fig2E1_ROC_GSE76293.pdf", XB, y, w3, "E")
s.pdf_panel(f"{P}/Fig2E2_ROC_GSE243066.pdf", XC, y, w3)
y += max(hb, hr) + ROW_GAP
s.pdf_panel(f"{P}/Fig2E3_ROC_GSE69528.pdf", XB, y, w3)
s.pdf_panel(f"{P}/Fig2E4_ROC_GSE100159.pdf", XC, y, w3)
y += hr + ROW_GAP
s.pdf_panel(f"{P}/Fig2F1_ML_AUC_barplot.pdf", X1, y, wF, "F")
s.pdf_panel(f"{P}/Fig2F2_SHAP_barplot.pdf", X2, y, wF)
s.save("Figure2_composite.pdf")

shutil.copy2(os.path.join(OUT, "Figure2_composite.pdf"),
             os.path.join(FIN, "Figure3_composite.pdf"))
print("synced to Figures_final/Figure3_composite.pdf")
