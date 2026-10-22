#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 65 (subject-reviewer revision) — rebuild Figure 4 and Figure 8 composites
after panel updates (step63 nomogram simplification, step64 group-color contrast).

Layout is byte-for-byte the same as step39_assemble_composites.py:
  Figure 4: A/B nomograms side by side; C/D ROC+calibration+DCA strips below.
  Figure 8: single-cell 8-panel grid (step39 internal "Figure 7" section).

Outputs:
  composites/Figure4_composite.pdf, composites/Figure8_composite.pdf
  Figures_final/{Figure4,Figure8}_composite.{pdf,tiff,png}  (TIFF 1000 dpi, PNG 300 dpi)
"""
import pymupdf, shutil, os

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission"
P = f"{BASE}/03_figures/Figures_ggplot2_NPG/main_panels"
SC = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/08_singlecell"
OUT = f"{BASE}/03_figures/Figures_ggplot2_NPG/composites"
FIN = f"{BASE}/03_figures/Figures_final"

MM = 72 / 25.4
PAGE_W = 178 * MM
ROW_GAP = 6.0
X1 = 8.0
TOP = 14.0
CW = (178 - X1 * 2 - 4) / 2
X2 = X1 + CW + 4

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
        src = self._pdf(path)
        p = src[0]
        h = w_mm * p.rect.height / p.rect.width
        self.page.show_pdf_page(
            pymupdf.Rect(mm(x_mm), mm(y_mm), mm(x_mm + w_mm), mm(y_mm + h)), src)
        if letter:
            self.letter(x_mm - 1.5, y_mm, letter)
        return h
    def save(self, name):
        out = os.path.join(OUT, name)
        self.doc.save(out)
        self.doc.close()
        print("saved:", out)

def export_final(name):
    import io
    from PIL import Image
    src = os.path.join(OUT, name + ".pdf")
    shutil.copy2(src, os.path.join(FIN, name + ".pdf"))
    d = pymupdf.open(src)
    pix = d[0].get_pixmap(dpi=300)
    pix.save(os.path.join(FIN, name + ".png"))
    pix = d[0].get_pixmap(dpi=1000)
    im = Image.open(io.BytesIO(pix.tobytes("png")))
    im.save(os.path.join(FIN, name + ".tiff"), format="TIFF", compression="tiff_lzw", dpi=(1000, 1000))
    d.close()
    print("Figures_final:", name, "(pdf + png300 + tiff1000)")

## ---------------- Figure 4 (identical layout to step39) ----------------
hN = pdf_h(f"{P}/Fig4A_nomogram_ARDS.pdf", CW)
hS = pdf_h(f"{P}/Fig4C_ROC_cal_DCA_ARDS.pdf", 162)
H = TOP + hN + ROW_GAP + hS + ROW_GAP + hS + 6
s = Sheet(H)
y = TOP
s.pdf_panel(f"{P}/Fig4A_nomogram_ARDS.pdf", X1, y, CW, "A")
s.pdf_panel(f"{P}/Fig4B_nomogram_Sepsis.pdf", X2, y, CW, "B")
y += hN + ROW_GAP
s.pdf_panel(f"{P}/Fig4C_ROC_cal_DCA_ARDS.pdf", X1, y, 162, "C")
y += hS + ROW_GAP
s.pdf_panel(f"{P}/Fig4D_ROC_cal_DCA_sepsis.pdf", X1, y, 162, "D")
s.save("Figure4_composite.pdf")

## ---------------- Figure 8 (single-cell; identical layout to step39 "Figure 7" section) ----------------
hU = pdf_h(f"{P}/Fig7A_UMAP_GSE151263.pdf", 76)
hCC = pdf_h(f"{SC}/cellchat/diff_heatmap.pdf", 84)
hT = pdf_h(f"{P}/Fig7C_trajectory_AUCell_GSE151263.pdf", 162)
hO = pdf_h(f"{P}/Fig7D_CellOracle_GSE151263.pdf", 84)
H = (TOP + max(hU, hCC) + ROW_GAP + hT + ROW_GAP + hO + ROW_GAP
     + max(hU, hCC) + ROW_GAP + hT + ROW_GAP + hO + 6)
s = Sheet(H)
y = TOP
s.pdf_panel(f"{P}/Fig7A_UMAP_GSE151263.pdf", X1, y, 76, "A")
s.pdf_panel(f"{SC}/cellchat/diff_heatmap.pdf", X1 + 80, y, 84, "B")
y += max(hU, hCC) + ROW_GAP
s.pdf_panel(f"{P}/Fig7C_trajectory_AUCell_GSE151263.pdf", X1, y, 162, "C")
y += hT + ROW_GAP
s.pdf_panel(f"{P}/Fig7D_CellOracle_GSE151263.pdf", X1, y, 84, "D")
y += hO + ROW_GAP
s.pdf_panel(f"{P}/Fig7E_UMAP_GSE167363.pdf", X1, y, 76, "E")
s.pdf_panel(f"{SC}/GSE167363/cellchat/diff_heatmap.pdf", X1 + 80, y, 84, "F")
y += max(hU, hCC) + ROW_GAP
s.pdf_panel(f"{P}/Fig7G_trajectory_AUCell_GSE167363.pdf", X1, y, 162, "G")
y += hT + ROW_GAP
s.pdf_panel(f"{P}/Fig7H_CellOracle_GSE167363.pdf", X1, y, 84, "H")
s.save("Figure8_composite.pdf")

export_final("Figure4_composite")
export_final("Figure8_composite")
print("Step 65 done.")
