#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 71 —— 重建 Figure 9 composite（交叉学科专家小问题：加 2D 相互作用图）
布局（在 step58 基础上加 D 面板）:
  A: docking site heatmap (PDF 矢量)
  B: PyMOL poses (PNG): VNN1-luteolin + GRB10-EGCG
  C: glucocorticoid repositioning (PDF 矢量)
  D: VNN1-luteolin 2D 相互作用图 (PLIP, PDF 矢量)
产出: composites/Figure9_composite.pdf → Figures_final/ (pdf + png@300 + tiff@1000)
"""
import pymupdf, shutil, os, io
from PIL import Image

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission"
P = f"{BASE}/03_figures/Figures_ggplot2_NPG/main_panels"
DK = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/09_docking_drug"
OUT = f"{BASE}/03_figures/Figures_ggplot2_NPG/composites"
FIN = f"{BASE}/03_figures/Figures_final"

MM = 72 / 25.4
PAGE_W = 178 * MM
ROW_GAP = 6.0
X1 = 8.0
TOP = 14.0

def mm(v): return v * MM

def pdf_h(path, w_mm):
    d = pymupdf.open(path); r = d[0].rect; d.close()
    return w_mm * r.height / r.width

def png_h(path, w_mm):
    pix = pymupdf.Pixmap(path)
    return w_mm * pix.height / pix.width

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
    def png_panel(self, path, x_mm, y_mm, w_mm, letter=None):
        pix = pymupdf.Pixmap(path)
        h = w_mm * pix.height / pix.width
        self.page.insert_image(pymupdf.Rect(mm(x_mm), mm(y_mm),
                                            mm(x_mm) + mm(w_mm), mm(y_mm) + h), pixmap=pix)
        if letter:
            self.letter(x_mm - 5.5, y_mm, letter)
        return h
    def save(self, name):
        out = os.path.join(OUT, name)
        self.doc.save(out)
        self.doc.close()
        print("saved:", out)
        return out

hDk = pdf_h(f"{P}/Fig8A_docking_site_heatmap.pdf", 100)
p1 = f"{DK}/poses_pymol/VNN1_luteolin_pymol.png"
p2 = f"{DK}/poses_pymol/GRB10_EGCG_pymol.png"
hPose = png_h(p1, 31)
hGc = pdf_h(f"{P}/Fig8C_glucocorticoid_repositioning.pdf", 92)
hPl = pdf_h(f"{DK}/plip_vnn1_luteolin_2D.pdf", 76)
H = TOP + max(hDk, hPose) + ROW_GAP + max(hGc, hPl) + 6
s = Sheet(H)
y = TOP
s.pdf_panel(f"{P}/Fig8A_docking_site_heatmap.pdf", X1, y, 100, "A")
s.png_panel(p1, X1 + 104, y, 31, "B")
s.png_panel(p2, X1 + 137, y, 31)
y += max(hDk, hPose) + ROW_GAP
s.pdf_panel(f"{P}/Fig8C_glucocorticoid_repositioning.pdf", X1, y, 92, "C")
s.pdf_panel(f"{DK}/plip_vnn1_luteolin_2D.pdf", X1 + 96, y, 76, "D")
out = s.save("Figure9_composite.pdf")

# 三格式导出到 Figures_final
shutil.copy2(out, os.path.join(FIN, "Figure9_composite.pdf"))
doc = pymupdf.open(out)
pix = doc[0].get_pixmap(dpi=300)
pix.save(os.path.join(FIN, "Figure9_composite.png"))
pix = doc[0].get_pixmap(dpi=1000)
im = Image.open(io.BytesIO(pix.tobytes("png")))
im.save(os.path.join(FIN, "Figure9_composite.tiff"), format="TIFF",
        compression="tiff_lzw", dpi=(1000, 1000))
print("exported pdf + png@300 + tiff@1000 to Figures_final")
