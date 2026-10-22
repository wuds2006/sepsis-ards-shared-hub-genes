#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 84 —— Figure 2 面板字母 B/D/C -> 常规 A/B/C：
  A = 火山图（上）
  B = Venn（中，居中）
  C = WGCNA（下）
同时同步 Figure_legends.md 中 Figure 2 的面板描述顺序。
"""
import os
import pymupdf
from PIL import Image

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission"
FIN = f"{BASE}/03_figures/Figures_final"
P = f"{BASE}/03_figures/Figures_ggplot2_NPG/main_panels"
MM = 72 / 25.4
PAGE_W = 178 * MM
MARGIN = 23.0
GAP = 12.0
CW = PAGE_W - 2 * MARGIN

def ph(path, w):
    d = pymupdf.open(path); r = d[0].rect; d.close()
    return w * r.height / r.width

panels = {
    "A1": f"{P}/Fig1B1_volcano_GSE76293.pdf",
    "A2": f"{P}/Fig1B2_volcano_GSE69528.pdf",
    "C1": f"{P}/Fig1C1_module_trait_GSE76293.pdf",
    "C2": f"{P}/Fig1C2_module_trait_GSE69528.pdf",
    "B":  f"{P}/Fig1D_venn_4sets.pdf",
}
half = (CW - GAP) / 2
hA = ph(panels["A1"], half)
hC = ph(panels["C1"], half)
B_W = 300.0
hB = ph(panels["B"], B_W)
TOP = 36.0
PAGE_H = TOP + hA + GAP + hB + GAP + hC + 18

doc = pymupdf.open()
page = doc.new_page(width=PAGE_W, height=PAGE_H)
def letter(x, y, ch):
    page.insert_text((x - 12, y + 3.4), ch, fontsize=9,
                     fontname="hebo", color=(0, 0, 0))
# Row 1: A1 A2 (volcano)
y = TOP
page.show_pdf_page(pymupdf.Rect(MARGIN, y, MARGIN + half, y + hA),
                   pymupdf.open(panels["A1"]), 0)
page.show_pdf_page(pymupdf.Rect(MARGIN + half + GAP, y, PAGE_W - MARGIN, y + hA),
                   pymupdf.open(panels["A2"]), 0)
letter(MARGIN, y, "A")
# Row 2: B centered (Venn)
y += hA + GAP
xB = (PAGE_W - B_W) / 2
page.show_pdf_page(pymupdf.Rect(xB, y, xB + B_W, y + hB),
                   pymupdf.open(panels["B"]), 0)
letter(xB, y, "B")
# Row 3: C1 C2 (WGCNA)
y += hB + GAP
page.show_pdf_page(pymupdf.Rect(MARGIN, y, MARGIN + half, y + hC),
                   pymupdf.open(panels["C1"]), 0)
page.show_pdf_page(pymupdf.Rect(MARGIN + half + GAP, y, PAGE_W - MARGIN, y + hC),
                   pymupdf.open(panels["C2"]), 0)
letter(MARGIN, y, "C")
doc.save(f"{FIN}/Figure2_composite.pdf", deflate=True)
doc.close()
print(f"Figure2 (A/B/C): {PAGE_W/72:.2f} x {PAGE_H/72:.2f} in | "
      f"A {half:.0f}x{hA:.0f} B {B_W:.0f}x{hB:.0f} C {half:.0f}x{hC:.0f}")

# raster: PNG @300 + TIFF @1000 LZW
d = pymupdf.open(f"{FIN}/Figure2_composite.pdf")
pg = d[0]
pix = pg.get_pixmap(matrix=pymupdf.Matrix(300 / 72, 300 / 72), alpha=False)
Image.frombytes("RGB", (pix.width, pix.height), pix.samples).save(
    f"{FIN}/Figure2_composite.png", dpi=(300, 300))
pix = pg.get_pixmap(matrix=pymupdf.Matrix(1000 / 72, 1000 / 72), alpha=False)
img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
img.save(f"{FIN}/Figure2_composite.tiff", format="TIFF",
         compression="tiff_lzw", dpi=(1000, 1000))
d.close()
for ext in ["pdf", "png", "tiff"]:
    print(f"  Figure2.{ext}", round(os.path.getsize(f"{FIN}/Figure2_composite.{ext}") / 1e6, 2), "MB")

# ---- update Figure_legends.md panel letters for Figure 2 ----
LG = f"{BASE}/03_figures/Figure_legends.md"
lg = open(LG, encoding="utf-8").read()
old = ("(B) Volcano plots of differentially expressed genes (DEGs) in GSE76293 "
       "(ARDS vs. control; 557 DEGs: 291 upregulated, 266 downregulated) and "
       "GSE69528 (sepsis vs. control; 1,458 DEGs: 658 upregulated, 800 downregulated). "
       "Red/pink points indicate upregulated genes and blue points indicate downregulated "
       "genes (|log2 fold change| > 1 and adjusted P < 0.05). (D) Four-set Venn diagram "
       "showing the intersection of DEGs and WGCNA module genes, yielding 179 shared "
       "candidate genes. (C) Weighted gene co-expression network analysis (WGCNA) "
       "module\u2013trait relationships in the two discovery cohorts, identifying 13 (ARDS) "
       "and 11 (sepsis) disease-associated modules.")
new = ("(A) Volcano plots of differentially expressed genes (DEGs) in GSE76293 "
       "(ARDS vs. control; 557 DEGs: 291 upregulated, 266 downregulated) and "
       "GSE69528 (sepsis vs. control; 1,458 DEGs: 658 upregulated, 800 downregulated). "
       "Red/pink points indicate upregulated genes and blue points indicate downregulated "
       "genes (|log2 fold change| > 1 and adjusted P < 0.05). (B) Four-set Venn diagram "
       "showing the intersection of DEGs and WGCNA module genes, yielding 179 shared "
       "candidate genes. (C) Weighted gene co-expression network analysis (WGCNA) "
       "module\u2013trait relationships in the two discovery cohorts, identifying 13 (ARDS) "
       "and 11 (sepsis) disease-associated modules.")
assert old in lg, "old legend text not found"
lg = lg.replace(old, new)
open(LG, "w", encoding="utf-8").write(lg)
print("Figure_legends.md Figure 2 letters updated: B/D/C -> A/B/C")
print("done")
