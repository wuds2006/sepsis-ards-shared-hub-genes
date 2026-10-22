#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 79 - 重建当前 Figure 2 composite（PPI + GO dotplot + sankey/chord 六面板，
即 step39 内部旧 "Figure 3" 节布局），并同步重栅格化 PNG@300 / TIFF@1000 LZW。
原因：Fig3A_PPI_network 面板字体修复（Times-Bold -> Helvetica 斜体基因名）。
布局参数与 step39 Figure 3 节逐行一致。"""
import pymupdf, os
from PIL import Image

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission"
P = f"{BASE}/03_figures/Figures_ggplot2_NPG/main_panels"
FIN = f"{BASE}/03_figures/Figures_final"

MM = 72 / 25.4
PAGE_W = 178 * MM
MARGIN_X = 8.0
ROW_GAP = 6.0
COL_GAP = 4.0
X1 = MARGIN_X
TOP = 8.0

def mm(v): return v * MM

def pdf_h(path, w_mm):
    d = pymupdf.open(path)
    r = d[0].rect
    d.close()
    return w_mm * r.height / r.width

hA3 = pdf_h(f"{P}/Fig3A_PPI_network.pdf", 84)
hB3 = pdf_h(f"{P}/Fig3B_GO_BP_dotplot.pdf", 74)
hC3 = pdf_h(f"{P}/Fig3C_GO_sankey_bubble_179.pdf", 100)
hD3 = pdf_h(f"{P}/Fig3D_KEGG_chord_179.pdf", 60)
hE3 = pdf_h(f"{P}/Fig3E_GO_sankey_bubble_6hub.pdf", 96)
hF3 = pdf_h(f"{P}/Fig3F_KEGG_chord_6hub.pdf", 60)
H = TOP + max(hA3, hB3) + ROW_GAP + max(hC3, hD3) + ROW_GAP + max(hE3, hF3) + 6

doc = pymupdf.open()
page = doc.new_page(width=PAGE_W, height=mm(H))
cache = {}

def panel(path, x_mm, y_mm, w_mm, letter=None):
    if path not in cache:
        cache[path] = pymupdf.open(path)
    src = cache[path]
    r = src[0].rect
    h_mm = w_mm * r.height / r.width
    dest = pymupdf.Rect(mm(x_mm), mm(y_mm), mm(x_mm + w_mm), mm(y_mm + h_mm))
    page.show_pdf_page(dest, src, 0)
    if letter:
        page.insert_text((mm(x_mm - 4.5), mm(y_mm) + 3.4), letter,
                         fontsize=9, fontname="hebo", color=(0, 0, 0))
    return h_mm

y = TOP
panel(f"{P}/Fig3A_PPI_network.pdf", X1, y, 84, "A")
panel(f"{P}/Fig3B_GO_BP_dotplot.pdf", X1 + 88, y, 74, "B")
y += max(hA3, hB3) + ROW_GAP
panel(f"{P}/Fig3C_GO_sankey_bubble_179.pdf", X1, y, 100, "C")
panel(f"{P}/Fig3D_KEGG_chord_179.pdf", X1 + 104, y, 60, "D")
y += max(hC3, hD3) + ROW_GAP
panel(f"{P}/Fig3E_GO_sankey_bubble_6hub.pdf", X1, y, 96, "E")
panel(f"{P}/Fig3F_KEGG_chord_6hub.pdf", X1 + 100, y, 60, "F")

out_pdf = f"{FIN}/Figure2_composite.pdf"
doc.save(out_pdf, deflate=True, garbage=3)
for d in cache.values():
    d.close()
doc.close()
print("saved", out_pdf)

# 栅格化：PNG @300 dpi + TIFF @1000 dpi LZW（与 rasterize_final.py 同一管线）
d = pymupdf.open(out_pdf)
pg = d[0]
print(f"page size: {pg.rect.width/72:.2f} x {pg.rect.height/72:.2f} in")
pix = pg.get_pixmap(matrix=pymupdf.Matrix(300 / 72, 300 / 72), alpha=False)
pix.save(f"{FIN}/Figure2_composite.png")
pix2 = pg.get_pixmap(matrix=pymupdf.Matrix(1000 / 72, 1000 / 72), alpha=False)
img = Image.frombytes("RGB", (pix2.width, pix2.height), pix2.samples)
img.save(f"{FIN}/Figure2_composite.tiff", format="TIFF", compression="tiff_lzw", dpi=(1000, 1000))
d.close()
print(f"png {pix.width}x{pix.height}px @300 | tiff {pix2.width}x{pix2.height}px @1000 "
      f"{os.path.getsize(f'{FIN}/Figure2_composite.tiff')/1e6:.1f}MB")
