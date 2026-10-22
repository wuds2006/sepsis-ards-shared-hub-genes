#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 81 —— 重建 Figure9_composite（分子对接 + 药物重定位）
修复原版布局缺陷：A/B/C/D 比例失调、大面积留白、B 构象图过小（31×23 pt）。
布局（178 mm 双栏宽，字号可读性优先）：
  左列 (332 pt)：A 对接热图 (332×147, 字号 5.0 pt)
                 C 糖皮质激素重定位条形图 (332×181, 字号 8.0 pt)
  右列 (119 pt)：B1 VNN1-luteolin 构象 (119×89)
                 B2 GRB10-EGCG 构象 (119×89)，各居其行
  底行：D PLIP 二维相互作用图 (360×288, 字号 5.0 pt) 居中
产出：Figures_final/Figure9_composite.{pdf,png(300dpi),tiff(1000dpi LZW)}
"""
import os
import pymupdf
from PIL import Image

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
P = f"{BASE}/AJRCCM_submission/03_figures/Figures_ggplot2_NPG/main_panels"
FIN = f"{BASE}/AJRCCM_submission/03_figures/Figures_final"
POSES = f"{BASE}/09_docking_drug/poses_pymol"

MM = 72 / 25.4
PAGE_W = 178 * MM
MARGIN = 7 * MM
TOP = 12 * MM
GAP = 12
LEFT_W = 332
RIGHT_W = PAGE_W - 2 * MARGIN - LEFT_W - GAP   # ~119 pt
D_W = 360

# --- flatten RGBA pose PNGs onto white RGB (keeps TIFF/PDF clean) ---
flat = {}
for key, src in [("B1", "VNN1_luteolin_pymol.png"), ("B2", "GRB10_EGCG_pymol.png")]:
    p = f"{POSES}/{src}"
    im = Image.open(p)
    if im.mode != "RGB":
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(im, mask=im.split()[-1] if im.mode in ("RGBA", "LA") else None)
        im = bg
    out = f"{POSES}/{src.replace('.png', '_flat.png')}"
    im.save(out, dpi=(300, 300))
    flat[key] = out

def panel_h(path, w):
    d = pymupdf.open(path)
    r = d[0].rect
    d.close()
    return w * r.height / r.width

A_PDF = f"{P}/Fig8A_docking_site_heatmap.pdf"
C_PDF = f"{P}/Fig8C_glucocorticoid_repositioning.pdf"
D_PDF = f"{BASE}/09_docking_drug/plip_vnn1_luteolin_2D.pdf"

hA = panel_h(A_PDF, LEFT_W)
hC = panel_h(C_PDF, LEFT_W)
hB = RIGHT_W * 1200 / 1600.0        # 4:3 pose images
hD = panel_h(D_PDF, D_W)

LEFT_COL_H = hA + GAP + hC
PAGE_H = TOP + LEFT_COL_H + 14 + hD + 10 * MM

doc = pymupdf.open()
page = doc.new_page(width=PAGE_W, height=PAGE_H)

def letter(x, y, ch):
    page.insert_text((x - 12, y + 3.4), ch, fontsize=9,
                     fontname="hebo", color=(0, 0, 0))

def caption(x, y_center, italic_part, roman_part, size=6.2):
    """One-line caption: italic gene/ligand part + roman part, centered."""
    wi = pymupdf.get_text_length(italic_part, fontname="heit", fontsize=size)
    wr = pymupdf.get_text_length(roman_part, fontname="helv", fontsize=size)
    x0 = x - (wi + wr) / 2
    page.insert_text((x0, y_center), italic_part, fontsize=size,
                     fontname="heit", color=(0, 0, 0))
    page.insert_text((x0 + wi, y_center), roman_part, fontsize=size,
                     fontname="helv", color=(0, 0, 0))

XA = MARGIN
XC = MARGIN + LEFT_W + GAP
XR = XA + LEFT_W + GAP

# --- Row 1 left: A heatmap ---
page.show_pdf_page(pymupdf.Rect(XA, TOP, XA + LEFT_W, TOP + hA),
                   pymupdf.open(A_PDF), 0)
letter(XA, TOP, "A")

# --- Row 1 right: B1 pose (centered in row height) ---
yB1 = TOP + (hA - hB) / 2
img = page.insert_image(pymupdf.Rect(XR, yB1, XR + RIGHT_W, yB1 + hB),
                        filename=flat["B1"])
letter(XR, yB1, "B")
cap_y = yB1 + hB + 8
caption(XR + RIGHT_W / 2, cap_y, "VNN1", "-luteolin (-10.01 kcal/mol)")

# --- Row 2 left: C barplot ---
yC = TOP + hA + GAP
page.show_pdf_page(pymupdf.Rect(XA, yC, XA + LEFT_W, yC + hC),
                   pymupdf.open(C_PDF), 0)
letter(XA, yC, "C")

# --- Row 2 right: B2 pose (centered) ---
yB2 = yC + (hC - hB) / 2 - 4
img = page.insert_image(pymupdf.Rect(XR, yB2, XR + RIGHT_W, yB2 + hB),
                        filename=flat["B2"])
cap_y2 = yB2 + hB + 8
caption(XR + RIGHT_W / 2, cap_y2, "GRB10", "-EGCG (-9.61 kcal/mol)")

# --- Bottom: D PLIP diagram centered ---
yD = TOP + LEFT_COL_H + 14
xD = (PAGE_W - D_W) / 2
page.show_pdf_page(pymupdf.Rect(xD, yD, xD + D_W, yD + hD),
                   pymupdf.open(D_PDF), 0)
letter(xD, yD, "D")

out_pdf = f"{FIN}/Figure9_composite.pdf"
doc.save(out_pdf, deflate=True)
doc.close()
print(f"saved: {out_pdf} ({PAGE_W/72:.2f} x {PAGE_H/72:.2f} in)")
print(f"  A {LEFT_W:.0f}x{hA:.0f}pt | B {RIGHT_W:.0f}x{hB:.0f}pt x2 | "
      f"C {LEFT_W:.0f}x{hC:.0f}pt | D {D_W}x{hD:.0f}pt")

# --- raster outputs ---
d = pymupdf.open(out_pdf)
pg = d[0]
pix = pg.get_pixmap(matrix=pymupdf.Matrix(300 / 72, 300 / 72), alpha=False)
Image.frombytes("RGB", (pix.width, pix.height), pix.samples).save(
    f"{FIN}/Figure9_composite.png", dpi=(300, 300))
pix = pg.get_pixmap(matrix=pymupdf.Matrix(1000 / 72, 1000 / 72), alpha=False)
Image.frombytes("RGB", (pix.width, pix.height), pix.samples).save(
    f"{FIN}/Figure9_composite.tiff", format="TIFF",
    compression="tiff_lzw", dpi=(1000, 1000))
d.close()
for ext in ["pdf", "png", "tiff"]:
    print(ext, round(os.path.getsize(f"{FIN}/Figure9_composite.{ext}") / 1e6, 2), "MB")
