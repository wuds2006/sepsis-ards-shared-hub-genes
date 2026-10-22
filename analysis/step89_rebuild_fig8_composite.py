#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 89 —— 重建 Figure8 composite（免疫微环境）
修复:
1) step87 中 F 面板高度按 W3 估算、实际按 half 放置, 页面高度少算 ~72 pt,
   导致 F(4 个相关热图 2x2)底部两个小图被裁掉 -> 本版按实际放置宽度精确计算行高。
2) G 与 H(两个 GSVA top20 热图)合并到同一行并排展示。
新布局(字母 A-H 空间顺序连续):
  Row1: A (hub6 GSVA violin, 全宽)
  Row2: B | C (CIBERSORT discovery)
  Row3: D | E (CIBERSORT validation)
  Row4: F (spearman corrplot 2x2, 居中, 完整显示)
  Row5: G | H (GSVA top20, 并排)
输出: 03_figures_RESULTS_NPG/composites/Figure8_composite.{pdf,png@300,tiff@1000 LZW}
"""
import os
import pymupdf
from PIL import Image

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission"
P  = f"{BASE}/03_figures/Figures_ggplot2_NPG/main_panels"
OUT = f"{BASE}/03_figures_RESULTS_NPG/composites"

MM = 72 / 25.4
PAGE_W = 178 * MM          # 504.57 pt
MARGIN = 8 * MM
CW = PAGE_W - 2 * MARGIN   # 459.2 pt
GAP_X = 4 * MM
GAP_Y = 5 * MM
TOP = 13 * MM
BOT = 6 * MM

def pdf_size(path):
    d = pymupdf.open(path); r = d[0].rect; d.close()
    return r.width, r.height

def aspect_h(path, w):
    pw, ph = pdf_size(path)
    return w * ph / pw

class Sheet:
    def __init__(self, height_pt):
        self.doc = pymupdf.open()
        self.page = self.doc.new_page(width=PAGE_W, height=height_pt)
    def letter(self, x, y, ch):
        self.page.insert_text((x - 10, y + 3.4), ch, fontsize=9,
                              fontname="hebo", color=(0, 0, 0))
    def pdf_panel(self, path, x, y, w, ch):
        src = pymupdf.open(path)
        pw, ph = pdf_size(path)
        h = w * ph / pw
        self.page.show_pdf_page(pymupdf.Rect(x, y, x + w, y + h), src, 0)
        src.close()
        if ch: self.letter(x, y, ch)
        return h
    def save(self, name):
        out = os.path.join(OUT, f"{name}_composite.pdf")
        self.doc.save(out, deflate=True)
        self.doc.close()
        return out

def export_rasters(pdf_path, name):
    d = pymupdf.open(pdf_path)
    pg = d[0]
    pix = pg.get_pixmap(matrix=pymupdf.Matrix(300/72, 300/72), alpha=False)
    Image.frombytes("RGB", (pix.width, pix.height), pix.samples).save(
        os.path.join(OUT, f"{name}_composite.png"), dpi=(300, 300))
    pix = pg.get_pixmap(matrix=pymupdf.Matrix(1000/72, 1000/72), alpha=False)
    Image.frombytes("RGB", (pix.width, pix.height), pix.samples).save(
        os.path.join(OUT, f"{name}_composite.tiff"),
        format="TIFF", compression="tiff_lzw", dpi=(1000, 1000))
    d.close()

half = (CW - GAP_X) / 2

# 面板路径
fA = f"{P}/Fig5A_hub6_GSVA.pdf"
fB = f"{P}/Fig5B_CIBERSORT_stacked_GSE76293.pdf"
fC = f"{P}/Fig5C_CIBERSORT_stacked_GSE69528.pdf"
fD = f"{P}/Fig5D_CIBERSORT_stacked_GSE243066.pdf"
fE = f"{P}/Fig5E_CIBERSORT_stacked_GSE100159.pdf"
fF = f"{P}/Fig5F_spearman_corrplot.pdf"
fG = f"{P}/Fig6A_GSVA_top20_GSE76293.pdf"
fH = f"{P}/Fig6B_GSVA_top20_GSE69528.pdf"

# 精确行高(按实际放置宽度)
F_W = 280                            # F 居中, 280 pt 宽, 保证 4 个热图完整且可读
hA = aspect_h(fA, CW)                # 全宽小提琴
hB = aspect_h(fB, half)              # CIBERSORT
hF = aspect_h(fF, F_W)               # 相关热图(完整)
hG = aspect_h(fG, half)              # GSVA top20

H = (TOP + hA + GAP_Y + hB + GAP_Y + hB + GAP_Y
     + hF + GAP_Y + hG + BOT)
s = Sheet(H)

y = TOP
s.pdf_panel(fA, MARGIN, y, CW, "A")
y += hA + GAP_Y
s.pdf_panel(fB, MARGIN, y, half, "B")
s.pdf_panel(fC, MARGIN + half + GAP_X, y, half, "C")
y += hB + GAP_Y
s.pdf_panel(fD, MARGIN, y, half, "D")
s.pdf_panel(fE, MARGIN + half + GAP_X, y, half, "E")
y += hB + GAP_Y
s.pdf_panel(fF, (PAGE_W - F_W) / 2, y, F_W, "F")   # 居中, 完整 2x2
y += hF + GAP_Y
s.pdf_panel(fG, MARGIN, y, half, "G")              # G 与 H 同一行并排
s.pdf_panel(fH, MARGIN + half + GAP_X, y, half, "H")

out = s.save("Figure8")
export_rasters(out, "Figure8")
print("Figure8 rebuilt:", round(PAGE_W/72, 2), "x", round(H/72, 2), "in ->", out)
