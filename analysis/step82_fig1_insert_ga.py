#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 82 —— 将 Figure1_composite 的面板 A（旧 pipeline 示意图）替换为
04_graphical_abstract/graphical_abstract.pdf（三区布局 graphical abstract）。
做法：GA 全宽置顶（保持纵横比），旧页面 B/C/D 区域（y>=118pt）整体下移，
页高相应增加；面板字母 A 按原样式（Helvetica-Bold 9pt）重标。
产出：Figures_final/Figure1_composite.{pdf,png(300dpi),tiff(LZW)}
TIFF 若超 20 MB 上限则降为 600 dpi（ATS 组合图合规下限）。
"""
import os
import pymupdf
from PIL import Image

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission"
FIN = f"{BASE}/03_figures/Figures_final"
GA = f"{BASE}/04_graphical_abstract/graphical_abstract.pdf"
OLD = f"{FIN}/Figure1_composite.pdf"

old = pymupdf.open(OLD)
opg = old[0]
PAGE_W = opg.rect.width           # 504.6 pt
OLD_H = opg.rect.height           # 792.7 pt
CLIP_TOP = 118.0                  # 面板 A 底界（B 字母 y=125 上方）

ga = pymupdf.open(GA)
gr = ga[0].rect                   # 1080 x 756 pt

X = 23.0                          # 与旧图内容左边距一致
GA_Y = 20.0
GA_W = PAGE_W - 2 * X
GA_H = GA_W * gr.height / gr.width
GAP = 12.0
OFFSET = GA_Y + GA_H + GAP - CLIP_TOP
NEW_H = OLD_H + OFFSET

doc = pymupdf.open()
page = doc.new_page(width=PAGE_W, height=NEW_H)
# 新面板 A：graphical abstract
page.show_pdf_page(pymupdf.Rect(X, GA_Y, X + GA_W, GA_Y + GA_H), ga, 0)
# 旧 B/C/D 区域整体下移
page.show_pdf_page(
    pymupdf.Rect(0, CLIP_TOP + OFFSET, PAGE_W, OLD_H + OFFSET),
    old, 0, clip=pymupdf.Rect(0, CLIP_TOP, PAGE_W, OLD_H))
# 新面板字母 A（样式与旧图一致：Helvetica-Bold 9pt）
page.insert_text((10, GA_Y + 7.3), "A", fontsize=9, fontname="hebo", color=(0, 0, 0))

out_pdf = f"{FIN}/Figure1_composite.pdf"
doc.save(out_pdf, deflate=True, garbage=3)
doc.close(); old.close(); ga.close()
print(f"saved: {out_pdf} ({PAGE_W/72:.2f} x {NEW_H/72:.2f} in), GA panel {GA_W:.0f}x{GA_H:.0f}pt")

# --- raster outputs ---
d = pymupdf.open(out_pdf)
pg = d[0]
pix = pg.get_pixmap(matrix=pymupdf.Matrix(300 / 72, 300 / 72), alpha=False)
Image.frombytes("RGB", (pix.width, pix.height), pix.samples).save(
    f"{FIN}/Figure1_composite.png", dpi=(300, 300))

dpi_tiff = 1000
pix = pg.get_pixmap(matrix=pymupdf.Matrix(dpi_tiff / 72, dpi_tiff / 72), alpha=False)
img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
img.save(f"{FIN}/Figure1_composite.tiff", format="TIFF",
         compression="tiff_lzw", dpi=(dpi_tiff, dpi_tiff))
sz = os.path.getsize(f"{FIN}/Figure1_composite.tiff") / 1e6
if sz > 20:
    dpi_tiff = 600
    pix = pg.get_pixmap(matrix=pymupdf.Matrix(dpi_tiff / 72, dpi_tiff / 72), alpha=False)
    img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    img.save(f"{FIN}/Figure1_composite.tiff", format="TIFF",
             compression="tiff_lzw", dpi=(dpi_tiff, dpi_tiff))
    sz = os.path.getsize(f"{FIN}/Figure1_composite.tiff") / 1e6
d.close()
for ext in ["pdf", "png", "tiff"]:
    print(ext, round(os.path.getsize(f"{FIN}/Figure1_composite.{ext}") / 1e6, 2), "MB")
print("TIFF dpi:", dpi_tiff, f"({sz:.2f} MB)")
