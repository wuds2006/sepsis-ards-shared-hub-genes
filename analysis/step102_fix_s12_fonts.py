#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 102 —— 09 包 FigureS12_CellChat_pathway_comparison 字体修复
问题: 09 包内 S12 为早期 PIL 栅格整图（单张 3170x3153 位图、无文本层），字母不可编辑。
修复: 回源 4 个矢量面板（rankNet x2 + diff_LR_bubble x2），按现行 2x2 版式矢量重拼；
     面板字母 A-D 置于各面板左上方，嵌入式 Helvetica-Bold 9pt 黑色可编辑文本；
     Ghostscript 重蒸馏展平 Form XObject；重出 PNG/TIFF@300dpi LZW。
"""
import os, shutil, subprocess
import pymupdf
from fontTools.ttLib import TTCollection
from PIL import Image

R   = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
O9  = f"{R}/AJRCCM_submission/09_review_deAI_2026-09-25/figures/supplementary"
BK  = f"{O9}/_backup_pre_fontfix"
TMP = "/tmp/step102_s12"
os.makedirs(BK, exist_ok=True)
os.makedirs(TMP, exist_ok=True)
NAME = "FigureS12_CellChat_pathway_comparison"

# ---- Helvetica-Bold TTF（与 step98 同一来源）----
BOLD_TTF = f"{TMP}/Helvetica-Bold.ttf"
if not os.path.exists(BOLD_TTF):
    coll = TTCollection("/System/Library/Fonts/Helvetica.ttc")
    for f in coll.fonts:
        if f["name"].getDebugName(2) == "Bold" and f["name"].getDebugName(1) == "Helvetica":
            f.save(BOLD_TTF); break
assert os.path.exists(BOLD_TTF)

CC1 = f"{R}/08_singlecell/cellchat"            # GSE151263 ARDS vs Sepsis
CC2 = f"{R}/08_singlecell/GSE167363/cellchat"  # GSE167363 HC vs Sepsis

PANELS = [
    (f"{CC1}/rankNet_pathways.pdf", "A"),
    (f"{CC2}/rankNet_pathways.pdf", "B"),
    (f"{CC1}/diff_LR_bubble.pdf",   "C"),
    (f"{CC2}/diff_LR_bubble.pdf",   "D"),
]
for p, _ in PANELS:
    assert os.path.exists(p), p

# 几何：沿用旧 3050px@300dpi 网格 -> 页宽 760.8pt，双列 360pt
PX = 72.0 / 300.0
PAGE_W = (3050 + 2 * 60) * PX      # 760.8
MARGIN = 60 * PX                   # 14.4
GAP    = 50 * PX                   # 12.0
CELL   = (3050 - 50) / 2 * PX      # 360.0

def pdf_size(path):
    d = pymupdf.open(path); r = d[0].rect; d.close()
    return r.width, r.height

rows = [PANELS[:2], PANELS[2:]]
row_hs = []
for row in rows:
    row_hs.append(max(CELL * pdf_size(p)[1] / pdf_size(p)[0] for p, _ in row))

PAGE_H = MARGIN + row_hs[0] + GAP + row_hs[1] + MARGIN

doc = pymupdf.open()
page = doc.new_page(width=PAGE_W, height=PAGE_H)
page.insert_font(fontname="HelvBold", fontfile=BOLD_TTF)

y = MARGIN
for row, rh in zip(rows, row_hs):
    x = MARGIN
    for path, ch in row:
        pw, ph = pdf_size(path)
        h = CELL * ph / pw
        yy = y + (rh - h) / 2
        src = pymupdf.open(path)
        page.show_pdf_page(pymupdf.Rect(x, yy, x + CELL, yy + h), src, 0)
        src.close()
        # 字母在面板左上方（不遮挡面板内容），Helvetica-Bold 9pt 黑色
        page.insert_text((x, yy - 2.5), ch, fontsize=9,
                         fontname="HelvBold", color=(0, 0, 0))
        x += CELL + GAP
    y += rh + GAP

raw = f"{TMP}/{NAME}_raw.pdf"
doc.save(raw, deflate=True)
doc.close()

final = f"{TMP}/{NAME}.pdf"
subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite",
                "-dCompatibilityLevel=1.7", "-dEmbedAllFonts=true",
                "-dSubsetFonts=true", "-dNOPLATFONTS",
                f"-sOutputFile={final}", raw], check=True)

# 备份并替换
for ext in ["pdf", "png", "tiff"]:
    old = f"{O9}/{NAME}.{ext}"
    bkd = f"{BK}/{NAME}.{ext}"
    if os.path.exists(old) and not os.path.exists(bkd):
        shutil.copy(old, bkd)
shutil.copy(final, f"{O9}/{NAME}.pdf")

# 重出位图 @300dpi
d = pymupdf.open(final)
pix = d[0].get_pixmap(matrix=pymupdf.Matrix(300/72, 300/72), alpha=False)
im = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
im.save(f"{O9}/{NAME}.png", dpi=(300, 300))
im.save(f"{O9}/{NAME}.tiff", format="TIFF", compression="tiff_lzw", dpi=(300, 300))
d.close()

# 验证
d = pymupdf.open(f"{O9}/{NAME}.pdf"); pg = d[0]
letters = []
for b in pg.get_text("dict")["blocks"]:
    if b["type"] != 0: continue
    for l in b["lines"]:
        for s in l["spans"]:
            if len(s["text"]) == 1 and s["text"] in "ABCD" and round(s["size"], 1) == 9.0:
                letters.append((s["text"], s["font"], round(s["size"], 1), s["color"]))
print("page:", pg.rect)
print("fonts:", sorted(set(f[3] for f in pg.get_fonts())))
print("letters:", sorted(letters))
d.close()
print("DONE ->", O9)
