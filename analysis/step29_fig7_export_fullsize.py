#!/usr/bin/env python
"""Step29 (v2): Figure 7 (GSVA top20) 全尺寸导出 — 不做 120mm 缩放。

输入（step28 生成的全尺寸矢量 PDF，Helvetica 文字，面板文字 8pt）：
  07_hub_immune/gsva_top20/vector/A_GSVA_top20_ARDS_full.pdf   (792 x 648 pt = 11x9in)
  07_hub_immune/gsva_top20/vector/B_GSVA_top20_Sepsis_full.pdf  (792 x 648 pt)

输出（覆盖 AJRCCM_submission/07_figure_exports/Figure 7/ 下 9 个文件）：
  A_GSVA_top20_ARDS.{pdf,png,tiff}      11 x 9 in 原尺寸
  B_GSVA_top20_Sepsis.{pdf,png,tiff}    11 x 9 in 原尺寸
  Figure5_composite.{pdf,png,tiff}       A/B 并排 + 字母标签（17 in 宽）

规范：PDF 为矢量（pypdf 变换不重蒸馏，文字保持 Helvetica 可编辑，8pt）；
      PNG/TIFF 1000 dpi 高清（满足彩色>=300 / 黑白>=500 / line art>=1000），
      RGB 白底，TIFF LZW 压缩。
"""
import os
import shutil
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf import PageObject
from PIL import Image

Image.MAX_IMAGE_PIXELS = None   # 23in x 1000dpi 合成图约 2.3 亿像素，解除 PIL 解压限制

ROOT = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
SRC = os.path.join(ROOT, "07_hub_immune/gsva_top20/vector")
FIG7 = os.path.join(ROOT, "AJRCCM_submission/07_figure_exports/Figure 7")

DPI_RASTER = 1000

for fp in ("/System/Library/Fonts/Helvetica.ttc",):
    if os.path.exists(fp):
        font_manager.fontManager.addfont(fp)
plt.rcParams["font.family"] = ["Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42

# ---------------------------------------------------------------- 1) A/B 独立面板：原尺寸直接复制
a_src = os.path.join(SRC, "A_GSVA_top20_ARDS_full.pdf")
b_src = os.path.join(SRC, "B_GSVA_top20_Sepsis_full.pdf")
shutil.copy(a_src, os.path.join(FIG7, "A_GSVA_top20_ARDS.pdf"))
shutil.copy(b_src, os.path.join(FIG7, "B_GSVA_top20_Sepsis.pdf"))
print("A/B panel PDFs copied at full size (11 x 9 in, 8pt text)")

# ---------------------------------------------------------------- 2) 组合图（A/B 并排 + 字母，全尺寸）
pa = PdfReader(a_src).pages[0]
pb = PdfReader(b_src).pages[0]
SW, SH = float(pa.mediabox.width), float(pa.mediabox.height)   # 792 x 648 pt

MARGIN = 24.0       # pt
GAP = 24.0          # pt
LETTER_ROW = 26.0   # 字母行高（含间距）
panel_w = SW                                    # 原尺寸 792 pt (11 in)
page_w = MARGIN * 2 + panel_w * 2 + GAP         # 1656 pt (23 in)
page_h = MARGIN + LETTER_ROW + SH + MARGIN      # 24+26+648+24 = 722 pt (~10 in)

comp = PageObject.create_blank_page(width=page_w, height=page_h)
y0 = MARGIN
comp.merge_transformed_page(pa, Transformation().scale(1.0).translate(MARGIN, y0))
comp.merge_transformed_page(pb, Transformation().scale(1.0).translate(MARGIN + panel_w + GAP, y0))

# 字母标签覆盖层（Helvetica Bold 18pt，位于各面板左上角上方）
ov = "/tmp/fig7_letters_overlay.pdf"
fig = plt.figure(figsize=(page_w / 72, page_h / 72))
for x_pt, ch in ((MARGIN, "A"), (MARGIN + panel_w + GAP, "B")):
    fig.text(x_pt / page_w, (page_h - MARGIN - 13.0) / page_h, ch, fontsize=18,
             fontweight="bold", ha="left", va="baseline", color="black")
fig.savefig(ov, transparent=True)
plt.close(fig)
lp = PdfReader(ov).pages[0]
comp.merge_page(lp)

w = PdfWriter()
w.add_page(comp)
comp_pdf = os.path.join(FIG7, "Figure5_composite.pdf")
with open(comp_pdf, "wb") as f:
    w.write(f)
print(f"composite: {page_w/72:.1f} x {page_h/72:.1f} in -> {comp_pdf}")

# ---------------------------------------------------------------- 3) 光栅导出（gs 1000dpi + PIL 规范化）
for stem in ("A_GSVA_top20_ARDS", "B_GSVA_top20_Sepsis", "Figure5_composite"):
    pdf = os.path.join(FIG7, stem + ".pdf")
    png_tmp = f"/tmp/{stem}_gs.png"
    tif_tmp = f"/tmp/{stem}_gs.tiff"
    subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=png16m",
                    f"-r{DPI_RASTER}", f"-sOutputFile={png_tmp}", pdf], check=True)
    subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=tiff24nc",
                    f"-r{DPI_RASTER}", f"-sOutputFile={tif_tmp}", pdf], check=True)
    for tmp, final, is_tiff in ((png_tmp, os.path.join(FIG7, stem + ".png"), False),
                                (tif_tmp, os.path.join(FIG7, stem + ".tiff"), True)):
        im = Image.open(tmp)
        if im.mode != "RGB":
            bg = Image.new("RGB", im.size, (255, 255, 255))
            if im.mode == "RGBA":
                bg.paste(im, mask=im.split()[3])
            else:
                bg.paste(im.convert("RGB"))
            im = bg
        if is_tiff:
            im.save(final, dpi=(DPI_RASTER, DPI_RASTER), compression="tiff_lzw")
        else:
            im.save(final, dpi=(DPI_RASTER, DPI_RASTER))
        print(f"  raster: {final}  {im.size}  @{DPI_RASTER}dpi RGB")
print("done.")
