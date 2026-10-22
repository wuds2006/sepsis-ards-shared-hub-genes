#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 99 —— step98 复查修复
1) FigureS8: figS3_parts C/D/E 面板内嵌字母前缀 "C "/"D "/"E " 与拼装字母重复
   -> 红action 清除前缀后重建。
2) FigureS18: B/E 面板标题紧贴左上角，字母白衬底遮挡标题首字符
   -> 字母改为面板左侧间隙/页边放置 (x-8.5, 9pt Helvetica-Bold, 黑色, 无衬底)。
"""
import os, shutil, subprocess
import pymupdf

R   = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
O9  = f"{R}/AJRCCM_submission/09_review_deAI_2026-09-25/figures/supplementary"
TMP = "/tmp/step99_fontfix"
os.makedirs(TMP, exist_ok=True)
BOLD_TTF = "/tmp/step98_fontfix/Helvetica-Bold.ttf"
assert os.path.exists(BOLD_TTF)

MR = f"{R}/04_MR_SMR"
SC = f"{R}/08_singlecell"

def pdf_size(path):
    d = pymupdf.open(path); r = d[0].rect; d.close()
    return r.width, r.height

def strip_prefix(src, dst, letter, xmax, ymax):
    """清除面板左上角内嵌单字母前缀"""
    d = pymupdf.open(src); pg = d[0]
    for w in pg.get_text("words"):
        text, x0, y0, x1, y1 = w[4], w[0], w[1], w[2], w[3]
        if text == letter and x1 < xmax and y1 < ymax:
            pg.add_redact_annot(pymupdf.Rect(x0 - 1, y0 - 1, x1 + 2, y1 + 1),
                                fill=(1, 1, 1))
    pg.apply_redactions()
    d.save(dst, deflate=True); d.close()
    return dst

def gs_flatten(raw, final):
    subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite",
                    "-dCompatibilityLevel=1.7", "-dEmbedAllFonts=true",
                    "-dSubsetFonts=true", "-dNOPLATFONTS",
                    f"-sOutputFile={final}", raw], check=True)

def export_rasters(pdf_path, stem):
    from PIL import Image
    d = pymupdf.open(pdf_path); pg = d[0]
    pix = pg.get_pixmap(matrix=pymupdf.Matrix(300/72, 300/72), alpha=False)
    Image.frombytes("RGB", (pix.width, pix.height), pix.samples).save(
        f"{stem}.png", dpi=(300, 300))
    Image.frombytes("RGB", (pix.width, pix.height), pix.samples).save(
        f"{stem}.tiff", format="TIFF", compression="tiff_lzw", dpi=(300, 300))
    d.close()

# ---------------- Figure S8 ----------------
pC = strip_prefix(f"{MR}/figS3_parts/FigureS3C_sensitivity_forest.pdf",
                  f"{TMP}/S8C_clean.pdf", "C", 140, 30)
pD = strip_prefix(f"{MR}/figS3_parts/FigureS3D_SMR_HEIDI.pdf",
                  f"{TMP}/S8D_clean.pdf", "D", 90, 30)
pE = strip_prefix(f"{MR}/figS3_parts/FigureS3E_LDSC.pdf",
                  f"{TMP}/S8E_clean.pdf", "E", 70, 35)

PX5 = 72.0 / 500.0
rows = [
    [(f"{MR}/MR_sepsis4980_scatter.pdf", "A", 1550 * PX5),
     (f"{MR}/MR_sepsis4980_forest.pdf", "B", 1550 * PX5)],
    [(pC, "C", 1010 * PX5), (pD, "D", 1010 * PX5), (pE, "E", 1010 * PX5)],
]
S8_W = 3470 * PX5
gap = 30 * PX5
row_hs = []
for row in rows:
    row_hs.append(max(cw * pdf_size(p)[1] / pdf_size(p)[0] for p, _c, cw in row))
y1 = (150 + 170) * PX5
H = y1 + row_hs[0] + 200 * PX5 + row_hs[1] + 170 * PX5
doc = pymupdf.open()
page = doc.new_page(width=S8_W, height=H)
page.insert_font(fontname="HelvBold", fontfile=BOLD_TTF)
def panel_s8(path, ch, cw, x, yy):
    pw, ph = pdf_size(path); h = cw * ph / pw
    src = pymupdf.open(path)
    page.show_pdf_page(pymupdf.Rect(x, yy, x + cw, yy + h), src, 0)
    src.close()
    page.draw_rect(pymupdf.Rect(x, yy, x + 11, yy + 11), fill=(1, 1, 1), color=None)
    page.insert_text((x + 2.2, yy + 8.6), ch, fontsize=9,
                     fontname="HelvBold", color=(0, 0, 0))
x = 170 * PX5
for path, ch, cw in rows[0]:
    panel_s8(path, ch, cw, x, y1 + (row_hs[0] - cw * pdf_size(path)[1] / pdf_size(path)[0]) / 2)
    x += cw + gap
y2 = y1 + row_hs[0] + 200 * PX5
w2 = sum(cw for _p, _c, cw in rows[1]) + gap * 2
x = (S8_W - w2) / 2
for path, ch, cw in rows[1]:
    panel_s8(path, ch, cw, x, y2 + (row_hs[1] - cw * pdf_size(path)[1] / pdf_size(path)[0]) / 2)
    x += cw + gap
raw = f"{TMP}/FigureS8_MR_validation_raw.pdf"
doc.save(raw, deflate=True); doc.close()
gs_flatten(raw, f"{TMP}/FigureS8_MR_validation.pdf")

# ---------------- Figure S18 (字母放面板左侧间隙) ----------------
PX = 72.0 / 300.0
MARGIN = 60 * PX; GAP = 50 * PX
CELL3 = (3050 - 100) / 3 * PX
PAGE_W = 3170 * PX
rows18 = [
    [(f"{SC}/knockout_stratified/GSE151263_KO_stratified_GO_BP_UPP1.pdf", "A", CELL3),
     (f"{SC}/knockout_stratified/GSE151263_KO_stratified_top10_UPP1.pdf", "B", CELL3),
     (f"{SC}/celloracle_stratified/GSE151263_CellOracle_stratified_hubdelta.pdf", "C", CELL3)],
    [(f"{SC}/knockout_stratified/GSE167363_KO_stratified_GO_BP_CDK5RAP2.pdf", "D", CELL3),
     (f"{SC}/knockout_stratified/GSE167363_KO_stratified_top10_CDK5RAP2.pdf", "E", CELL3),
     (f"{SC}/celloracle_stratified/GSE167363_CellOracle_stratified_hubdelta.pdf", "F", CELL3)],
]
row_hs = [max(cw * pdf_size(p)[1] / pdf_size(p)[0] for p, _c, cw in row) for row in rows18]
H = MARGIN + sum(row_hs) + GAP + MARGIN
doc = pymupdf.open()
page = doc.new_page(width=PAGE_W, height=H)
page.insert_font(fontname="HelvBold", fontfile=BOLD_TTF)
y = MARGIN
for row, rh in zip(rows18, row_hs):
    row_w = sum(cw for _p, _c, cw in row) + GAP * (len(row) - 1)
    x = (PAGE_W - row_w) / 2
    for path, ch, cw in row:
        pw, ph = pdf_size(path); h = cw * ph / pw
        yy = y + (rh - h) / 2
        src = pymupdf.open(path)
        page.show_pdf_page(pymupdf.Rect(x, yy, x + cw, yy + h), src, 0)
        src.close()
        # 字母: 面板左侧间隙, 9pt Helvetica-Bold 黑, 无衬底
        page.insert_text((x - 8.5, yy + 8.6), ch, fontsize=9,
                         fontname="HelvBold", color=(0, 0, 0))
        x += cw + GAP
    y += rh + GAP
raw = f"{TMP}/FigureS18_stratified_knockout_raw.pdf"
doc.save(raw, deflate=True); doc.close()
gs_flatten(raw, f"{TMP}/FigureS18_stratified_knockout.pdf")

# ---------------- 替换 + 重出位图 ----------------
for name in ["FigureS8_MR_validation", "FigureS18_stratified_knockout"]:
    shutil.copy(f"{TMP}/{name}.pdf", f"{O9}/{name}.pdf")
    export_rasters(f"{TMP}/{name}.pdf", f"{O9}/{name}")
    print("fixed", name)
print("DONE")
