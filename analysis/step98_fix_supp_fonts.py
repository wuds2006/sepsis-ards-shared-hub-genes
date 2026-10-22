#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 98 —— 09_review_deAI_2026-09-25 补充图字体修复（7 张）
问题: S1/S8/S11/S13/S14/S15/S18 为 PIL 栅格拼装，无文本层，字母不可编辑。
修复: 回源矢量面板 PDF，按原版式矢量重拼；字母改为嵌入式 Helvetica-Bold 9pt 黑色
      可编辑文本（白底小衬底，沿用 step25 原样式）；Ghostscript 重蒸馏展平
      Form XObject 保证 Illustrator 可编辑；重新导出 PNG@300 + TIFF@300 LZW。
"""
import os, shutil, subprocess
import pymupdf
from fontTools.ttLib import TTCollection

R   = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
SUB = f"{R}/AJRCCM_submission"
O9  = f"{SUB}/09_review_deAI_2026-09-25/figures/supplementary"
BK  = f"{O9}/_backup_pre_fontfix"
TMP = "/tmp/step98_fontfix"
os.makedirs(BK, exist_ok=True)
os.makedirs(TMP, exist_ok=True)

# ---- 1. 从 Helvetica.ttc 提取 Helvetica-Bold (index 1) 为独立 TTF 以便嵌入 ----
BOLD_TTF = f"{TMP}/Helvetica-Bold.ttf"
if not os.path.exists(BOLD_TTF):
    coll = TTCollection("/System/Library/Fonts/Helvetica.ttc")
    for f in coll.fonts:
        if f["name"].getDebugName(2) == "Bold" and f["name"].getDebugName(1) == "Helvetica":
            f.save(BOLD_TTF); break
assert os.path.exists(BOLD_TTF), "Helvetica-Bold extraction failed"

DPI_OLD = 300.0          # step25/26 栅格底稿 dpi: 1px = 72/300 pt
PX = 72.0 / DPI_OLD
MARGIN = 60 * PX         # 14.4
GAP    = 50 * PX         # 12.0
FULL   = 3050
CELL3  = (FULL - 2 * 50) / 3 * PX     # 236.0
CELL4  = (FULL - 3 * 50) / 4 * PX     # 174.0
PAGE_W = (FULL + 2 * 60) * PX         # 760.8

WG  = f"{R}/02_WGCNA/figS1_parts"
MR  = f"{R}/04_MR_SMR"
G167= f"{R}/08_singlecell/GSE167363"
SC  = f"{R}/08_singlecell"

def pdf_size(path):
    d = pymupdf.open(path); r = d[0].rect; d.close()
    return r.width, r.height

def build(name, rows, page_w=PAGE_W, margin=MARGIN, gap=GAP,
          row_pad_top=None, row_pad_bot=None):
    """rows: [[(pdf_path, letter, cell_w_pt), ...], ...]  行高=该行最高面板; 行内垂直居中; 行水平居中"""
    row_hs = []
    for row in rows:
        hs = []
        for path, _ch, cw in row:
            pw, ph = pdf_size(path)
            hs.append(cw * ph / pw)
        row_hs.append(max(hs))
    if row_pad_top is None: row_pad_top = [0.0] * len(rows)
    if row_pad_bot is None: row_pad_bot = [0.0] * len(rows)
    H = margin + sum(row_hs) + gap * (len(rows) - 1) + sum(row_pad_top) + sum(row_pad_bot) + margin
    doc = pymupdf.open()
    page = doc.new_page(width=page_w, height=H)
    page.insert_font(fontname="HelvBold", fontfile=BOLD_TTF)
    y = margin
    for row, rh, pt, pb in zip(rows, row_hs, row_pad_top, row_pad_bot):
        y += pt
        row_w = sum(cw for _p, _c, cw in row) + gap * (len(row) - 1)
        x = (page_w - row_w) / 2
        for path, ch, cw in row:
            pw, ph = pdf_size(path)
            h = cw * ph / pw
            yy = y + (rh - h) / 2
            src = pymupdf.open(path)
            page.show_pdf_page(pymupdf.Rect(x, yy, x + cw, yy + h), src, 0)
            src.close()
            # 字母: 白底衬底 + Helvetica-Bold 9pt 黑色 (cell 左上角内)
            page.draw_rect(pymupdf.Rect(x, yy, x + 11, yy + 11), fill=(1, 1, 1), color=None)
            page.insert_text((x + 2.2, yy + 8.6), ch, fontsize=9,
                             fontname="HelvBold", color=(0, 0, 0))
            x += cw + gap
        y += rh + pb + gap
    out_raw = f"{TMP}/{name}_raw.pdf"
    doc.save(out_raw, deflate=True)
    doc.close()
    # Ghostscript 重蒸馏: 展平 Form XObject + 嵌入字体 -> Illustrator 可编辑
    out_final = f"{TMP}/{name}.pdf"
    subprocess.run([
        "gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite",
        "-dCompatibilityLevel=1.7", "-dEmbedAllFonts=true", "-dSubsetFonts=true",
        "-dNOPLATFONTS", f"-sOutputFile={out_final}", out_raw], check=True)
    return out_final

def export_rasters(pdf_path, stem):
    from PIL import Image
    d = pymupdf.open(pdf_path); pg = d[0]
    pix = pg.get_pixmap(matrix=pymupdf.Matrix(300/72, 300/72), alpha=False)
    Image.frombytes("RGB", (pix.width, pix.height), pix.samples).save(
        f"{stem}.png", dpi=(300, 300))
    Image.frombytes("RGB", (pix.width, pix.height), pix.samples).save(
        f"{stem}.tiff", format="TIFF", compression="tiff_lzw", dpi=(300, 300))
    d.close()

SPECS = {
    "FigureS1_WGCNA_soft_threshold_dendrograms": [
        [(f"{WG}/a_sft_GSE76293.pdf", "A", CELL3),
         (f"{WG}/b_conn_GSE76293.pdf", "B", CELL3),
         (f"{WG}/c_dendro_GSE76293.pdf", "C", CELL3)],
        [(f"{WG}/d_sft_GSE69528.pdf", "D", CELL3),
         (f"{WG}/e_conn_GSE69528.pdf", "E", CELL3),
         (f"{WG}/f_dendro_GSE69528.pdf", "F", CELL3)],
    ],
    "FigureS11_single_cell_atlas_GSE167363": [
        [(f"{G167}/GSE167363_UMAP.pdf", "A", CELL3),
         (f"{G167}/GSE167363_celltype_prop.pdf", "B", CELL3),
         (f"{G167}/GSE167363_marker_dotplot.pdf", "C", CELL3)],
        [(f"{G167}/GSE167363_hub_dotplot.pdf", "D", CELL3),
         (f"{G167}/GSE167363_hub_featureplot.pdf", "E", CELL3),
         (f"{G167}/GSE167363_hub_violin_group.pdf", "F", CELL3)],
    ],
    "FigureS13_pseudotime_AUCell_CellOracle": [
        [(f"{G167}/slingshot/GSE167363_slingshot_overview.pdf", "A", CELL3),
         (f"{G167}/slingshot/GSE167363_slingshot_trajectory_PCA.pdf", "B", CELL3),
         (f"{G167}/slingshot/GSE167363_slingshot_gene_trend.pdf", "C", CELL3)],
        [(f"{G167}/aucell/GSE167363_AUCell_violin.pdf", "D", CELL3),
         (f"{G167}/aucell/GSE167363_AUCell_heatmap.pdf", "E", CELL3),
         (f"{G167}/aucell/GSE167363_AUCell_umap.pdf", "F", CELL3)],
        [(f"{G167}/celloracle/CellOracle_KO_effect_UPP1.pdf", "G", CELL4),
         (f"{G167}/celloracle/CellOracle_KO_effect_CDK5RAP2.pdf", "H", CELL4),
         (f"{G167}/celloracle/CellOracle_KO_effect_EXOSC4.pdf", "I", CELL4),
         (f"{G167}/celloracle/CellOracle_KO_effect_VNN1.pdf", "J", CELL4)],
    ],
    "FigureS14_monocle3_trajectories": [
        [(f"{SC}/monocle3/GSE151263_monocle3_trajectory_celltype.pdf", "A", CELL3),
         (f"{SC}/monocle3/GSE151263_monocle3_pseudotime_group.pdf", "B", CELL3),
         (f"{SC}/monocle3/GSE151263_monocle3_hub_moransI_bar.pdf", "C", CELL3)],
        [(f"{SC}/monocle3/GSE167363_monocle3_trajectory_celltype.pdf", "D", CELL3),
         (f"{SC}/monocle3/GSE167363_monocle3_pseudotime_group.pdf", "E", CELL3),
         (f"{SC}/monocle3/GSE167363_monocle3_hub_moransI_bar.pdf", "F", CELL3)],
    ],
    "FigureS15_metabolic_pathways": [
        [(f"{SC}/metabolism/GSE151263_metab_pathway_heatmap.pdf", "A", CELL3),
         (f"{SC}/metabolism/GSE151263_metab_group_difference.pdf", "B", CELL3),
         (f"{SC}/metabolism/GSE151263_metab_hub_correlation.pdf", "C", CELL3)],
        [(f"{SC}/metabolism/GSE167363_metab_pathway_heatmap.pdf", "D", CELL3),
         (f"{SC}/metabolism/GSE167363_metab_group_difference.pdf", "E", CELL3),
         (f"{SC}/metabolism/GSE167363_metab_hub_correlation.pdf", "F", CELL3)],
    ],
    "FigureS18_stratified_knockout": [
        [(f"{SC}/knockout_stratified/GSE151263_KO_stratified_GO_BP_UPP1.pdf", "A", CELL3),
         (f"{SC}/knockout_stratified/GSE151263_KO_stratified_top10_UPP1.pdf", "B", CELL3),
         (f"{SC}/celloracle_stratified/GSE151263_CellOracle_stratified_hubdelta.pdf", "C", CELL3)],
        [(f"{SC}/knockout_stratified/GSE167363_KO_stratified_GO_BP_CDK5RAP2.pdf", "D", CELL3),
         (f"{SC}/knockout_stratified/GSE167363_KO_stratified_top10_CDK5RAP2.pdf", "E", CELL3),
         (f"{SC}/celloracle_stratified/GSE167363_CellOracle_stratified_hubdelta.pdf", "F", CELL3)],
    ],
}

# ---- S8: MR validation (step21 版式, 500dpi 像素坐标 -> pt) ----
PX5 = 72.0 / 500.0
S8_ROWS = [
    [(f"{MR}/MR_sepsis4980_scatter.pdf", "A", 1550 * PX5),
     (f"{MR}/MR_sepsis4980_forest.pdf", "B", 1550 * PX5)],
    [(f"{MR}/figS3_parts/FigureS3C_sensitivity_forest.pdf", "C", 1010 * PX5),
     (f"{MR}/figS3_parts/FigureS3D_SMR_HEIDI.pdf", "D", 1010 * PX5),
     (f"{MR}/figS3_parts/FigureS3E_LDSC.pdf", "E", 1010 * PX5)],
]
S8_W = 3470 * PX5          # 499.68
S8_ROWS_PAD_TOP = [(150 + 170) * PX5 - 24.48, 200 * PX5]   # 首行 LBL+PAD, 次行 PAD 间隙
# 简化: row1 顶距 = (LBL+PAD)*PX5; margin 参数并入 pad
def build_s8():
    rows = S8_ROWS
    gap = 30 * PX5
    row_hs = []
    for row in rows:
        hs = []
        for path, _ch, cw in row:
            pw, ph = pdf_size(path)
            hs.append(cw * ph / pw)
        row_hs.append(max(hs))
    y1 = (150 + 170) * PX5                     # row1 顶
    H = y1 + row_hs[0] + 200 * PX5 + row_hs[1] + 170 * PX5
    doc = pymupdf.open()
    page = doc.new_page(width=S8_W, height=H)
    page.insert_font(fontname="HelvBold", fontfile=BOLD_TTF)
    # row1: A 左对齐 PAD, B 紧随其后 (与原图一致)
    x = 170 * PX5
    for path, ch, cw in rows[0]:
        pw, ph = pdf_size(path); h = cw * ph / pw
        yy = y1 + (row_hs[0] - h) / 2
        src = pymupdf.open(path)
        page.show_pdf_page(pymupdf.Rect(x, yy, x + cw, yy + h), src, 0)
        src.close()
        page.draw_rect(pymupdf.Rect(x, yy, x + 11, yy + 11), fill=(1, 1, 1), color=None)
        page.insert_text((x + 2.2, yy + 8.6), ch, fontsize=9,
                         fontname="HelvBold", color=(0, 0, 0))
        x += cw + gap
    # row2: C/D/E 水平居中
    y2 = y1 + row_hs[0] + 200 * PX5
    w2 = sum(cw for _p, _c, cw in rows[1]) + gap * 2
    x = (S8_W - w2) / 2
    for path, ch, cw in rows[1]:
        pw, ph = pdf_size(path); h = cw * ph / pw
        yy = y2 + (row_hs[1] - h) / 2
        src = pymupdf.open(path)
        page.show_pdf_page(pymupdf.Rect(x, yy, x + cw, yy + h), src, 0)
        src.close()
        page.draw_rect(pymupdf.Rect(x, yy, x + 11, yy + 11), fill=(1, 1, 1), color=None)
        page.insert_text((x + 2.2, yy + 8.6), ch, fontsize=9,
                         fontname="HelvBold", color=(0, 0, 0))
        x += cw + gap
    raw = f"{TMP}/FigureS8_MR_validation_raw.pdf"
    doc.save(raw, deflate=True); doc.close()
    final = f"{TMP}/FigureS8_MR_validation.pdf"
    subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite",
                    "-dCompatibilityLevel=1.7", "-dEmbedAllFonts=true",
                    "-dSubsetFonts=true", "-dNOPLATFONTS",
                    f"-sOutputFile={final}", raw], check=True)
    return final

# ---- 主流程 ----
built = {}
for name, rows in SPECS.items():
    built[name] = build(name, rows)
    print("built", name)
built["FigureS8_MR_validation"] = build_s8()
print("built FigureS8_MR_validation")

# 备份原件并替换
for name, newpdf in built.items():
    for ext in ["pdf", "png", "tiff"]:
        old = f"{O9}/{name}.{ext}"
        bkd = f"{BK}/{name}.{ext}"
        if os.path.exists(old) and not os.path.exists(bkd):
            shutil.copy(old, bkd)
    stem = f"{O9}/{name}"
    shutil.copy(newpdf, f"{stem}.pdf")
    export_rasters(newpdf, stem)
    print("replaced", name)

print("ALL DONE ->", O9)
