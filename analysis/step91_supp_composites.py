#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 91 —— 05_supplementary_NPG 补充图拼装 (FigureS1-S14)
- 16 个新 R 面板 (step90) + 复用 NPG 面板 + 2 张整体拷贝 (S6/S12)
- 布局全 pt 单位; 字母由本脚本统一绘制
- FigS1A/FigS1B/FigS14 源面板内嵌 'A1'/'B1'/'A' 标题前缀 -> 拷贝到本地 panels/ 后红action 清除
输出: composites/FigureS*_composite.{pdf,png@300,tiff@1000 LZW}
"""
import os, shutil
import pymupdf
from PIL import Image

R   = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
SUB = f"{R}/AJRCCM_submission"
S   = f"{SUB}/03_figures/Figures_ggplot2_NPG/supp_panels"
W   = f"{R}/02_WGCNA/figS1_parts"
M   = f"{R}/04_MR_SMR"
C   = f"{R}/08_singlecell"
O   = f"{SUB}/05_supplementary/Supplementary_Figures_Tables/Figures"
NP  = f"{SUB}/05_supplementary_NPG/panels"
OUT = f"{SUB}/05_supplementary_NPG/composites"
os.makedirs(NP, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

MM = 72 / 25.4
PAGE_W = 178 * MM
MARGIN = 8 * MM
CW = PAGE_W - 2 * MARGIN
GAP_X = 4 * MM
GAP_Y = 5 * MM
TOP = 13 * MM
BOT = 6 * MM
half = (CW - GAP_X) / 2

def pdf_size(path):
    d = pymupdf.open(path); r = d[0].rect; d.close()
    return r.width, r.height

def aspect_h(path, w):
    pw, ph = pdf_size(path)
    return w * ph / pw

def strip_words(src, dst, words_xmax_ymax):
    """拷贝面板并红action 清除左上角内嵌字母前缀。
    words_xmax_ymax: list of (word, xmax, ymax) 匹配规则"""
    d = pymupdf.open(src)
    pg = d[0]
    for w in pg.get_text("words"):
        text, x0, y0, x1, y1 = w[4], w[0], w[1], w[2], w[3]
        for pat, xm, ym in words_xmax_ymax:
            if text == pat and x1 < xm and y1 < ym:
                pg.add_redact_annot(pymupdf.Rect(x0 - 1, y0 - 1, x1 + 2.5, y1 + 1),
                                    fill=(1, 1, 1))
    pg.apply_redactions()
    d.save(dst, deflate=True)
    d.close()
    return dst

# 预处理: 清除内嵌前缀
strip_words(f"{S}/FigS1A_soft_threshold_GSE76293.pdf",
            f"{NP}/FigS1A_soft_threshold_GSE76293.pdf", [("A1", 90, 30)])
strip_words(f"{S}/FigS1B_soft_threshold_GSE69528.pdf",
            f"{NP}/FigS1B_soft_threshold_GSE69528.pdf", [("B1", 90, 30)])
strip_words(f"{S}/FigS14_donor_sensitivity.pdf",
            f"{NP}/FigS14_donor_sensitivity.pdf", [("A", 90, 30)])

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

def build(name, rows):
    """rows: list of row-specs; 每个 row = [(path, letter, width_pt or 'half'), ...]
    单面板 row 自动居中; 双面板 row 按 half 排列。行高取该行最高面板。"""
    row_hs = []
    for row in rows:
        hs = []
        for path, _ch, wspec in row:
            w = half if wspec == "half" else wspec
            hs.append(aspect_h(path, w))
        row_hs.append(max(hs))
    H = TOP + sum(row_hs) + GAP_Y * (len(rows) - 1) + BOT
    s = Sheet(H)
    y = TOP
    for row, rh in zip(rows, row_hs):
        if len(row) == 1:
            path, ch, wspec = row[0]
            w = half if wspec == "half" else wspec
            s.pdf_panel(path, (PAGE_W - w) / 2, y, w, ch)
        else:
            x = MARGIN
            for path, ch, wspec in row:
                w = half if wspec == "half" else wspec
                s.pdf_panel(path, x, y, w, ch)
                x += w + GAP_X
        y += rh + GAP_Y
    out = s.save(name)
    export_rasters(out, name)
    print(name, round(H/72, 2), "in")

# ---------- Figure S1: WGCNA network diagnostics ----------
build("FigureS1", [
    [(f"{NP}/FigS1A_soft_threshold_GSE76293.pdf", "A", "half"),
     (f"{NP}/FigS1B_soft_threshold_GSE69528.pdf", "B", "half")],
    [(f"{W}/b_conn_GSE76293.pdf", "C", "half"),
     (f"{W}/e_conn_GSE69528.pdf", "D", "half")],
    [(f"{W}/c_dendro_GSE76293.pdf", "E", "half"),
     (f"{W}/f_dendro_GSE69528.pdf", "F", "half")],
])

# ---------- Figure S2: hub-gene co-expression evidence ----------
build("FigureS2", [
    [(f"{NP}/RS2A_hub_kME.pdf", "A", "half"),
     (f"{NP}/RS2B_hub_degree_179CG.pdf", "B", "half")],
])

# ---------- Figure S3: enrichment supplements ----------
build("FigureS3", [
    [(f"{NP}/RS3A_KEGG_179_bar.pdf", "A", "half"),
     (f"{NP}/RS3B_Reactome_hub_dot.pdf", "B", "half")],
    [(f"{NP}/RS3C_KEGG2021_hub_dot.pdf", "C", "half"),
     (f"{NP}/RS3D_GOBP2023_hub_dot.pdf", "D", "half")],
])

# ---------- Figure S4: ML OOF ROC + confusion ----------
build("FigureS4", [
    [(f"{S}/FigS2A_OOF_ROC_GSE76293.pdf", "A", "half"),
     (f"{S}/FigS2B_OOF_ROC_GSE69528.pdf", "B", "half")],
    [(f"{S}/FigS2C_confusion_GSE76293.pdf", "C", "half"),
     (f"{S}/FigS2D_confusion_GSE69528.pdf", "D", "half")],
])

# ---------- Figure S5: ML robustness ----------
build("FigureS5", [
    [(f"{NP}/RS5A_apparent_vs_cv_GSE76293.pdf", "A", "half"),
     (f"{NP}/RS5B_apparent_vs_cv_GSE69528.pdf", "B", "half")],
    [(f"{NP}/RS5C_lasso_bootstrap_ARDS.pdf", "C", "half"),
     (f"{NP}/RS5D_lasso_bootstrap_Sepsis.pdf", "D", "half")],
    [(f"{NP}/RS5E_jaccard_ARDS.pdf", "E", "half"),
     (f"{NP}/RS5F_jaccard_Sepsis.pdf", "F", "half")],
    [(f"{NP}/RS5G_permutation_null.pdf", "G", "half")],
])

# ---------- Figure S6: single-gene GSEA (整体拷贝) ----------
for ext in ["pdf", "png", "tiff"]:
    shutil.copy(f"{O}/FigureS4_single_gene_gsea.{ext}",
                f"{OUT}/FigureS6_composite.{ext}")
print("FigureS6 copied from FigureS4_single_gene_gsea")

# ---------- Figure S7: external validation diagnostics ----------
build("FigureS7", [
    [(f"{S}/FigS16_cohort_PCA.pdf", "A", "half"),
     (f"{NP}/RS7B_MUC1_sensitivity.pdf", "B", "half")],
    [(f"{NP}/RS7C_AUC1_exact_lower_bounds.pdf", "C", "half"),
     (f"{NP}/RS7D_sepsis_candidates_in_ARDS.pdf", "D", "half")],
])

# ---------- Figure S8: MR sensitivity ----------
build("FigureS8", [
    [(f"{M}/MR_scatter.pdf", "A", "half"),
     (f"{M}/MR_sepsis4980_scatter.pdf", "B", "half")],
    [(f"{S}/FigS3_MR_method_forest.pdf", "C", "half"),
     (f"{M}/figS3_parts/FigureS3C_sensitivity_forest.pdf", "D", "half")],
    [(f"{M}/figS3_parts/FigureS3D_SMR_HEIDI.pdf", "E", "half"),
     (f"{M}/figS3_parts/FigureS3E_LDSC.pdf", "F", "half")],
])

# ---------- Figure S9: ssGSEA immune infiltration ----------
build("FigureS9", [
    [(f"{S}/FigS15_box_GSE76293.pdf", "A", "half"),
     (f"{S}/FigS15_box_GSE69528.pdf", "B", "half")],
    [(f"{S}/FigS15_box_GSE243066.pdf", "C", "half"),
     (f"{S}/FigS15_box_GSE100159.pdf", "D", "half")],
    [(f"{S}/FigS15_cor_GSE76293.pdf", "E", "half"),
     (f"{S}/FigS15_cor_GSE69528.pdf", "F", "half")],
    [(f"{S}/FigS15_cor_GSE243066.pdf", "G", "half"),
     (f"{S}/FigS15_cor_GSE100159.pdf", "H", "half")],
])

# ---------- Figure S10: hub x CIBERSORT correlations ----------
build("FigureS10", [
    [(f"{S}/FigS4A_hub_ciber_cor_GSE76293.pdf", "A", "half"),
     (f"{S}/FigS4B_hub_ciber_cor_GSE69528.pdf", "B", "half")],
    [(f"{S}/FigS4C_hub_ciber_cor_GSE243066.pdf", "C", "half"),
     (f"{S}/FigS4D_hub_ciber_cor_GSE100159.pdf", "D", "half")],
])

# ---------- Figure S11: single-cell atlases ----------
build("FigureS11", [
    [(f"{C}/GSE151263_UMAP.pdf", "A", "half"),
     (f"{C}/GSE151263_celltype_prop.pdf", "B", "half")],
    [(f"{C}/GSE151263_marker_dotplot.pdf", "C", "half"),
     (f"{C}/GSE151263_hub_dotplot.pdf", "D", "half")],
    [(f"{C}/GSE167363/GSE167363_hub_dotplot.pdf", "E", "half"),
     (f"{S}/FigS9_hub_expr_heatmap_GSE167363.pdf", "F", "half")],
])

# ---------- Figure S12: CellChat comparison (整体拷贝) ----------
for ext in ["pdf", "png", "tiff"]:
    shutil.copy(f"{O}/FigureS10_CellChat_pathway_comparison.{ext}",
                f"{OUT}/FigureS12_composite.{ext}")
print("FigureS12 copied from FigureS10_CellChat_pathway_comparison")

# ---------- Figure S13: pseudotime / metabolism / knockout ----------
build("FigureS13", [
    [(f"{S}/FigS10A_AUCell_heatmap_GSE167363.pdf", "A", "half"),
     (f"{S}/FigS10B_CellOracle_TF_GSE167363.pdf", "B", "half")],
    [(f"{S}/FigS12_metabolic_delta.pdf", "C", "half"),
     (f"{S}/FigS13A_KO_concordance.pdf", "D", "half")],
    [(f"{S}/FigS13B_CEBPB_UPP1_gradient.pdf", "E", "half"),
     (f"{NP}/FigS14_donor_sensitivity.pdf", "F", "half")],
])

# ---------- Figure S14: protein & therapeutic validation ----------
build("FigureS14", [
    [(f"{S}/FigS7_HPA_infection_logFC.pdf", "A", "half"),
     (f"{O}/FigureS18_CMap_null_distribution.pdf", "B", "half")],
    [(f"{O}/FigureS19_Redocking_validation.pdf", "C", "half"),
     (f"{NP}/RS14D_alphafold_pLDDT.pdf", "D", "half")],
])

print("all supplementary composites done ->", OUT)
