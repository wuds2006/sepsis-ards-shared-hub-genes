#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 88 —— 清除 main_panels 中面板标题的旧字母前缀（"A  "、"E2  "、"C1  " 等）
就地编辑 PDF：定位 Helvetica-Bold 标题 span → 删除前缀 → 重写标题文本。
只处理 RESULTS 图组用到的面板。
"""
import re, sys
import pymupdf

P = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission/03_figures/Figures_ggplot2_NPG/main_panels"

FILES = [
    "Fig2A_boxplot_GSE76293.pdf", "Fig2B_boxplot_GSE243066.pdf",
    "Fig2C_boxplot_GSE69528.pdf", "Fig2D_boxplot_GSE100159.pdf",
    "Fig2E1_ROC_GSE76293.pdf", "Fig2E2_ROC_GSE243066.pdf",
    "Fig2E3_ROC_GSE69528.pdf", "Fig2E4_ROC_GSE100159.pdf",
    "Fig2F1_ML_AUC_barplot.pdf", "Fig2F2_SHAP_barplot.pdf",
    "Fig4C_ROC_cal_DCA_ARDS.pdf", "Fig4D_ROC_cal_DCA_sepsis.pdf",
    "Fig6A_GSVA_top20_GSE76293.pdf", "Fig6B_GSVA_top20_GSE69528.pdf",
    "Fig7A_UMAP_GSE151263.pdf", "Fig7C_trajectory_AUCell_GSE151263.pdf",
    "Fig7D_CellOracle_GSE151263.pdf", "Fig7E_UMAP_GSE167363.pdf",
    "Fig7G_trajectory_AUCell_GSE167363.pdf", "Fig7H_CellOracle_GSE167363.pdf",
    "Fig8A_docking_site_heatmap.pdf", "Fig8C_glucocorticoid_repositioning.pdf",
]

PREFIX = re.compile(r"^[A-H][0-9]?\s\s+")

total = 0
for fname in FILES:
    path = f"{P}/{fname}"
    d = pymupdf.open(path)
    page = d[0]
    jobs = []
    for b in page.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            for s in l["spans"]:
                if "Bold" in s["font"] and PREFIX.match(s["text"]):
                    new_text = PREFIX.sub("", s["text"]).replace("\u2212", "-")
                    jobs.append((s["bbox"], s["origin"], s["size"], new_text, s["text"]))
    if not jobs:
        d.close()
        print(f"  {fname}: no prefixed titles")
        continue
    for bbox, origin, size, new_text, old_text in jobs:
        # 红action 删除旧标题（白色填充）
        page.add_redact_annot(pymupdf.Rect(bbox), fill=(1, 1, 1))
    page.apply_redactions()
    for bbox, origin, size, new_text, old_text in jobs:
        page.insert_text(pymupdf.Point(origin), new_text, fontsize=size,
                         fontname="hebo", color=(0, 0, 0))
        total += 1
        print(f"  {fname}: '{old_text}' -> '{new_text}'")
    tmp = path + ".tmp"
    d.save(tmp, deflate=True)
    d.close()
    import os, shutil
    shutil.move(tmp, path)
print(f"step88 done: {total} titles rewritten")
