#!/usr/bin/env python
"""Assemble supplementary composite figures S11 (Monocle3), S12 (metabolism),
S13 (group-stratified knockout + CellOracle).
Style follows scripts/step25_assemble_supplementary.py: PIL, white bg, letter labels,
300dpi triple export.
"""
from PIL import Image, ImageDraw, ImageFont
import os

ROOT = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
SC = os.path.join(ROOT, "08_singlecell")
MOC = os.path.join(SC, "monocle3")
MET = os.path.join(SC, "metabolism")
KO = os.path.join(SC, "knockout_stratified")
COS = os.path.join(SC, "celloracle_stratified")
SUB = os.path.join(ROOT, "AJRCCM_submission/05_supplementary/Supplementary_Figures_Tables/Figures")

GAP = 50
MARGIN = 60
LETTER_SZ = 80
font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", LETTER_SZ)


def load(src, target_w):
    im = Image.open(src)
    if im.mode in ("RGBA", "LA", "P"):
        im = im.convert("RGBA")
        bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
        im = Image.alpha_composite(bg, im).convert("RGB")
    else:
        im = im.convert("RGB")
    w, h = im.size
    return im.resize((target_w, int(round(h * target_w / w))), Image.LANCZOS)


def assemble(rows_spec, out_stem):
    composed = []
    for row in rows_spec:
        cells = []
        for src, letter, cw in row:
            cells.append((load(src, cw), letter))
        composed.append(cells)
    W = max(sum(c.size[0] for c, _ in row) + GAP * (len(row) - 1) for row in composed)
    Hh = MARGIN * 2 + sum(max(c.size[1] for c, _ in row) for row in composed) + GAP * (len(composed) - 1)
    canvas = Image.new("RGB", (W + MARGIN * 2, Hh), "white")
    d = ImageDraw.Draw(canvas)
    y = MARGIN
    for row in composed:
        h = max(c.size[1] for c, _ in row)
        row_w = sum(c.size[0] for c, _ in row) + GAP * (len(row) - 1)
        x = MARGIN + (W - row_w) // 2
        for c, letter in row:
            canvas.paste(c, (x, y + (h - c.size[1]) // 2))
            d.rectangle([x, y, x + 105, y + 105], fill="white")
            d.text((x + 12, y + 6), letter, fill="black", font=font)
            x += c.size[0] + GAP
        y += h + GAP
    for ext, kw in [("png", dict(dpi=(300, 300))),
                    ("tiff", dict(compression="tiff_lzw", dpi=(300, 300))),
                    ("pdf", dict(resolution=300))]:
        out = f"{out_stem}.{ext}"
        canvas.save(out, **kw)
        print("saved", out, canvas.size)


FULL = 3050
CELL3 = (FULL - 2 * GAP) // 3

# ---------- Figure S11: Monocle3 trajectories, both cohorts ----------
assemble([
    [(os.path.join(MOC, "GSE151263_monocle3_trajectory_celltype.png"), "A", CELL3),
     (os.path.join(MOC, "GSE151263_monocle3_pseudotime_group.png"), "B", CELL3),
     (os.path.join(MOC, "GSE151263_monocle3_hub_moransI_bar.png"), "C", CELL3)],
    [(os.path.join(MOC, "GSE167363_monocle3_trajectory_celltype.png"), "D", CELL3),
     (os.path.join(MOC, "GSE167363_monocle3_pseudotime_group.png"), "E", CELL3),
     (os.path.join(MOC, "GSE167363_monocle3_hub_moransI_bar.png"), "F", CELL3)],
], os.path.join(SUB, "FigureS11_monocle3_trajectories"))

# ---------- Figure S12: metabolic pathway activity, both cohorts ----------
assemble([
    [(os.path.join(MET, "GSE151263_metab_pathway_heatmap.png"), "A", CELL3),
     (os.path.join(MET, "GSE151263_metab_group_difference.png"), "B", CELL3),
     (os.path.join(MET, "GSE151263_metab_hub_correlation.png"), "C", CELL3)],
    [(os.path.join(MET, "GSE167363_metab_pathway_heatmap.png"), "D", CELL3),
     (os.path.join(MET, "GSE167363_metab_group_difference.png"), "E", CELL3),
     (os.path.join(MET, "GSE167363_metab_hub_correlation.png"), "F", CELL3)],
], os.path.join(SUB, "FigureS12_metabolic_pathways"))

# ---------- Figure S13: group-stratified KO + CellOracle ----------
assemble([
    [(os.path.join(KO, "GSE151263_KO_stratified_GO_BP_UPP1.png"), "A", CELL3),
     (os.path.join(KO, "GSE151263_KO_stratified_top10_UPP1.png"), "B", CELL3),
     (os.path.join(COS, "GSE151263_CellOracle_stratified_hubdelta.png"), "C", CELL3)],
    [(os.path.join(KO, "GSE167363_KO_stratified_GO_BP_CDK5RAP2.png"), "D", CELL3),
     (os.path.join(KO, "GSE167363_KO_stratified_top10_CDK5RAP2.png"), "E", CELL3),
     (os.path.join(COS, "GSE167363_CellOracle_stratified_hubdelta.png"), "F", CELL3)],
], os.path.join(SUB, "FigureS13_stratified_knockout"))

print("DONE")
