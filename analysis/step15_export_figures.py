#!/usr/bin/env python3
# AJRCCM figure export: enforce resolution and multi-format output.
# Color >=300 dpi, grayscale >=500 dpi, line art >=1000 dpi; emit PNG + TIFF + PDF.
import os, subprocess, glob
from PIL import Image

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT = f"{BASE}/AJRCCM_submission/07_figure_exports"

# (figure, panel, source_file, kind)  kind in {color, lineart, grayscale}
FIGS = [
    # Figure 1
    ("Figure1", "A_volcano_ARDS",       "01_DEG/GSE76293_volcano.pdf",              "color"),
    ("Figure1", "B_volcano_Sepsis",     "01_DEG/GSE69528_volcano.pdf",              "color"),
    ("Figure1", "C_moduleTrait_ARDS",   "02_WGCNA/GSE76293_module_trait.pdf",       "color"),
    ("Figure1", "D_moduleTrait_Sepsis", "02_WGCNA/GSE69528_module_trait.pdf",       "color"),
    ("Figure1", "E_venn4sets",          "03_intersection/venn_4sets.pdf",           "color"),
    ("Figure1", "F_sharedGO_dot",       "03_intersection/Shared_CGs_GO_dot.pdf",    "color"),
    # Figure 2
    ("Figure2", "A_ML_AUC_barplot",     "05_ML_SHAP/ML_AUC_barplot.pdf",            "color"),
    ("Figure2", "B_SHAP_beeplot",       "05_ML_SHAP/SHAP_beeplot.pdf",              "color"),
    ("Figure2", "C_ROC_Sepsis",         "06_external_validation/GSE100159_ROC.pdf", "color"),
    ("Figure2", "D_ROC_ARDS",           "06_external_validation/GSE243066_ROC.pdf", "color"),
    # Figure 3
    ("Figure3", "A_PPI_network",        "11_network_pharma/PPI_network.pdf",        "lineart"),
    ("Figure3", "B_GO_BP_dotplot",      "11_network_pharma/GO_BP_dotplot.pdf",      "color"),
    # Figure 4
    ("Figure4", "A_immune_box_ARDS",    "07_hub_immune/GSE76293_immune_boxplot.pdf","color"),
    ("Figure4", "B_immune_box_Sepsis",  "07_hub_immune/GSE69528_immune_boxplot.pdf","color"),
    ("Figure4", "C_hubHeat_ARDS",       "07_hub_immune/GSE76293_hub_immune_heatmap.pdf","color"),
    ("Figure4", "D_hubHeat_Sepsis",     "07_hub_immune/GSE69528_hub_immune_heatmap.pdf","color"),
    # Figure 5
    ("Figure5", "A_GSVA_top20_ARDS",    "07_hub_immune/gsva_top20/GSE76293_GSVA_heatmap_top20.png","color"),
    ("Figure5", "B_GSVA_top20_Sepsis",  "07_hub_immune/gsva_top20/GSE69528_GSVA_heatmap_top20.png","color"),
    # Figure 6
    ("Figure6", "A_UMAP",               "08_singlecell/GSE151263_UMAP.pdf",         "color"),
    ("Figure6", "B_cellchat_diff",      "08_singlecell/cellchat/diff_heatmap.pdf",  "color"),
    ("Figure6", "C_slingshot",          "08_singlecell/slingshot/slingshot_trajectory_UMAP.png","color"),
    # Figure 7
    ("Figure7", "A_docking_heatmap",    "09_docking_drug/docking_site_heatmap.pdf", "color"),
    ("Figure7", "B_VNN1_luteolin",      "09_docking_drug/poses_pymol/VNN1_luteolin_pymol.png","color"),
    ("Figure7", "C_VNN1_quercetin",     "09_docking_drug/poses_pymol/VNN1_quercetin_pymol.png","color"),
    ("Figure7", "D_GRB10_EGCG",         "09_docking_drug/poses_pymol/GRB10_EGCG_pymol.png","color"),
    ("Figure7", "E_UPP1_curcumin",      "09_docking_drug/poses_pymol/UPP1_curcumin_pymol.png","color"),
]

DPI = {"color": 300, "grayscale": 500, "lineart": 1000}

def render_pdf(pdf, out_prefix, dpi):
    """Render PDF to PNG and TIFF at given dpi via ghostscript."""
    png = f"{out_prefix}.png"
    tiff = f"{out_prefix}.tiff"
    subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH",
                    "-sDEVICE=png16m", f"-r{dpi}",
                    "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4",
                    f"-sOutputFile={png}", pdf], check=True)
    subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH",
                    "-sDEVICE=tiff24nc", f"-r{dpi}",
                    f"-sOutputFile={tiff}", pdf], check=True)
    # keep the vector PDF too
    subprocess.run(["cp", pdf, f"{out_prefix}.pdf"], check=True)
    return png

def process_png(src, out_prefix, dpi):
    """Open PNG, set DPI metadata, write PNG + TIFF + embedded PDF."""
    im = Image.open(src)
    if im.mode in ("RGBA", "P", "LA"):
        im = im.convert("RGB")
    png = f"{out_prefix}.png"
    tiff = f"{out_prefix}.tiff"
    pdf = f"{out_prefix}.pdf"
    im.save(png, dpi=(dpi, dpi))
    im.save(tiff, dpi=(dpi, dpi))
    im.save(pdf, "PDF", resolution=dpi)
    return png

count = 0
for fig, panel, src, kind in FIGS:
    srcp = f"{BASE}/{src}"
    if not os.path.exists(srcp):
        print(f"[SKIP] missing {src}")
        continue
    figdir = f"{OUT}/{fig}"
    os.makedirs(figdir, exist_ok=True)
    dpi = DPI[kind]
    out_prefix = f"{figdir}/{panel}"
    try:
        if src.endswith(".pdf"):
            render_pdf(srcp, out_prefix, dpi)
        else:
            process_png(srcp, out_prefix, dpi)
        print(f"[OK] {fig}/{panel} @ {dpi} dpi ({kind})")
        count += 1
    except Exception as e:
        print(f"[ERR] {fig}/{panel}: {e}")

print(f"\nDone: {count} panels exported to {OUT}")
