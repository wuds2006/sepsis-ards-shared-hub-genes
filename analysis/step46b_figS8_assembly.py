#!/usr/bin/env python3
"""Wrap the rebuilt four-cohort single-gene ssGSEA grid into FigureS8_composite.pdf
(178 mm page, vector-preserving embed), replacing the original two-cohort render."""
import pymupdf

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission"
PANEL = f"{BASE}/03_figures/Figures_ggplot2_NPG/supp_panels/FigureS8_4cohorts.pdf"
OUT = f"{BASE}/03_figures/Figures_ggplot2_NPG/composites/FigureS8_composite.pdf"

MM = 72 / 25.4
MARGIN_X = 8.0
TOP = 8.0
W = 178 - 2 * MARGIN_X  # 162 mm usable

src = pymupdf.open(PANEL)
r = src[0].rect
h_mm = W * r.height / r.width

doc = pymupdf.open()
page = doc.new_page(width=178 * MM, height=(TOP + h_mm + 6) * MM)
dest = pymupdf.Rect(MARGIN_X * MM, TOP * MM, (MARGIN_X + W) * MM, (TOP + h_mm) * MM)
page.show_pdf_page(dest, src, 0)
doc.save(OUT, deflate=True, garbage=3)
src.close()
doc.close()
print("saved", OUT, f"(178x{TOP + h_mm + 6:.0f}mm)")
