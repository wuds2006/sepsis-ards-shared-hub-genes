#!/usr/bin/env python3
"""step112: export the missing PNG / TIFF companions for figures that currently exist
only as PDF, matching the resolution convention already used in this package
(PNG 300 dpi, TIFF 1000 dpi, LZW compression).

Targets (verified 2026-10-19):
  figures/main/          Figure1_composite, Figure2_composite, Figure10_composite
  figures/supplementary/ FigureS2_*, FigureS4_*, FigureS5_*, FigureS17_*
"""
import glob
import os
import sys

import fitz  # PyMuPDF
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
FIGDIR = os.path.join(ROOT, "AJRCCM_submission/09_review_deAI_2026-09-25/figures")

TARGETS = [
    ("main", "Figure1_composite"),
    ("main", "Figure2_composite"),
    ("main", "Figure10_composite"),
    ("supplementary", "FigureS2_"),
    ("supplementary", "FigureS4_"),
    ("supplementary", "FigureS5_"),
    ("supplementary", "FigureS17_"),
]


def resolve(sub, stem):
    hits = sorted(glob.glob(os.path.join(FIGDIR, sub, stem + "*.pdf")))
    return hits[0] if hits else None


def export(pdf_path, sub):
    base = os.path.splitext(pdf_path)[0]
    doc = fitz.open(pdf_path)
    n = len(doc)
    for res, ext, save in ((300, "png", "png"), (1000, "tiff", "tiff_lzw")):
        out = f"{base}.{ext}"
        if os.path.exists(out):
            print("  exists, skipped:", os.path.basename(out))
            continue
        zoom = res / 72.0
        mat = fitz.Matrix(zoom, zoom)
        if n == 1:
            pix = doc[0].get_pixmap(matrix=mat)
            mode = "RGB" if pix.n < 4 else "RGBA"
            img = Image.frombytes(mode, (pix.width, pix.height), pix.samples)
            if save == "tiff_lzw":
                img.save(out, format="TIFF", dpi=(res, res), compression="tiff_lzw")
            else:
                img.save(out, format="PNG", dpi=(res, res))
            print(f"  wrote {os.path.basename(out)}  {img.size}  @{res} dpi  "
                  f"{os.path.getsize(out)/1e6:.2f} MB")
        else:
            # multi-page PDF: export page 1 with a suffix, and note the page count
            pix = doc[0].get_pixmap(matrix=mat)
            mode = "RGB" if pix.n < 4 else "RGBA"
            img = Image.frombytes(mode, (pix.width, pix.height), pix.samples)
            out_p1 = f"{base}_p1.{ext}"
            if save == "tiff_lzw":
                img.save(out_p1, format="TIFF", dpi=(res, res), compression="tiff_lzw")
            else:
                img.save(out_p1, format="PNG", dpi=(res, res))
            print(f"  [multi-page PDF: {n} pages] wrote {os.path.basename(out_p1)}  "
                  f"{img.size}  @{res} dpi")


for sub, stem in TARGETS:
    pdf = resolve(sub, stem)
    print(f"== {sub} / {stem}")
    if not pdf:
        print("  !! no PDF found")
        continue
    export(pdf, sub)

print("\ndone")
