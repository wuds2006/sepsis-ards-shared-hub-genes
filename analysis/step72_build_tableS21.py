#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 72 —— 合并 MM/GBSA 结果为 Supplementary Table S21 (xlsx)
输入: 09_docking_drug/mmgb_{GENE}.csv (step70 产出, 含 vina 列)
输出: AJRCCM_submission/05_supplementary/Supplementary_Figures_Tables/Tables/
      TableS21_MMGBSA_free_energies.xlsx
  Sheet1 "MMGBSA": 全部 50 组合逐行 + 按靶点/配体的 ΔG 排名
  Sheet2 "Summary": 每靶点统计 (配体数, ΔG 范围/中位数, Spearman ρ 与 Vina 的一致性)
  Sheet3 "Method note": 协议与限制说明 (含 openforcefield 供应链事件披露)
"""
import csv, os
import numpy as np
from scipy.stats import spearmanr
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/09_docking_drug"
OUT = ("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/"
       "AJRCCM_submission/05_supplementary/Supplementary_Figures_Tables/Tables/"
       "TableS21_MMGBSA_free_energies.xlsx")

GENES = ["VNN1", "UPP1", "CD177", "EXOSC4", "GRB10"]  # CDK5RAP2 按 pLDDT 73.5 排除

rows = []
for g in GENES:
    p = f"{BASE}/mmgb_{g}.csv"
    if not os.path.exists(p):
        print("missing:", p); continue
    for r in csv.DictReader(open(p)):
        r["gene"] = g
        rows.append(r)

if not rows:
    raise SystemExit("no mmgb CSVs found")

for r in rows:
    for k in ("E_complex", "E_protein", "E_ligand", "dSASA", "deltaG", "vina"):
        r[k] = float(r[k])

# 全局排名 (越负越强)
order = sorted(rows, key=lambda r: r["deltaG"])
for rank, r in enumerate(order, 1):
    r["rank_global"] = rank

wb = Workbook()
ws = wb.active
ws.title = "MMGBSA"
hdr = ["rank_global", "gene", "ligand", "E_complex_kcal_mol", "E_protein_kcal_mol",
       "E_ligand_kcal_mol", "dSASA_A2", "deltaG_kcal_mol", "Vina_kcal_mol"]
ws.append(hdr)
for c in range(1, len(hdr) + 1):
    ws.cell(1, c).font = Font(bold=True)
for r in sorted(rows, key=lambda x: x["deltaG"]):
    ws.append([r["rank_global"], r["gene"], r["ligand"], r["E_complex"],
               r["E_protein"], r["E_ligand"], r["dSASA"], r["deltaG"], r["vina"]])
for c, w in zip(range(1, 10), [10, 10, 18, 16, 16, 15, 11, 14, 13]):
    ws.column_dimensions[get_column_letter(c)].width = w

ws2 = wb.create_sheet("Summary")
ws2.append(["target", "n_ligands", "dG_best", "dG_best_ligand", "dG_median",
            "dG_range", "Spearman_rho_dG_vs_Vina", "Spearman_P"])
for c in range(1, 9):
    ws2.cell(1, c).font = Font(bold=True)
summary_all = []
for g in GENES:
    sub = [r for r in rows if r["gene"] == g]
    if not sub:
        continue
    dg = np.array([r["deltaG"] for r in sub])
    vi = np.array([r["vina"] for r in sub])
    best = sub[int(np.argmin(dg))]
    if len(sub) >= 3 and np.std(dg) > 0 and np.std(vi) > 0:
        rho, pv = spearmanr(-dg, -vi)   # 均转为"越正越强"方向
    else:
        rho, pv = float("nan"), float("nan")
    ws2.append([g, len(sub), dg.min(), best["ligand"], float(np.median(dg)),
                f"{dg.min():.1f} to {dg.max():.1f}",
                round(float(rho), 3) if rho == rho else "NA",
                f"{pv:.3g}" if pv == pv else "NA"])
    summary_all.extend([(r["deltaG"], r["vina"]) for r in sub])
dgA = np.array([x[0] for x in summary_all]); viA = np.array([x[1] for x in summary_all])
rhoA, pvA = spearmanr(-dgA, -viA)
ws2.append([])
ws2.append(["ALL (pooled)", len(dgA), dgA.min(), "", float(np.median(dgA)),
            f"{dgA.min():.1f} to {dgA.max():.1f}", round(float(rhoA), 3), f"{pvA:.3g}"])
ws2.cell(ws2.max_row, 1).font = Font(bold=True)
for c, w in zip(range(1, 9), [12, 10, 9, 16, 11, 14, 22, 12]):
    ws2.column_dimensions[get_column_letter(c)].width = w

ws3 = wb.create_sheet("Method_note")
notes = [
    "Supplementary Table S21. Single-point MM/GBSA free-energy rescoring of docked complexes.",
    "",
    "Protocol: OpenMM 8.5.1; protein ff14SB with GBn2 (OBC-II) implicit solvent, NoCutoff;",
    "ligands with Gasteiger partial charges and UFF Lennard-Jones parameters (Lennard-Jones",
    "terms of polar hydrogens zeroed, following the Amber convention; ff14SB HO type has",
    "epsilon = 0); generalized Born radii from the standard OBC parameter assignment.",
    "Each complex was relaxed by local energy minimization with the ligand held by positional",
    "restraints (5,000 kJ/mol/nm^2); a short-range repulsive wall (15 kcal/mol x (1.5 A / r)^12,",
    "ligand-involving pairs only, applied during minimization only) prevented unphysical",
    "Coulombic collapse of charge pairs in the absence of a bonded ligand force field.",
    "DeltaG = E(complex) - E(protein) - E(ligand) + 0.005 kcal/mol/A^2 x dSASA",
    "(Shrake-Rupley, 1.4 A probe), single-trajectory convention on the minimized coordinates.",
    "",
    "Scope: 5 targets x 10 natural-product ligands = 50 combinations. CDK5RAP2 was excluded",
    "because its docking-box model confidence was moderate (mean pLDDT 73.5; Supplementary",
    "Table S19); its Vina docking energies are reported as exploratory only.",
    "",
    "Limitations: GAFF/AM1-BCC ligand parameterization was unavailable in this environment",
    "(the openforcefield toolchain was removed from PyPI and conda-forge in 2026 following a",
    "supply-chain incident; AmberTools could not be retrieved), so UFF/Gasteiger parameters",
    "were substituted. These are single-point estimates on a single docked pose per",
    "combination, without conformational sampling or entropy; values are intended for",
    "internal ranking and consistency checks with the Vina scores, not as calibrated",
    "absolute binding free energies. E_ligand contains no bonded terms and is reported only",
    "for accounting; the single-trajectory convention cancels all bonded contributions.",
]
for line in notes:
    ws3.append([line])
ws3.column_dimensions["A"].width = 100

os.makedirs(os.path.dirname(OUT), exist_ok=True)
wb.save(OUT)
print(f"saved {OUT}: {len(rows)} combinations across {len(set(r['gene'] for r in rows))} targets")
print(f"pooled Spearman rho(dG vs Vina, direction-aligned) = {rhoA:.3f}, P = {pvA:.3g}")
for g in GENES:
    sub = [r for r in rows if r["gene"] == g]
    if sub:
        best = min(sub, key=lambda r: r["deltaG"])
        print(f"  {g}: best {best['ligand']} dG={best['deltaG']:.1f}")
