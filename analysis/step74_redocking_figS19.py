#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 74 —— Figure S19: Redocking 验证图（编辑决策包 R4 验收标准）
- 重算 UPP1(3EUF/BAU) 与 VNN1(7SLV/9TC) 回对接各模式重原子 RMSD
  (原位、对称校正, 与 step56_redocking.py 同一算法)
- 两面板条形图: 各模式 RMSD + 2.0 A 阈值线 + 每模式亲和力标注
- 输出 PDF + PNG@300dpi + TIFF@1000dpi (LZW)
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit import Chem
from PIL import Image

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/09_docking_drug/redocking"
OUTDIR = ("/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/"
          "AJRCCM_submission/05_supplementary/Supplementary_Figures_Tables/Figures")

C_MAIN = "#6A5C95"   # Nature 5 色面板
C_ACC = "#D99AAE"
C_THR = "#E64B35"


# ---------- 与 step56_redocking.py 一致的 RMSD 计算 ----------
def _topo_heavy(mol):
    m = Chem.RemoveHs(Chem.Mol(mol))
    rw = Chem.RWMol(m)
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in rw.GetAtoms():
        a.SetIsAromatic(False)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
    out = rw.GetMol()
    Chem.SanitizeMol(out, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SYMMRINGS)
    return out


def pose_mols_from_pdbqt(path):
    from meeko import PDBQTMolecule, RDKitMolCreate
    pq = PDBQTMolecule.from_file(path, skip_typing=True)
    mols = RDKitMolCreate.from_pdbqt_mol(pq)
    if isinstance(mols, tuple):
        mols = mols[0]
    assert mols and mols[0] is not None, "姿态重建失败"
    return mols[0]


def rmsd_inplace(pose_mol, crystal_mol, conf_id):
    pt, ct = _topo_heavy(pose_mol), _topo_heavy(crystal_mol)
    assert pt.GetNumAtoms() == ct.GetNumAtoms()
    matches = ct.GetSubstructMatches(pt, uniquify=False, useChirality=False,
                                     useQueryQueryMatches=False, maxMatches=100000)
    assert matches, "拓扑图不匹配"
    keep_pose = [a.GetIdx() for a in pose_mol.GetAtoms() if a.GetSymbol() != "H"]
    keep_crys = [a.GetIdx() for a in crystal_mol.GetAtoms() if a.GetSymbol() != "H"]
    P = np.array(pose_mol.GetConformer(conf_id).GetPositions())[keep_pose]
    C = np.array(crystal_mol.GetConformer().GetPositions())[keep_crys]
    best = np.inf
    for match in matches:
        Cr = C[list(match)]
        r = float(np.sqrt(np.mean(np.sum((P - Cr) ** 2, axis=1))))
        best = min(best, r)
    return best


def affinities_from_pdbqt(path):
    affs = []
    with open(path) as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                affs.append(float(line.split()[3]))
    return affs


def per_mode_rmsd(case):
    pose = pose_mols_from_pdbqt(f"{BASE}/{case}_out.pdbqt")
    crystal = Chem.MolFromMolFile(f"{BASE}/{case}_crystal_lig.sdf", removeHs=False)
    if crystal is None:
        crystal = Chem.MolFromPDBFile(f"{BASE}/{case}_crystal_lig.pdb", removeHs=False)
    affs = affinities_from_pdbqt(f"{BASE}/{case}_out.pdbqt")
    n = pose.GetNumConformers()
    assert n == len(affs), f"{case}: conformers {n} vs affinities {len(affs)}"
    rmsds = [rmsd_inplace(pose, crystal, cid) for cid in range(n)]
    return rmsds, affs


# ---------- 绘图 ----------
CASES = [
    ("UPP1_3EUF", "UPP1 (PDB 3EUF) - BAU", "A"),
    ("VNN1_7SLV", "VNN1 (PDB 7SLV) - 9TC", "B"),
]

fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.1))
fig.subplots_adjust(left=0.09, right=0.98, bottom=0.14, top=0.86, wspace=0.28)

for ax, (case, title, panel) in zip(axes, CASES):
    rmsds, affs = per_mode_rmsd(case)
    modes = np.arange(1, len(rmsds) + 1)
    colors = [C_MAIN if i == 0 else C_ACC for i in range(len(rmsds))]
    bars = ax.bar(modes, rmsds, width=0.62, color=colors, edgecolor="black",
                  linewidth=0.6, zorder=3)
    ax.axhline(2.0, color=C_THR, linestyle="--", linewidth=1.1, zorder=4)
    ax.text(len(rmsds) + 0.42, 2.0, "2.0 " + chr(197), color=C_THR,
            fontsize=7.5, va="center", ha="left")
    for m, r, a in zip(modes, rmsds, affs):
        ax.text(m, r + 0.06, f"{r:.2f}", ha="center", fontsize=7.5)
        ax.text(m, -0.28, f"{a:.2f}", ha="center", fontsize=6.8, color="0.35")
    ax.set_xticks(modes)
    ax.set_xticklabels([f"Mode {m}" for m in modes], fontsize=7.5)
    ax.tick_params(axis="y", labelsize=7.5)
    top = max(max(rmsds) * 1.32, 2.5)
    ax.set_ylim(0, top)
    ax.set_xlim(0.3, len(rmsds) + 1.15)
    ax.set_ylabel("Heavy-atom RMSD (" + chr(197) + ")", fontsize=8.5)
    ax.set_title(f"{panel}  {title}", fontsize=9, loc="left", pad=16)
    ax.text(0.0, 1.04, "top-ranked pose", transform=ax.transAxes, fontsize=7,
            color=C_MAIN)
    ax.text(0.0, -0.20, "Vina affinity (kcal/mol) below each mode",
            transform=ax.transAxes, fontsize=6.8, color="0.35")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    print(f"{case}: RMSD per mode = "
          + ", ".join(f"{r:.3f}" for r in rmsds))

os.makedirs(OUTDIR, exist_ok=True)
stem = f"{OUTDIR}/FigureS19_Redocking_validation"
fig.savefig(f"{stem}.pdf")
fig.savefig(f"{stem}.png", dpi=300)
fig.savefig(f"{stem}.tiff", dpi=1000)
plt.close(fig)

# TIFF LZW 压缩 (PIL 中转)
img = Image.open(f"{stem}.tiff")
img.save(f"{stem}.tiff", compression="tiff_lzw")
print("saved:", stem + ".{pdf,png,tiff}")
