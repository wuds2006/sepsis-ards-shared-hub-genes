#!/usr/bin/env python3
"""真 PyMOL 渲染对接姿态图（用 conda 环境的 PyMOL）"""
import os, sys
import numpy as np

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/09_docking_drug"
DOCK = f"{BASE}/docking_site"
OUT = f"{BASE}/poses_pymol"
os.makedirs(OUT, exist_ok=True)

import pymol
from pymol import cmd

TOP = [
    ("VNN1", "luteolin", "cocrystal/7SLV.pdb", "A", (18.043, -28.315, -26.014)),
    ("VNN1", "quercetin", "cocrystal/7SLV.pdb", "A", (18.043, -28.315, -26.014)),
    ("UPP1", "curcumin", "cocrystal/3EUF.pdb", "A", (34.725, 55.011, -32.831)),
    ("GRB10", "EGCG", "receptors/GRB10_Q13322.pdb", None, None),
]

def pdbqt_ligand_to_pdb(pdbqt_path, out_pdb, model=1):
    """vina 配体输出 -> PDB"""
    lines = []
    cur = 0
    in_model = False
    n = 0
    for line in open(pdbqt_path):
        if line.startswith("MODEL"):
            cur = int(line.split()[1]); in_model = (cur == model); continue
        if line.startswith("ENDMDL"):
            if cur == model: break
            in_model = False; continue
        if in_model and (line.startswith("ATOM") or line.startswith("HETATM")):
            n += 1
            atom_name = line[12:16].strip()
            elem = atom_name[0] if atom_name else "C"
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            lines.append(f"HETATM{n:5d}  {atom_name:>3s} LIG A   1    "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}")
    with open(out_pdb, "w") as f:
        f.write("\n".join(lines) + "\nEND\n")

for gene, lig, pdb_path, chain, center in TOP:
    lig_file = f"{DOCK}/{gene}_{lig}_out.pdbqt"
    if not os.path.exists(lig_file):
        print(f"跳过 {gene}-{lig}: 缺对接输出"); continue
    # 配体转 PDB
    lig_pdb = f"{OUT}/{gene}_{lig}_ligand.pdb"
    pdbqt_ligand_to_pdb(lig_file, lig_pdb)

    # PyMOL 渲染
    cmd.reinitialize()
    cmd.load(f"{BASE}/{pdb_path}", "receptor")
    cmd.load(lig_pdb, "ligand")

    # 受体显示：cartoon + 结合口袋 stick
    cmd.hide("everything")
    cmd.show("cartoon", "receptor")
    cmd.color("lightblue", "receptor")
    # 结合口袋（配体周围 6Å 的受体残基）
    cmd.select("pocket", "receptor within 6 of ligand")
    cmd.show("sticks", "pocket")
    cmd.color("white", "pocket")
    cmd.util.cnc("pocket")
    # 配体显示
    cmd.show("sticks", "ligand")
    cmd.color("green", "ligand")
    cmd.set("stick_radius", 0.25, "ligand")
    cmd.set("sphere_scale", 0.3, "ligand")
    cmd.show("spheres", "ligand")

    # 视图 + 渲染
    cmd.zoom("ligand", 6)
    cmd.bg_color("white")
    cmd.ray(1600, 1200)
    out_png = f"{OUT}/{gene}_{lig}_pymol.png"
    cmd.png(out_png, dpi=300)
    print(f"  PyMOL 渲染: {out_png}")

cmd.quit()
print("===== PyMOL 渲染完成 =====")
