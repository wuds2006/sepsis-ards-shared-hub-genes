#!/usr/bin/env python3
"""py3Dmol 交互式 3D 对接姿态图（替代 PyMOL，浏览器可旋转查看）"""
import os, re
import numpy as np

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/09_docking_drug"
DOCK = f"{BASE}/docking_site"
OUT = f"{BASE}/poses_3dmol"
os.makedirs(OUT, exist_ok=True)

import py3Dmol

TOP = [
    ("VNN1", "luteolin"),
    ("VNN1", "quercetin"),
    ("GRB10", "EGCG"),
    ("UPP1", "curcumin"),
]

def pdbqt_to_pdb_ligand(pdbqt_path, model=1):
    """把 vina 输出的配体 PDBQT 转成 PDB 字符串（MODEL 1）"""
    lines = []
    cur = 0
    in_model = False
    for line in open(pdbqt_path):
        if line.startswith("MODEL"):
            cur = int(line.split()[1])
            in_model = (cur == model)
            continue
        if line.startswith("ENDMDL"):
            if cur == model:
                break
            in_model = False
            continue
        if in_model and (line.startswith("ATOM") or line.startswith("HETATM")):
            # PDBQT -> PDB：改 HETATM 保留坐标
            elem = line[77:79].strip() or line[12:16].strip()[0]
            atom_name = line[12:16].strip()
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            # 标准 PDB ATOM 格式
            lines.append(f"HETATM{len(lines)+1:5d}  {atom_name:>3s} LIG A   1    "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}")
    return "\n".join(lines)

def extract_pocket_pdb(pdb_path, chain, center, radius=8.0):
    """从共结晶 PDB 提取结合口袋蛋白原子（配体质心半径内）"""
    lines = []
    for line in open(pdb_path):
        if line.startswith("ATOM") and line[21:22].strip() == chain:
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            d = np.sqrt((x-center[0])**2 + (y-center[1])**2 + (z-center[2])**2)
            if d < radius:
                lines.append(line.rstrip())
    return "\n".join(lines)

def make_html(protein_pdb, ligand_pdb, title, out_html):
    view = py3Dmol.view(width=800, height=600)
    view.addModel(protein_pdb, "pdb")
    view.setStyle({"model": 0}, {"cartoon": {"color": "#CFE3E6"}})
    view.addStyle({"model": 0}, {"stick": {"radius": 0.15, "colorscheme": "whiteCarbon"}})
    view.addModel(ligand_pdb, "pdb")
    view.setStyle({"model": 1}, {"stick": {"radius": 0.3, "colorscheme": "greenCarbon"}})
    view.addStyle({"model": 1}, {"sphere": {"scale": 0.35, "colorscheme": "greenCarbon"}})
    view.zoomTo({"model": 1})
    view.setBackgroundColor("white")
    html = view._make_html().replace("</head>", f"<title>{title}</title></head>")
    # 加标题
    html = html.replace("<body>", f"<body><h3 style='font-family:sans-serif'>{title}</h3>")
    with open(out_html, "w") as f:
        f.write(html)
    print(f"  交互式 3D: {out_html}")

COC_PDB = {
    "VNN1": ("cocrystal/7SLV.pdb", "A", (18.043, -28.315, -26.014)),
    "UPP1": ("cocrystal/3EUF.pdb", "A", (34.725, 55.011, -32.831)),
}

for gene, lig in TOP:
    lig_file = f"{DOCK}/{gene}_{lig}_out.pdbqt"
    if not os.path.exists(lig_file):
        print(f"跳过 {gene}-{lig}: 缺对接输出"); continue
    ligand_pdb = pdbqt_to_pdb_ligand(lig_file)
    if gene in COC_PDB:
        pdb, chain, center = COC_PDB[gene]
        protein_pdb = extract_pocket_pdb(f"{BASE}/{pdb}", chain, center)
    else:
        # 非共结晶基因，用 AlphaFold 受体 PDBQT 转 PDB 的简化口袋
        protein_pdb = ""
    title = f"{gene} – {lig} 结合姿态（共结晶位点，交互式）"
    make_html(protein_pdb, ligand_pdb, title, f"{OUT}/{gene}_{lig}_3dmol.html")

print("===== py3Dmol 交互式 3D 姿态图完成 =====")
