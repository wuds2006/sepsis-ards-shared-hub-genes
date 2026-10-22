#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 69 —— VNN1-luteolin 相互作用分析（交叉学科专家 W2(b) 建议 + 小问题裁定）
1. 从 docking_site/VNN1_luteolin_out.pdbqt 的 MODEL 1 (-10.006 kcal/mol) 重构配体姿态（meeko）
2. 7SLV 链 A 提取蛋白 + pdbfixer 加氢 (pH 7.4)
3. PLIP 全类型相互作用检测（氢键/疏水/π-堆积/π-阳离子/盐桥/水桥）
4. LigPlot 风格 2D 相互作用示意图（matplotlib, 矢量 PDF）
输出: 09_docking_drug/plip_vnn1_luteolin.csv + plip_vnn1_luteolin_2D.pdf/.png
"""
import os, csv
import numpy as np

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/09_docking_drug"
SITE = f"{BASE}/docking_site"
CO   = f"{BASE}/cocrystal/7SLV.pdb"
OUT  = BASE

# ---------- 1) 配体姿态重构 ----------
from meeko import PDBQTMolecule, RDKitMolCreate
pq = PDBQTMolecule.from_file(f"{SITE}/VNN1_luteolin_out.pdbqt", skip_typing=True)
mols = RDKitMolCreate.from_pdbqt_mol(pq)
lig = mols[0]  # MODEL 1 = 最优姿态
conf = lig.GetConformer()
lig_pos = np.array(conf.GetPositions())
print("ligand atoms:", lig.GetNumAtoms(), "| heavy:", lig.GetNumHeavyAtoms())

# ---------- 2) 蛋白提取 + 加氢 ----------
from pdbfixer import PDBFixer
from openmm.app import PDBFile
prot_lines = []
for line in open(CO):
    if line.startswith(("ATOM", "HETATM")) and line[21] == "A":
        resn = line[17:20].strip()
        if resn in {"HOH", "9TC"}:  # 去水 + 去共晶配体（与对接受体一致）
            continue
        prot_lines.append(line)
chainA = f"{OUT}/_vnn1_chainA_raw.pdb"
with open(chainA, "w") as f:
    f.writelines(prot_lines)
fixer = PDBFixer(filename=chainA)
fixer.findMissingResidues()
fixer.missingResidues = {}          # 不补末端缺失（保持晶体状态）
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.addMissingHydrogens(7.4)
prot_pdb = f"{OUT}/_vnn1_chainA_H.pdb"
with open(prot_pdb, "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
print("protein prepared:", prot_pdb)

# ---------- 3) 组装复合物 PDB 字符串 ----------
complex_lines = [l for l in open(prot_pdb) if l.startswith(("ATOM", "HETATM"))]
elements = [a.GetSymbol().upper() for a in lig.GetAtoms()]
for i, p in enumerate(lig_pos):
    el = elements[i]
    nm = {"C": " C", "N": " N", "O": " O", "H": " H"}.get(el, f" {el}")
    complex_lines.append(
        "HETATM%5d %4s LUT B%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
        % (9000 + i, ("H" if el == "H" else nm.strip()), 1,
           p[0], p[1], p[2], el.rjust(2)))
pdb_string = "END\n".join([""] * 0) + "".join(complex_lines) + "END\n"

# ---------- 4) PLIP 检测 ----------
from plip.structure.preparation import PDBComplex
pc = PDBComplex()
pc.load_pdb(pdb_string, as_string=True)
pc.analyze()
keys = list(pc.interaction_sets.keys())
print("PLIP ligand keys:", keys)
lut_keys = [k for k in keys if k.startswith("LUT")]
assert lut_keys, "LUT ligand not found by PLIP"
site = pc.interaction_sets[lut_keys[0]]

def dist_of(o):
    for f in ("distance", "distance_ad", "distance_ah", "dist", "center_distance"):
        v = getattr(o, f, None)
        if v is not None:
            return round(v, 2)
    return ""

records = []
for src_list, side_label in [(site.hbonds_ldon, "ligand donor"), (site.hbonds_pdon, "ligand acceptor")]:
    for hb in src_list:
        records.append(dict(type="hydrogen bond", residue=hb.restype + str(hb.resnr),
                            distance=dist_of(hb), side=side_label, detail=""))
for hy in site.hydrophobic_contacts:
    records.append(dict(type="hydrophobic", residue=hy.restype + str(hy.resnr),
                        distance=dist_of(hy), side="", detail=""))
for ps in site.pistacking:
    records.append(dict(type="pi-stacking", residue=ps.restype + str(ps.resnr),
                        distance=dist_of(ps), side=ps.type, detail=""))
for src_list, side_label in [(site.pication_laro, "protein cation"), (site.pication_paro, "ligand cation")]:
    for pcn in src_list:
        records.append(dict(type="pi-cation", residue=pcn.restype + str(pcn.resnr),
                            distance=dist_of(pcn), side=side_label, detail=""))
for sb in site.saltbridge_lneg + site.saltbridge_pneg:
    records.append(dict(type="salt bridge", residue=sb.restype + str(sb.resnr),
                        distance=dist_of(sb), side="", detail=""))
for wb in site.water_bridges:
    records.append(dict(type="water bridge", residue=wb.restype + str(wb.resnr),
                        distance=dist_of(wb), side="", detail=""))

with open(f"{OUT}/plip_vnn1_luteolin.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["type", "residue", "distance", "side", "detail"])
    w.writeheader(); w.writerows(records)
from collections import Counter
print("interaction summary:", dict(Counter(r["type"] for r in records)))
for r in records:
    print("  ", r["type"], r["residue"], r["distance"])

# ---------- 5) LigPlot 风格 2D 相互作用图 ----------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from rdkit.Chem import AllChem, rdDepictor
from rdkit import Chem

rdDepictor.SetPreferCoordGen(True)
lig2d = Chem.Mol(lig)
AllChem.Compute2DCoords(lig2d)
c2 = lig2d.GetConformer()
pos2d = np.array([list(c2.GetAtomPosition(i))[:2] for i in range(lig2d.GetNumAtoms())])
# 与 3D 一致的原子索引（副本拓扑相同）
elem2d = [a.GetSymbol() for a in lig2d.GetAtoms()]

TYPE_COLOR = {
    "hydrogen bond": "#00A087", "hydrophobic": "#7E6148", "pi-stacking": "#E64B35",
    "pi-cation": "#8491B4", "salt bridge": "#3C5488", "water bridge": "#4DBBD5",
}
TYPE_Z = {"salt bridge": 0, "pi-stacking": 1, "pi-cation": 2, "hydrogen bond": 3, "water bridge": 4, "hydrophobic": 5}

# 每个残基取最高优先级类型 + 涉及的配体原子
# 需要从 PLIP 对象拿配体原子索引；hbond 等记录无 idx，改用最近原子法：
prot_atoms = {}
for line in complex_lines:
    if line.startswith("ATOM"):
        res = line[17:20].strip() + line[22:26].strip()
        prot_atoms.setdefault(res, []).append(
            (float(line[30:38]), float(line[38:46]), float(line[46:54])))
centroid = lig_pos.mean(axis=0)
res_anchor = {}
for r in records:
    if r["residue"] not in prot_atoms:
        continue
    pts = np.array(prot_atoms[r["residue"]])
    d = np.linalg.norm(pts[:, None, :] - lig_pos[None, :, :], axis=2)
    amin = np.unravel_index(np.argmin(d), d.shape)
    la, dist = amin[1], d[amin]
    if r["residue"] not in res_anchor or dist < res_anchor[r["residue"]][1]:
        res_anchor[r["residue"]] = (la, dist)
    if r["residue"] in res_anchor:
        res_anchor[r["residue"]] = (
            la if dist < res_anchor[r["residue"]][1] else res_anchor[r["residue"]][0],
            min(dist, res_anchor[r["residue"]][1]))

# 残基类型优先级着色
res_type = {}
for r in records:
    if r["residue"] not in res_type or TYPE_Z[r["type"]] < TYPE_Z[res_type[r["residue"]]]:
        res_type[r["residue"]] = r["type"]

fig, ax = plt.subplots(figsize=(7.0, 5.6))
plt.rcParams.update({"font.family": "Arial", "font.size": 8})
# 画配体键
bondset = set()
for b in lig2d.GetBonds():
    i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
    bondset.add((min(i, j), max(i, j)))
for i, j in bondset:
    ax.plot(*zip(pos2d[i], pos2d[j]), color="#515151", lw=1.4, zorder=2, solid_capstyle="round")
# 画配体重原子
for i, (x, y) in enumerate(pos2d):
    if elem2d[i] == "H":
        continue
    color = {"O": "#E64B35", "N": "#3C5488", "C": "#F5F5F5"}.get(elem2d[i], "#F5F5F5")
    edge = "#515151" if elem2d[i] != "C" else "#9E9E9E"
    ax.add_patch(plt.Circle((x, y), 0.14, facecolor=color, edgecolor=edge, lw=1.0, zorder=3))
    if elem2d[i] != "C":
        ax.text(x, y, elem2d[i], ha="center", va="center", fontsize=7,
                color="white" if elem2d[i] != "C" else "#333", zorder=4, fontweight="bold")

# 画残基辐条（标签置于分子外围均匀分布，避免重叠）
cen2d = pos2d.mean(axis=0)
span = max(pos2d[:, 0].max() - pos2d[:, 0].min(), pos2d[:, 1].max() - pos2d[:, 1].min())
R = span / 2 + 1.1        # 标签半径
spokes = []
for res, (la, _) in res_anchor.items():
    if res not in res_type:
        continue
    a = pos2d[la]
    v = a - cen2d
    ang = np.arctan2(v[1], v[0])
    spokes.append((ang, res, la))
spokes.sort()
# 角度去拥挤：相邻标签至少 0.5 rad
MIN_A = 0.55
angs = [s[0] for s in spokes]
for i in range(1, len(spokes)):
    if angs[i] - angs[i - 1] < MIN_A:
        angs[i] = angs[i - 1] + MIN_A
for i in range(len(spokes) - 2, -1, -1):   # 反向再校一次（环绕 2pi）
    gap = angs[i + 1] - angs[i]
    if gap < MIN_A:
        angs[i] = angs[i + 1] - MIN_A
for (ang0, res, la), ang in zip(spokes, angs):
    t = res_type[res]
    a = pos2d[la]
    v = np.array([np.cos(ang), np.sin(ang)])
    p0 = a + v * 0.28
    p1 = cen2d + v * R
    ax.plot(*zip(p0, p1), ls=(0, (3, 3)), color=TYPE_COLOR[t], lw=1.0, zorder=1)
    ax.text(p1[0], p1[1], res, ha="center", va="center", fontsize=7.5,
            color=TYPE_COLOR[t], fontweight="bold", zorder=5)

legend_elems = [Line2D([0], [0], color=c, lw=1.4, ls="--", label=t)
                for t, c in TYPE_COLOR.items() if t in set(res_type.values())]
ax.legend(handles=legend_elems, loc="upper left", bbox_to_anchor=(-0.06, 1.02),
          fontsize=7, frameon=False, title="interaction", title_fontsize=7)
ax.set_title("VNN1 (7SLV) – luteolin: predicted interactions", fontsize=9, loc="left")
ax.set_aspect("equal"); ax.axis("off")
pad = R + 1.4
ax.set_xlim(cen2d[0] - pad, cen2d[0] + pad)
ax.set_ylim(cen2d[1] - pad, cen2d[1] + pad)
plt.tight_layout()
plt.savefig(f"{OUT}/plip_vnn1_luteolin_2D.pdf")
plt.savefig(f"{OUT}/plip_vnn1_luteolin_2D.png", dpi=300)
print("saved 2D interaction figure")
