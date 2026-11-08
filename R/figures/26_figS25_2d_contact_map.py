#!/usr/bin/env python3
"""把 PLIP 接触锚定到配体重原子，并导出 RDKit 2D 坐标供 R 绘图。
输出: contact_atoms.csv (residue, type, lig_atom_idx0, dist), atoms2d.csv, bonds.csv
"""
import os, json
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor
from rdkit.Geometry import Point3D

BASE = os.environ.get("DOCK_ROOT", "09_docking_drug")
POSE = os.path.join(BASE, "docking_site", "VNN1_luteolin_out.pdbqt")
SDF  = os.path.join(BASE, "ligands", "luteolin.sdf")
REC  = os.path.join(BASE, "cocrystal", "7SLV.pdb")
OUT  = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_figS25_work")
os.makedirs(OUT, exist_ok=True)

# ---------- 1. 参考分子（含 3D 坐标 + 正确键级） ----------
ref = Chem.MolFromMolFile(SDF, removeHs=True, sanitize=True)
assert ref is not None, "SDF 读取失败"
ref = Chem.RemoveHs(ref)
print("参考分子:", Chem.rdMolDescriptors.CalcMolFormula(ref),
      "重原子", ref.GetNumAtoms())

# ---------- 2. 读出对接姿态 MODEL 1 ----------
AD2EL = {"A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N", "OA": "O", "OS": "O",
         "O": "O", "SA": "S", "S": "S", "H": "H", "HD": "H", "HS": "H", "P": "P",
         "F": "F", "Cl": "Cl", "CL": "Cl", "BR": "Br", "Br": "Br", "I": "I"}


def read_pose_pdbqt(path, model=1):
    """pdbqt 的 77-79 列是 AutoDock 原子类型(A/OA/HD)，需映射为元素"""
    atoms, cur, inmodel = [], None, False
    for line in open(path):
        if line.startswith("MODEL"):
            cur = int(line.split()[1]); inmodel = (cur == model); continue
        if line.startswith("ENDMDL"):
            if inmodel: break
            continue
        if inmodel and line.startswith(("ATOM", "HETATM")):
            ad = line[77:79].strip()
            el = AD2EL.get(ad, AD2EL.get(ad.capitalize(), None))
            if el is None:
                el = line[12:16].strip()[0]
            atoms.append((el, float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return atoms

pose_all = read_pose_pdbqt(POSE)
heavy = [a for a in pose_all if a[0] != "H"]
print(f"姿态原子: 总 {len(pose_all)}, 重原子 {len(heavy)}")
assert len(heavy) == ref.GetNumAtoms(), "重原子数不一致"

# ---------- 3. 建立 pose→ref 的图同构映射 ----------
pdb_lines = []
for i, (el, x, y, z) in enumerate(pose_all, 1):
    aname = f"{el}{i}" if len(el) == 1 else f"{el}{i}"
    pdb_lines.append(
        f"HETATM{i:5d} {aname:<4s} LIG A   1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}")
pdb_block = "\n".join(pdb_lines) + "\nEND\n"
pose_mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False)
pose_mol.UpdatePropertyCache(strict=False)
# 用距离判定连接
rdDetermine = None
try:
    from rdkit.Chem import rdDetermineBonds
    rdDetermineBonds.DetermineConnectivity(pose_mol)
except Exception as e:
    print("DetermineConnectivity 失败:", e)

_qp = Chem.AdjustQueryParameters()
_qp.makeBondsGeneric = True
query = Chem.AdjustQueryProperties(ref, _qp)
match = pose_mol.GetSubstructMatch(query)
print("子结构匹配原子数:", len(match), "/", ref.GetNumAtoms())
assert len(match) == ref.GetNumAtoms(), "图同构匹配失败"

# pose 原子索引（含 H 的原始序）→ ref 原子索引
pose2ref = {p: r for r, p in enumerate(match)}
# 重原子在 pose_all 中的序号
heavy_all_idx = [i for i, a in enumerate(pose_all) if a[0] != "H"]
pose_xyz = {pose2ref[p]: np.array(pose_all[p][1:]) for p in match}
ref_xyz = ref.GetConformer().GetPositions()
# 验证：Kabsch 后 RMSD
P = np.array([pose_xyz[r] for r in range(ref.GetNumAtoms())])
Q = ref_xyz
Pc, Qc = P - P.mean(0), Q - Q.mean(0)
U, S, Vt = np.linalg.svd(Pc.T @ Qc)
d = np.sign(np.linalg.det(U @ Vt))
R = U @ np.diag([1, 1, d]) @ Vt
rmsd = np.sqrt((((Pc @ R) - Qc) ** 2).sum(1).mean())
print(f"pose vs 晶体/SDF 构象 RMSD(去平移旋转) = {rmsd:.3f} Å")

# ---------- 4. 受体原子 ----------
rec_atoms = []
for line in open(REC):
    if line.startswith(("ATOM", "HETATM")):
        if line[21] != "A":  # 只要 A 链
            continue
        el = line[76:78].strip() or line[12:16].strip()[0]
        if el == "H":
            continue
        rec_atoms.append((int(line[22:26]), line[17:20].strip(), line[12:16].strip(),
                          el, np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
print("受体 A 链重原子:", len(rec_atoms))

TARGETS = [("GLU164", "hydrogen bond"), ("LYS157", "hydrogen bond"),
           ("MET218", "hydrogen bond"), ("TRP64", "hydrophobic"),
           ("PHE161", "hydrophobic"), ("PHE191", "hydrophobic"),
           ("TRP217", "pi-stacking"), ("PHE161", "pi-stacking")]

rows = []
for resn, itype in TARGETS:
    best = None
    for (resnum, resname, aname, el, xyz) in rec_atoms:
        if f"{resname}{resnum}" != resn:
            continue
        for r in range(ref.GetNumAtoms()):
            dmin = np.linalg.norm(pose_xyz[r] - xyz)
            if best is None or dmin < best[0]:
                best = (dmin, r, aname, el)
    rows.append(dict(residue=resn, interaction=itype, lig_atom_idx0=best[1],
                     lig_symbol=ref.GetAtomWithIdx(best[1]).GetSymbol(),
                     rec_atom=best[2], dist=round(best[0], 2)))
    print(f"{resn:8s} {itype:14s} → 配体原子 #{best[1]} ({ref.GetAtomWithIdx(best[1]).GetSymbol()}) "
          f"↔ {best[2]:4s}  {best[0]:.2f} Å")

import csv
with open(f"{OUT}/contact_atoms.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

# ---------- 5. 2D 坐标 ----------
mol2d = Chem.Mol(ref)
rdDepictor.SetPreferCoordGen(True)
rdDepictor.Compute2DCoords(mol2d)
conf = mol2d.GetConformer()
with open(f"{OUT}/atoms2d.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["idx0", "symbol", "x", "y", "nH", "arom"])
    for a in mol2d.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        w.writerow([a.GetIdx(), a.GetSymbol(), f"{p.x:.4f}", f"{p.y:.4f}",
                    a.GetTotalNumHs(), int(a.GetIsAromatic())])
with open(f"{OUT}/bonds.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["i", "j", "order", "arom"])
    for b in mol2d.GetBonds():
        bt = b.GetBondTypeAsDouble()
        w.writerow([b.GetBeginAtomIdx(), b.GetEndAtomIdx(),
                    int(bt) if bt in (1, 2, 3) else 1, int(b.GetIsAromatic())])
print("已写出", OUT)
print("SMILES:", Chem.MolToSmiles(mol2d))
