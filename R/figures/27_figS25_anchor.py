#!/usr/bin/env python3
"""按作用类型重建 VNN1–luteolin 的接触原子锚定，并核对 PLIP CSV 距离。

规则（与 PLIP 定义一致）：
  hydrogen bond : 配体 O/N  ↔ 受体 N/O/S
  hydrophobic   : 配体 C     ↔ 受体 C/S
  pi-stacking   : 配体芳环   ↔ 受体芳环（质心—质心距离）
输出 contact_anchors.json 供 R 绘图使用。
"""
import os, json, itertools
import numpy as np
from rdkit import Chem

BASE = os.environ.get("DOCK_ROOT", "09_docking_drug")
POSE = os.path.join(BASE, "docking_site", "VNN1_luteolin_out.pdbqt")
SDF = os.path.join(BASE, "ligands", "luteolin.sdf")
REC = os.path.join(BASE, "cocrystal", "7SLV.pdb")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_figS25_work")
os.makedirs(OUT, exist_ok=True)

AD2EL = {"A": "C", "C": "C", "N": "N", "NA": "N", "OA": "O", "OS": "O", "O": "O",
         "SA": "S", "S": "S", "H": "H", "HD": "H", "HS": "H"}


def read_pose(path, model=1):
    out, cur, inm = [], None, False
    for line in open(path):
        if line.startswith("MODEL"):
            cur = int(line.split()[1]); inm = (cur == model); continue
        if line.startswith("ENDMDL"):
            if inm: break
            continue
        if inm and line.startswith(("ATOM", "HETATM")):
            ad = line[77:79].strip()
            el = AD2EL.get(ad, line[12:16].strip()[0])
            out.append((el, np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    return out


ref = Chem.RemoveHs(Chem.MolFromMolFile(SDF, removeHs=True))
pose_all = read_pose(POSE)
heavy = [a for a in pose_all if a[0] != "H"]

# --- pose -> ref 图同构映射（同 26_ 脚本） ---
lines = []
for i, (el, xyz) in enumerate(pose_all, 1):
    an = f"{el}{i}"
    lines.append(f"HETATM{i:5d} {an:<4s} LIG A   1    "
                 f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {el:>2s}")
pm = Chem.MolFromPDBBlock("\n".join(lines) + "\nEND\n", removeHs=False, sanitize=False)
pm.UpdatePropertyCache(strict=False)
from rdkit.Chem import rdDetermineBonds
rdDetermineBonds.DetermineConnectivity(pm)
qp = Chem.AdjustQueryParameters(); qp.makeBondsGeneric = True
match = pm.GetSubstructMatch(Chem.AdjustQueryProperties(ref, qp))
pose2ref = {p: r for r, p in enumerate(match)}
lig_xyz = {pose2ref[p]: pose_all[p][1] for p in match}          # ref 索引 -> 3D 坐标
lig_el = {i: ref.GetAtomWithIdx(i).GetSymbol() for i in range(ref.GetNumAtoms())}
lig_arom = {i: ref.GetAtomWithIdx(i).GetIsAromatic() for i in range(ref.GetNumAtoms())}

# --- 受体 A 链 ---
rec = {}
for line in open(REC):
    if line.startswith(("ATOM", "HETATM")) and line[21] == "A":
        el = line[76:78].strip() or line[12:16].strip()[0]
        if el == "H":
            continue
        key = f"{line[17:20].strip()}{int(line[22:26])}"
        rec.setdefault(key, []).append(
            (line[12:16].strip(), el,
             np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
print("受体残基原子组:", len(rec))

# --- 配体芳环（SSSR 中的芳香环） ---
ri = ref.GetRingInfo()
arom_rings = [list(r) for r in ri.AtomRings()
              if all(ref.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
print("配体芳香环:", arom_rings)

TARGETS = [("GLU164", "hydrogen bond", 2.87), ("LYS157", "hydrogen bond", 2.90),
           ("MET218", "hydrogen bond", 3.42), ("TRP64", "hydrophobic", 3.12),
           ("PHE161", "hydrophobic", 3.72), ("PHE191", "hydrophobic", 3.38),
           ("TRP217", "pi-stacking", 4.57), ("PHE161", "pi-stacking", 3.87)]

ALLOWED = {"hydrogen bond": ({"O", "N"}, {"N", "O", "S"}),
           "hydrophobic":   ({"C"}, {"C", "S"}),
           "pi-stacking":   ({"C"}, {"C", "N"})}

anchors, detail = [], []
for resn, itype, rep in TARGETS:
    la_ok, ra_ok = ALLOWED[itype]
    if itype == "pi-stacking":
        # 配体环质心 × 受体环质心（受体六元芳环：PHE/TRP 侧链）
        rec_rings, ring_names = [], []
        names = {"PHE": [["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]],
                 "TRP": [["CG", "CD1", "NE1", "CE2", "CD2"],
                         ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
                         ["CG", "CD2", "CE2", "CZ3", "CH2", "CZ2"]]}
        rn = resn[:3]
        for rnames in names.get(rn, []):
            pts = [xyz for (an, el, xyz) in rec[resn] if an in rnames]
            if len(pts) == len(rnames):
                rec_rings.append(np.mean(pts, axis=0)); ring_names.append(rnames)
        best = None
        for k, lidx in enumerate(arom_rings):
            lc = np.mean([lig_xyz[i] for i in lidx], axis=0)
            for j, rc in enumerate(rec_rings):
                d = np.linalg.norm(lc - rc)
                if best is None or abs(d - rep) < abs(best[0] - rep):
                    best = (d, lidx, lc, rc, j)
        d, lidx, lc, rc, j = best
        anchors.append(dict(residue=resn, interaction=itype, reported=rep,
                            kind="ring", lig_atoms=lidx,
                            lig_xy=None, dist=round(float(d), 2),
                            rec_atoms=ring_names[j]))
        print(f"{resn:8s} {itype:13s} 配体环{k} 质心距离 {d:.2f} Å (PLIP {rep})")
        continue

    best = None
    for (an, el, xyz) in rec[resn]:
        if el not in ra_ok:
            continue
        for li in range(ref.GetNumAtoms()):
            if lig_el[li] not in la_ok:
                continue
            d = float(np.linalg.norm(lig_xyz[li] - xyz))
            if best is None or abs(d - rep) < abs(best[0] - rep):
                best = (d, li, an, el)
    d, li, an, el = best
    anchors.append(dict(residue=resn, interaction=itype, reported=rep, kind="atom",
                        lig_atoms=[li], dist=round(d, 2), rec_atoms=[an],
                        rec_element=el, lig_element=lig_el[li]))
    print(f"{resn:8s} {itype:13s} 配体 #{li:2d} {lig_el[li]:2s} ↔ {an:4s} {el:2s}  "
          f"{d:.2f} Å (PLIP {rep})  偏差 {abs(d-rep):.2f}")

    # 记录该残基全部近距原子对，供核对
    for (an, el, xyz) in rec[resn]:
        for li in range(ref.GetNumAtoms()):
            if lig_el[li] not in la_ok or el not in ra_ok:
                continue
            d = float(np.linalg.norm(lig_xyz[li] - xyz))
            if d < 4.2:
                detail.append(dict(residue=resn, interaction=itype, lig=li,
                                   lig_el=lig_el[li], rec_atom=an, rec_el=el, dist=round(d, 2)))

with open(f"{OUT}/contact_anchors.json", "w") as f:
    json.dump(anchors, f, indent=1, ensure_ascii=False)
import csv
with open(f"{OUT}/candidate_pairs.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(detail[0].keys())); w.writeheader(); w.writerows(detail)
print("\n写出 contact_anchors.json / candidate_pairs.csv")
