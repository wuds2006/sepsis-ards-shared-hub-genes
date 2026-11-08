#!/usr/bin/env python3
"""导出 2D 骨架（Kekulé 键级）+ 接触锚点的 2D 坐标，供 R 绘制 Figure S25A。
输出 atoms2d.csv / bonds2d.csv / contacts2d.csv
"""
import os, csv, json
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor

BASE = os.environ.get("DOCK_ROOT", "09_docking_drug")
SDF = os.path.join(BASE, "ligands", "luteolin.sdf")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_figS25_work")

mol = Chem.RemoveHs(Chem.MolFromMolFile(SDF, removeHs=True))
rdDepictor.SetPreferCoordGen(True)
rdDepictor.Compute2DCoords(mol)
Chem.Kekulize(mol, clearAromaticFlags=False)
conf = mol.GetConformer()

atoms = []
for a in mol.GetAtoms():
    p = conf.GetAtomPosition(a.GetIdx())
    atoms.append(dict(idx0=a.GetIdx(), symbol=a.GetSymbol(),
                      x=round(p.x, 4), y=round(p.y, 4),
                      nH=a.GetTotalNumHs(), aromatic=int(a.GetIsAromatic()),
                      degree=a.GetDegree()))
with open(f"{OUT}/atoms2d.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(atoms[0].keys())); w.writeheader(); w.writerows(atoms)

bonds = []
for b in mol.GetBonds():
    bonds.append(dict(i=b.GetBeginAtomIdx(), j=b.GetEndAtomIdx(),
                      order=int(b.GetBondTypeAsDouble()), arom=int(b.GetIsAromatic())))
with open(f"{OUT}/bonds2d.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(bonds[0].keys())); w.writeheader(); w.writerows(bonds)

# --- 环信息（用于把双键的第二条线画在环内侧） ---
rings = []
for rid, ring in enumerate(mol.GetRingInfo().AtomRings()):
    for idx in ring:
        rings.append(dict(ring=rid, idx0=idx))
with open(f"{OUT}/rings2d.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["ring", "idx0"]); w.writeheader(); w.writerows(rings)

# --- 接触锚点的 2D 坐标 ---
anchors = json.load(open(f"{OUT}/contact_anchors.json"))
xy = {a["idx0"]: np.array([a["x"], a["y"]]) for a in atoms}
rows = []
for a in anchors:
    if a["kind"] == "ring":
        pts = np.array([xy[i] for i in a["lig_atoms"]])
        c = pts.mean(0)
        rows.append(dict(residue=a["residue"], interaction=a["interaction"],
                         dist=a["dist"], x=round(c[0], 4), y=round(c[1], 4), kind="ring"))
    else:
        i = a["lig_atoms"][0]
        rows.append(dict(residue=a["residue"], interaction=a["interaction"],
                         dist=a["dist"], x=round(xy[i][0], 4), y=round(xy[i][1], 4),
                         kind="atom", lig_idx0=i))
with open(f"{OUT}/contacts2d.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=sorted({k for r in rows for k in r}))
    w.writeheader()
    for r in rows: w.writerow(r)

print("原子 2D 坐标：")
for a in atoms:
    print(f"  #{a['idx0']:2d} {a['symbol']:2s} ({a['x']:6.2f},{a['y']:6.2f}) deg={a['degree']} nH={a['nH']}")
print("\n接触锚点：")
for r in rows:
    print(f"  {r['residue']:7s} {r['interaction']:13s} {r['dist']:5.2f}  @ ({r['x']:6.2f},{r['y']:6.2f})  {r['kind']}")
print("\nSMILES:", Chem.MolToSmiles(mol))
