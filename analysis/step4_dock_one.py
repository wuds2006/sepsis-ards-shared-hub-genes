#!/usr/bin/env python3
"""定点对接：按基因独立运行（避免内存累积）。用法: python step4_dock_one.py <GENE>"""
import sys, os, subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/09_docking_drug"
VINA = f"{BASE}/bin/vina"
DOCK = f"{BASE}/docking_site"
os.makedirs(DOCK, exist_ok=True)

GENE = sys.argv[1] if len(sys.argv) > 1 else "UPP1"

COCRYSTAL = {
    "UPP1":  {"pdb": "cocrystal/3EUF.pdb", "chain": "A", "center": (34.725, 55.011, -32.831), "box": 22.0},
    "VNN1":  {"pdb": "cocrystal/7SLV.pdb", "chain": "A", "center": (18.043, -28.315, -26.014), "box": 22.0},
}
AF2 = {
    "CDK5RAP2": "receptors/CDK5RAP2_Q96SN8.pdb",
    "CD177": "receptors/CD177_Q8N6Q3.pdb",
    "GRB10": "receptors/GRB10_Q13322.pdb",
    "EXOSC4": "receptors/EXOSC4_Q9NPD3.pdb",
}
LIGANDS = ["EGCG","quercetin","tetrahydroalstonine","resveratrol","curcumin",
           "baicalein","luteolin","kaempferol","berberine","emodin"]

def auto_type(atom):
    sym = atom.GetSymbol(); arom = atom.GetIsAromatic()
    if sym == "C": return "A" if arom else "C"
    if sym == "N":
        hasH = any(n.GetSymbol()=="H" for n in atom.GetNeighbors())
        return "N" if hasH else "NA"
    if sym == "O": return "OA"
    if sym == "H":
        nb = atom.GetNeighbors()[0] if atom.GetNeighbors() else None
        return "HD" if (nb and nb.GetSymbol() in ("O","N")) else "H"
    if sym == "S": return "SA" if arom else "S"
    return sym

def write_receptor(mol, out):
    conf = mol.GetConformer()
    with open(out, "w") as f:
        f.write("REMARK receptor\n")
        for i, atom in enumerate(mol.GetAtoms(), 1):
            pos = conf.GetAtomPosition(atom.GetIdx())
            f.write(f"HETATM{i:5d}  {atom.GetSymbol():<3s} UNL A   1    "
                    f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}  1.00  0.00           {auto_type(atom)}\n")

# 1) 受体
rec_pdbqt = f"{DOCK}/{GENE}_receptor.pdbqt"
if GENE in COCRYSTAL:
    cfg = COCRYSTAL[GENE]
    clean = f"{DOCK}/{GENE}_protein.pdb"
    with open(f"{BASE}/{cfg['pdb']}") as f, open(clean, "w") as w:
        for line in f:
            if line.startswith("ATOM") and line[21:22].strip() == cfg["chain"]:
                w.write(line)
            elif line.startswith("TER"):
                w.write("END\n"); break
    mol = Chem.MolFromPDBFile(clean, removeHs=True)
    mol = AllChem.AddHs(mol, addCoords=True)
    write_receptor(mol, rec_pdbqt)
    cx, cy, cz, box = cfg["center"][0], cfg["center"][1], cfg["center"][2], cfg["box"]
    site = "cocrystal"
else:
    mol = Chem.MolFromPDBFile(f"{BASE}/{AF2[GENE]}", removeHs=True)
    mol = AllChem.AddHs(mol, addCoords=True)
    write_receptor(mol, rec_pdbqt)
    # 质心
    conf = mol.GetConformer()
    xs = [conf.GetAtomPosition(a.GetIdx()).x for a in mol.GetAtoms()]
    ys = [conf.GetAtomPosition(a.GetIdx()).y for a in mol.GetAtoms()]
    zs = [conf.GetAtomPosition(a.GetIdx()).z for a in mol.GetAtoms()]
    cx, cy, cz = sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)
    box = 40.0
    site = "centroid"

# 2) 逐个配体对接
results = []
for lname in LIGANDS:
    sdf = f"{BASE}/ligands/{lname}.sdf"
    lig_pdbqt = f"{DOCK}/{GENE}_{lname}_lig.pdbqt"
    lig = Chem.MolFromMolFile(sdf, removeHs=False)
    if lig is None: print(f"  {lname} 读失败"); continue
    lig = AllChem.AddHs(lig, addCoords=True)
    prep = MoleculePreparation()
    if not prep.prepare(lig): print(f"  {lname} 准备失败"); continue
    open(lig_pdbqt, "w").write(prep.write_pdbqt_string())
    out_pdbqt = f"{DOCK}/{GENE}_{lname}_out.pdbqt"
    cmd = [VINA, "--receptor", rec_pdbqt, "--ligand", lig_pdbqt,
           "--center_x", f"{cx:.2f}", "--center_y", f"{cy:.2f}", "--center_z", f"{cz:.2f}",
           "--size_x", f"{box:.2f}", "--size_y", f"{box:.2f}", "--size_z", f"{box:.2f}",
           "--exhaustiveness", "16", "--num_modes", "5", "--out", out_pdbqt]
    r = subprocess.run(cmd, capture_output=True, text=True)
    aff = None
    for line in r.stdout.split("\n"):
        parts = line.split()
        if len(parts) >= 4 and parts[0].isdigit():
            try:
                if float(parts[1]) < 0: aff = float(parts[1]); break
            except: pass
    print(f"  {GENE} - {lname}: {aff if aff is not None else 'N/A'} kcal/mol")
    results.append({"gene": GENE, "ligand": lname, "affinity": aff, "site": site})

import pandas as pd
pd.DataFrame(results).to_csv(f"{DOCK}/{GENE}_results.csv", index=False)
print(f"===== {GENE} 定点对接完成 ({site}) =====")
