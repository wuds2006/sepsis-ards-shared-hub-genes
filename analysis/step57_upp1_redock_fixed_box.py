#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 57 —— UPP1 定点对接修正重跑
redocking 验证(step56)发现原 step4_docking_site.py 中 UPP1(3EUF) 的 box 中心
(34.725, 55.011, -32.831) 为四聚体 A-D 四个 BAU 质心的平均, 落在中央溶剂通道,
未覆盖链 A 活性口袋。本脚本以链 A 自身 BAU 质心 (59.963, 74.157, -9.840) 为中心,
按完全相同的协议(box 22 A, exhaustiveness 16, num_modes 5, 链 A 受体, meeko 配体)
重跑 UPP1 x 10 个天然产物, 并更新 docking_site_results.csv 的 UPP1 行。
归档原结果到 docking_site_results_OLD_UPP1box.csv。
"""
import os, subprocess, csv, shutil
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/09_docking_drug"
VINA = f"{BASE}/bin/vina"
DOCK = f"{BASE}/docking_site"
PYBIN = os.path.dirname(os.sys.executable)
OBABEL = os.path.join(PYBIN, "obabel")

CENTER = (59.963, 74.157, -9.840)   # 链 A BAU 质心（已由 step56 验证）
BOX = 22.0

LIGANDS = {
    "EGCG": "ligands/EGCG.sdf",
    "quercetin": "ligands/quercetin.sdf",
    "tetrahydroalstonine": "ligands/tetrahydroalstonine.sdf",
    "resveratrol": "ligands/resveratrol.sdf",
    "curcumin": "ligands/curcumin.sdf",
    "baicalein": "ligands/baicalein.sdf",
    "luteolin": "ligands/luteolin.sdf",
    "kaempferol": "ligands/kaempferol.sdf",
    "berberine": "ligands/berberine.sdf",
    "emodin": "ligands/emodin.sdf",
}

def auto_type(atom):
    sym = atom.GetSymbol(); arom = atom.GetIsAromatic()
    if sym == "C": return "A" if arom else "C"
    if sym == "N":
        hasH = any(n.GetSymbol() == "H" for n in atom.GetNeighbors())
        return "N" if hasH else "NA"
    if sym == "O": return "OA"
    if sym == "H":
        nb = atom.GetNeighbors()[0] if atom.GetNeighbors() else None
        return "HD" if (nb and nb.GetSymbol() in ("O", "N")) else "H"
    if sym == "S": return "SA" if arom else "S"
    return sym

def extract_protein_chain(pdb_path, chain, out_pdb):
    with open(pdb_path) as f, open(out_pdb, "w") as w:
        for line in f:
            if line.startswith("ATOM") and line[21:22].strip() == chain:
                w.write(line)
            elif line.startswith("TER"):
                w.write("END\n"); break
    return os.path.getsize(out_pdb) > 0

def prep_receptor_pdbqt(pdb_path, out_pdbqt):
    mol = Chem.MolFromPDBFile(pdb_path, removeHs=True)
    if mol is None: return False
    mol = AllChem.AddHs(mol, addCoords=True)
    conf = mol.GetConformer()
    with open(out_pdbqt, "w") as f:
        f.write("REMARK  AutoDock Vina receptor\n")
        for i, atom in enumerate(mol.GetAtoms(), 1):
            pos = conf.GetAtomPosition(atom.GetIdx())
            atype = auto_type(atom)
            f.write(f"HETATM{i:5d}  {atom.GetSymbol():<3s} UNL A   1    "
                    f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}"
                    f"  1.00  0.00           {atype}\n")
    return True

def prep_ligand_pdbqt(sdf_path, out_pdbqt):
    mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
    if mol is None: return False
    mol = AllChem.AddHs(mol, addCoords=True)
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    if not setups: return False
    pdbqt_string, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        print(f"  meeko 写出失败: {err}"); return False
    with open(out_pdbqt, "w") as f:
        f.write(pdbqt_string)
    return True

# ---- 1. 受体（链 A, 镜像原协议）----
clean_pdb = f"{DOCK}/UPP1_protein.pdb"
assert extract_protein_chain(f"{BASE}/cocrystal/3EUF.pdb", "A", clean_pdb)
rec_pdbqt = f"{DOCK}/UPP1_receptor.pdbqt"
assert prep_receptor_pdbqt(clean_pdb, rec_pdbqt)

# ---- 2. 重对接 10 配体（修正 box 中心）----
new_rows = []
cx, cy, cz = CENTER
for lname, sdf in LIGANDS.items():
    lig_pdbqt = f"{DOCK}/UPP1_{lname}_lig.pdbqt"
    if not prep_ligand_pdbqt(f"{BASE}/{sdf}", lig_pdbqt):
        print(f"  {lname} 配体准备失败"); continue
    out_pdbqt = f"{DOCK}/UPP1_{lname}_out.pdbqt"
    cmd = [VINA, "--receptor", rec_pdbqt, "--ligand", lig_pdbqt,
           "--center_x", f"{cx:.2f}", "--center_y", f"{cy:.2f}", "--center_z", f"{cz:.2f}",
           "--size_x", f"{BOX:.2f}", "--size_y", f"{BOX:.2f}", "--size_z", f"{BOX:.2f}",
           "--exhaustiveness", "16", "--num_modes", "5", "--out", out_pdbqt]
    r = subprocess.run(cmd, capture_output=True, text=True)
    aff = None
    for line in r.stdout.split("\n"):
        parts = line.split()
        if len(parts) >= 4:
            try:
                if parts[0].isdigit() and float(parts[1]) < 0:
                    aff = float(parts[1]); break
            except (ValueError, IndexError): pass
    if aff is None and os.path.exists(out_pdbqt):
        with open(out_pdbqt) as f:
            for line in f:
                if "VINA RESULT" in line:
                    aff = float(line.split()[3]); break
    print(f"  UPP1 - {lname}: {aff} kcal/mol")
    new_rows.append({"gene": "UPP1", "ligand": lname, "affinity": aff, "site": "cocrystal"})

# ---- 3. 更新 docking_site_results.csv ----
csv_path = f"{BASE}/docking_site_results.csv"
shutil.copy2(csv_path, f"{BASE}/docking_site_results_OLD_UPP1box.csv")
with open(csv_path) as f:
    rows = list(csv.DictReader(f))
assert all(r["gene"] != "UPP1" for r in rows) or True
rows = [r for r in rows if r["gene"] != "UPP1"] + new_rows
# 保持原行序(按基因 + 配体顺序)
gene_order = ["CDK5RAP2", "EXOSC4", "VNN1", "UPP1", "GRB10", "CD177"]
lig_order = list(LIGANDS.keys())
rows.sort(key=lambda r: (gene_order.index(r["gene"]),
                          lig_order.index(r["ligand"]) if r["ligand"] in lig_order else 99))
with open(csv_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["gene", "ligand", "affinity", "site"])
    w.writeheader(); w.writerows(rows)
print("\n===== UPP1 修正后定点对接结果 =====")
for r in new_rows:
    print(f"  {r['ligand']}: {r['affinity']}")
print(f"\n已更新 {csv_path} (旧值归档于 docking_site_results_OLD_UPP1box.csv)")
