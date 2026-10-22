#!/usr/bin/env python3
"""
Step 4 —— 分子对接: 6 个 Hub 基因蛋白 x 3 个天然产物 (AutoDock Vina)
"""
import os, subprocess, glob
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTReceptor

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/09_docking_drug"
VINA = f"{BASE}/bin/vina"

receptors = {
    "CDK5RAP2": "receptors/CDK5RAP2_Q96SN8.pdb",
    "CD177": "receptors/CD177_Q8N6Q3.pdb",
    "GRB10": "receptors/GRB10_Q13322.pdb",
    "UPP1": "receptors/UPP1_Q16831.pdb",
    "EXOSC4": "receptors/EXOSC4_Q9NPD3.pdb",
    "VNN1": "receptors/VNN1_O95497.pdb",
}
ligands = {
    "EGCG": "ligands/EGCG.sdf",
    "quercetin": "ligands/quercetin.sdf",
    "tetrahydroalstonine": "ligands/tetrahydroalstonine.sdf",
}

def prep_receptor(pdb_path, out_pdbqt):
    """受体准备: rdkit 去水加氢, 手动构造 PDBQT"""
    mol = Chem.MolFromPDBFile(pdb_path, removeHs=True)
    if mol is None:
        return False
    mol = AllChem.AddHs(mol, addCoords=True)
    conf = mol.GetConformer()
    # AutoDock 原子类型映射
    def auto_type(atom):
        sym = atom.GetSymbol()
        arom = atom.GetIsAromatic()
        if sym == "C":
            return "A" if arom else "C"
        if sym == "N":
            # 判断是否带 H(供体)
            hasH = any(n.GetSymbol() == "H" for n in atom.GetNeighbors())
            return "N" if hasH else "NA"
        if sym == "O":
            return "OA"
        if sym == "H":
            # 判断是否供体 H(连在 O/N)
            nb = atom.GetNeighbors()[0] if atom.GetNeighbors() else None
            if nb and nb.GetSymbol() in ("O", "N"):
                return "HD"
            return "H"
        if sym == "S":
            return "SA" if arom else "S"
        return sym
    with open(out_pdbqt, "w") as f:
        f.write("REMARK  AutoDock Vina receptor (rdkit-prepared)\n")
        for i, atom in enumerate(mol.GetAtoms(), 1):
            pos = conf.GetAtomPosition(atom.GetIdx())
            atype = auto_type(atom)
            f.write(f"HETATM{i:5d}  {atom.GetSymbol():<3s} UNL A   1    "
                    f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}"
                    f"  1.00  0.00           {atype}\n")
    return True

def prep_ligand(sdf_path, out_pdbqt):
    """配体准备: meeko"""
    mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
    if mol is None:
        return False
    mol = AllChem.AddHs(mol, addCoords=True)
    prep = MoleculePreparation()
    mol_setups = prep.prepare(mol)
    if not mol_setups:
        return False
    pdbqt_string = prep.write_pdbqt_string()
    with open(out_pdbqt, "w") as f:
        f.write(pdbqt_string)
    return True

def get_center_size(pdbqt_path, box=40.0):
    """从受体 PDBQT 计算质心 + 固定大小 box(避免超大蛋白 OOM)"""
    xs, ys, zs = [], [], []
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
    cx = sum(xs)/len(xs); cy = sum(ys)/len(ys); cz = sum(zs)/len(zs)
    return cx, cy, cz, box, box, box

results = []
for gname, pdb in receptors.items():
    print(f"\n=== {gname} ===")
    rec_pdbqt = f"{BASE}/docking/{gname}_receptor.pdbqt"
    os.makedirs(f"{BASE}/docking", exist_ok=True)
    if not prep_receptor(f"{BASE}/{pdb}", rec_pdbqt):
        print("  受体准备失败, 跳过"); continue
    cx, cy, cz, sx, sy, sz = get_center_size(rec_pdbqt)
    for lname, sdf in ligands.items():
        lig_pdbqt = f"{BASE}/docking/{gname}_{lname}_ligand.pdbqt"
        if not prep_ligand(f"{BASE}/{sdf}", lig_pdbqt):
            print(f"  {lname} 配体准备失败"); continue
        out_pdbqt = f"{BASE}/docking/{gname}_{lname}_out.pdbqt"
        cmd = [VINA, "--receptor", rec_pdbqt, "--ligand", lig_pdbqt,
               "--center_x", f"{cx:.2f}", "--center_y", f"{cy:.2f}", "--center_z", f"{cz:.2f}",
               "--size_x", f"{sx:.2f}", "--size_y", f"{sy:.2f}", "--size_z", f"{sz:.2f}",
               "--exhaustiveness", "16", "--num_modes", "5", "--out", out_pdbqt]
        r = subprocess.run(cmd, capture_output=True, text=True)
        # 解析结合能
        aff = None
        for line in r.stdout.split("\n"):
            if line.strip().startswith(("1", "2", "3", "4", "5")) and len(line.split()) >= 4:
                parts = line.split()
                try:
                    float(parts[0]); float(parts[1])
                    aff = float(parts[1]); break
                except: pass
        if aff is None and os.path.exists(out_pdbqt):
            # 从 out 文件解析
            with open(out_pdbqt) as f:
                for line in f:
                    if "VINA RESULT" in line:
                        aff = float(line.split()[3]); break
        print(f"  {gname} - {lname}: {aff if aff is not None else 'N/A'} kcal/mol")
        results.append({"gene": gname, "ligand": lname, "affinity": aff})

import pandas as pd
df = pd.DataFrame(results)
df.to_csv(f"{BASE}/docking_results.csv", index=False)
print("\n===== 对接结果 =====")
print(df.to_string(index=False))
print("\n完成")
