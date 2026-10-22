#!/usr/bin/env python3
"""
Step 4 补充 —— 定点对接（共结晶位点）+ 扩充配体
- UPP1: 用共结晶结构 3EUF (1.9Å), BAU 配体质心做 box 中心
- VNN1: 用共结晶结构 7SLV (2.13Å), 9TC 配体质心做 box 中心
- 其余 4 基因: 用 AlphaFold 结构 + 质心 box(缩小到 30Å)
- 配体: 10 个天然产物
"""
import os, subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/09_docking_drug"
VINA = f"{BASE}/bin/vina"
DOCK = f"{BASE}/docking_site"
os.makedirs(DOCK, exist_ok=True)

# 共结晶位点坐标（已从 PDB 提取）
COCRYSTAL = {
    "UPP1": {"pdb": "cocrystal/3EUF.pdb", "chain": "A", "center": (34.725, 55.011, -32.831), "box": 22.0},
    "VNN1": {"pdb": "cocrystal/7SLV.pdb", "chain": "A", "center": (18.043, -28.315, -26.014), "box": 22.0},
}
# 其余基因用 AlphaFold（无共结晶配体）
AF2 = {
    "CDK5RAP2": "receptors/CDK5RAP2_Q96SN8.pdb",
    "CD177": "receptors/CD177_Q8N6Q3.pdb",
    "GRB10": "receptors/GRB10_Q13322.pdb",
    "EXOSC4": "receptors/EXOSC4_Q9NPD3.pdb",
}

# 10 个天然产物
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
    """从共结晶 PDB 提取指定链的蛋白原子(去 HETATM/水)"""
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
    mol_setups = prep.prepare(mol)
    if not mol_setups: return False
    with open(out_pdbqt, "w") as f:
        f.write(prep.write_pdbqt_string())
    return True

def get_centroid(pdbqt_path, box=30.0):
    xs, ys, zs = [], [], []
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
    return sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs), box

def dock(gene, rec_pdbqt, cx, cy, cz, box):
    """对 10 个配体做定点对接, 返回结果列表"""
    res = []
    for lname, sdf in LIGANDS.items():
        lig_pdbqt = f"{DOCK}/{gene}_{lname}_lig.pdbqt"
        if not prep_ligand_pdbqt(f"{BASE}/{sdf}", lig_pdbqt):
            print(f"  {lname} 配体准备失败"); continue
        out_pdbqt = f"{DOCK}/{gene}_{lname}_out.pdbqt"
        cmd = [VINA, "--receptor", rec_pdbqt, "--ligand", lig_pdbqt,
               "--center_x", f"{cx:.2f}", "--center_y", f"{cy:.2f}", "--center_z", f"{cz:.2f}",
               "--size_x", f"{box:.2f}", "--size_y", f"{box:.2f}", "--size_z", f"{box:.2f}",
               "--exhaustiveness", "16", "--num_modes", "5", "--out", out_pdbqt]
        r = subprocess.run(cmd, capture_output=True, text=True)
        aff = None
        for line in r.stdout.split("\n"):
            parts = line.split()
            if len(parts) >= 4:
                try:
                    if parts[0].isdigit() and float(parts[1]) < 0:
                        aff = float(parts[1]); break
                except: pass
        if aff is None and os.path.exists(out_pdbqt):
            with open(out_pdbqt) as f:
                for line in f:
                    if "VINA RESULT" in line:
                        aff = float(line.split()[3]); break
        print(f"  {gene} - {lname}: {aff if aff is not None else 'N/A'} kcal/mol")
        res.append({"gene": gene, "ligand": lname, "affinity": aff,
                    "site": "cocrystal" if gene in COCRYSTAL else "centroid"})
    return res

all_res = []
# 1) 共结晶定点对接（UPP1 + VNN1）
for gene, cfg in COCRYSTAL.items():
    print(f"\n=== {gene} (共结晶定点对接, box={cfg['box']}Å) ===")
    clean_pdb = f"{DOCK}/{gene}_protein.pdb"
    if not extract_protein_chain(f"{BASE}/{cfg['pdb']}", cfg["chain"], clean_pdb):
        print("  蛋白链提取失败"); continue
    rec_pdbqt = f"{DOCK}/{gene}_receptor.pdbqt"
    if not prep_receptor_pdbqt(clean_pdb, rec_pdbqt):
        print("  受体准备失败"); continue
    cx, cy, cz = cfg["center"]
    all_res += dock(gene, rec_pdbqt, cx, cy, cz, cfg["box"])

# 2) 其余 4 基因质心对接（box 缩小到 30Å）
for gene, pdb in AF2.items():
    print(f"\n=== {gene} (质心对接, box=30Å) ===")
    rec_pdbqt = f"{DOCK}/{gene}_receptor.pdbqt"
    if not prep_receptor_pdbqt(f"{BASE}/{pdb}", rec_pdbqt):
        print("  受体准备失败"); continue
    cx, cy, cz, _ = get_centroid(rec_pdbqt)
    all_res += dock(gene, rec_pdbqt, cx, cy, cz, 30.0)

import pandas as pd
df = pd.DataFrame(all_res)
df.to_csv(f"{BASE}/docking_site_results.csv", index=False)
print("\n===== 定点对接结果汇总 =====")
piv = df.pivot_table(index="gene", columns="ligand", values="affinity")
print(piv.to_string())
print("\n完成")
