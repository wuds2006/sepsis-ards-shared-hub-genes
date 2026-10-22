#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 56 —— Redocking 验证（审稿意见 B7：分子对接协议可靠性）
- UPP1: 3EUF (1.9 A) 共晶配体 BAU 重对接, 计算重原子 RMSD
- VNN1: 7SLV (2.13 A) 共晶配体 9TC 重对接, 计算重原子 RMSD
- 协议完全镜像 step4_docking_site.py：链 A 受体（rdkit 去水加氢 + 手动 PDBQT）、
  共晶配体质心 box 22 A、exhaustiveness 16、num_modes 5、meeko 配体准备
- 判定：最佳姿态重原子 RMSD <= 2.0 A 视为协议可复现共晶构象
输出：09_docking_drug/redocking_results.csv
"""
import os, subprocess, csv, sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

PYBIN  = os.path.dirname(sys.executable)
OBABEL = os.path.join(PYBIN, "obabel")
BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/09_docking_drug"
VINA = f"{BASE}/bin/vina"
OUTDIR = f"{BASE}/redocking"
os.makedirs(OUTDIR, exist_ok=True)

CASES = {
    # 注意：3EUF 为四聚体(A-D 各含一个 BAU)。原 step4_docking_site.py 存档中心
    # (34.725, 55.011, -32.831) 为四个 BAU 质心的平均, 落在四聚体中央溶剂通道;
    # 正确定点 box 应取链 A 自身 BAU 质心 (59.963, 74.157, -9.840)。
    "UPP1_3EUF": {"pdb": f"{BASE}/cocrystal/3EUF.pdb", "chain": "A", "resname": "BAU",
                  "center": (59.963, 74.157, -9.840), "box": 22.0},
    "VNN1_7SLV": {"pdb": f"{BASE}/cocrystal/7SLV.pdb", "chain": "A", "resname": "9TC",
                  "center": (18.043, -28.315, -26.014), "box": 22.0},
}

# ---------- 与 step4_docking_site.py 完全一致的受体准备 ----------
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

# ---------- 共晶配体提取 ----------
def extract_cocrystal_ligand(pdb_path, chain, resname, out_pdb):
    """提取指定链上的共晶配体重原子 + CONECT 键信息"""
    serials = []
    with open(pdb_path) as f, open(out_pdb, "w") as w:
        for line in f:
            if line.startswith("HETATM") and line[17:20].strip() == resname \
               and line[21:22].strip() == chain:
                serials.append(int(line[6:11]))
                w.write(line)
    # 追加涉及这些原子的 CONECT 记录（键级用于 obabel 感知）
    want = set(serials)
    with open(pdb_path) as f, open(out_pdb, "a") as w:
        for line in f:
            if line.startswith("CONECT"):
                refs = [int(line[6:11])] if line[6:11].strip() else []
                for k in range(11, 71, 5):
                    if k + 5 <= len(line) and line[k:k+5].strip():
                        refs.append(int(line[k:k+5]))
                if refs and refs[0] in want:
                    w.write(line)
    return len(serials)

# ---------- pose 解析与 RMSD ----------
def _topo_heavy(mol):
    """去氢 + 全部键压平为单键的拓扑副本（消除芳香/键级表示差异, 保留构象）"""
    m = Chem.RemoveHs(Chem.Mol(mol))
    rw = Chem.RWMol(m)
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in rw.GetAtoms():
        a.SetIsAromatic(False)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
    out = rw.GetMol()
    Chem.SanitizeMol(out, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SYMMRINGS)
    return out

def pose_mols_from_pdbqt(path):
    """meeko 官方姿态重建：返回带 5 个 conformer 的 RDKit mol"""
    from meeko import PDBQTMolecule, RDKitMolCreate
    pq = PDBQTMolecule.from_file(path, skip_typing=True)
    mols = RDKitMolCreate.from_pdbqt_mol(pq)
    if isinstance(mols, tuple):
        mols = mols[0]
    assert mols and mols[0] is not None, "姿态重建失败"
    return mols[0]

def rmsd_inplace(pose_mol, crystal_mol):
    """原位(不对齐)、对称校正的重原子 RMSD, 取所有对称映射的最小值"""
    pt, ct = _topo_heavy(pose_mol), _topo_heavy(crystal_mol)
    assert pt.GetNumAtoms() == ct.GetNumAtoms(), \
        f"重原子数不一致: pose {pt.GetNumAtoms()} vs crystal {ct.GetNumAtoms()}"
    matches = ct.GetSubstructMatches(pt, uniquify=False, useChirality=False,
                                     useQueryQueryMatches=False, maxMatches=100000)
    assert matches, "拓扑图不匹配, 无法建立原子对应"
    # 坐标只取重原子（拓扑副本保留原子顺序, 但 idx 与含氢 mol 不同）
    keep_pose = [a.GetIdx() for a in pose_mol.GetAtoms() if a.GetSymbol() != "H"]
    keep_crys = [a.GetIdx() for a in crystal_mol.GetAtoms() if a.GetSymbol() != "H"]
    assert len(keep_pose) == pt.GetNumAtoms() and len(keep_crys) == ct.GetNumAtoms()
    P = np.array(pose_mol.GetConformer().GetPositions())[keep_pose]
    C = np.array(crystal_mol.GetConformer().GetPositions())[keep_crys]
    rmsds = []
    for match in matches:
        # match[i] = 晶体拓扑原子索引, 对应 pose 拓扑原子 i
        Cr = C[list(match)]
        rmsds.append(float(np.sqrt(np.mean(np.sum((P - Cr) ** 2, axis=1)))))
    return min(rmsds)

# ---------- 主流程 ----------
results = []
for case, cfg in CASES.items():
    tag = case
    print(f"\n===== {tag}: {cfg['resname']} redocking =====")
    # 1. 提取共晶配体（链 A, 重原子）
    lig_pdb = f"{OUTDIR}/{tag}_crystal_lig.pdb"
    n_atoms = extract_cocrystal_ligand(cfg["pdb"], cfg["chain"], cfg["resname"], lig_pdb)
    assert n_atoms > 0, f"{tag} 未提取到 {cfg['resname']}"
    print(f"  共晶配体重原子: {n_atoms}")
    # 2. 共晶配体质心核对（应与 step4_docking_site 一致）
    _, _, _ = 0, 0, 0
    xs, ys, zs = [], [], []
    for line in open(lig_pdb):
        if line.startswith("HETATM"):
            xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
    cx, cy, cz = sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)
    print(f"  配体质心: ({cx:.3f}, {cy:.3f}, {cz:.3f})  存档值: {cfg['center']}")
    assert abs(cx-cfg['center'][0]) < 1.0 and abs(cy-cfg['center'][1]) < 1.0 \
        and abs(cz-cfg['center'][2]) < 1.0, "质心与存档 box 中心不一致!"
    # 3. obabel PDB -> SDF（键级感知, 保持原子顺序, 不加氢）
    lig_sdf = f"{OUTDIR}/{tag}_crystal_lig.sdf"
    r = subprocess.run([OBABEL, lig_pdb, "-O", lig_sdf], capture_output=True, text=True)
    assert os.path.exists(lig_sdf) and os.path.getsize(lig_sdf) > 0, \
        f"obabel 失败: {r.stderr[:300]}"
    # 4. 晶体参照（obabel 键级感知的 SDF, 坐标即共晶构象）
    crystal_mol = Chem.MolFromMolFile(lig_sdf, removeHs=True, sanitize=True)
    assert crystal_mol is not None and crystal_mol.GetNumConformers() == 1
    print(f"  SDF 重原子: {crystal_mol.GetNumAtoms()}  "
          f"SMILES: {Chem.MolToSmiles(crystal_mol)}")
    # 5. 受体准备（镜像 step4_docking_site）
    rec_pdb = f"{OUTDIR}/{tag}_receptor_chainA.pdb"
    assert extract_protein_chain(cfg["pdb"], cfg["chain"], rec_pdb)
    rec_pdbqt = f"{OUTDIR}/{tag}_receptor.pdbqt"
    assert prep_receptor_pdbqt(rec_pdb, rec_pdbqt)
    # 6. 配体 pdbqt（镜像: AddHs + meeko）
    mol = Chem.MolFromMolFile(lig_sdf, removeHs=False)
    assert mol is not None
    mol = AllChem.AddHs(mol, addCoords=True)
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    assert setups, "meeko 配体准备失败"
    lig_pdbqt = f"{OUTDIR}/{tag}_lig.pdbqt"
    pdbqt_string, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    assert ok, f"meeko PDBQT 写出失败: {err}"
    with open(lig_pdbqt, "w") as f:
        f.write(pdbqt_string)
    # 7. vina redocking（镜像: box 22 A, exhaustiveness 16, 5 modes）
    out_pdbqt = f"{OUTDIR}/{tag}_out.pdbqt"
    cmd = [VINA, "--receptor", rec_pdbqt, "--ligand", lig_pdbqt,
           "--center_x", f"{cx:.2f}", "--center_y", f"{cy:.2f}", "--center_z", f"{cz:.2f}",
           "--size_x", f"{cfg['box']:.2f}", "--size_y", f"{cfg['box']:.2f}",
           "--size_z", f"{cfg['box']:.2f}",
           "--exhaustiveness", "16", "--num_modes", "5", "--out", out_pdbqt,
           "--seed", "2026"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert os.path.exists(out_pdbqt), f"vina 失败: {r.stderr[:400]}"
    affs = [float(l.split()[1]) for l in r.stdout.splitlines()
            if len(l.split()) >= 2 and l.split()[0].isdigit()]
    # 8. 各姿态 RMSD（meeko 姿态重建 + 图匹配原位对称校正 RMSD）
    pose_mol = pose_mols_from_pdbqt(out_pdbqt)
    n_pose = pose_mol.GetNumConformers()
    assert n_pose == len(affs), f"姿态数 {n_pose} != 亲和力数 {len(affs)}"
    pose_rmsds = []
    for cid in range(n_pose):
        pm = Chem.Mol(pose_mol)
        pm.RemoveAllConformers()
        conf = Chem.Conformer(pose_mol.GetConformer(cid))
        pm.AddConformer(conf, assignId=True)
        pose_rmsds.append((cid + 1, affs[cid], rmsd_inplace(pm, crystal_mol)))
    best = min(pose_rmsds, key=lambda t: t[2])
    top1 = pose_rmsds[0]
    print(f"  姿态 RMSD: " + ", ".join(
        f"#{m}({a}):{d:.3f}" for m, a, d in pose_rmsds))
    print(f"  Top-1 RMSD = {top1[2]:.3f} A | 最佳 RMSD = {best[2]:.3f} A"
          f" | 判定 {'PASS' if best[2] <= 2.0 else 'FAIL'} (<= 2.0 A)")
    results.append({"case": tag, "pdb": os.path.basename(cfg["pdb"]),
                    "cocrystal_ligand": cfg["resname"], "n_heavy_atoms": crystal_mol.GetNumAtoms(),
                    "top1_affinity_kcal_mol": top1[1], "top1_rmsd_A": round(top1[2], 3),
                    "best_affinity_kcal_mol": best[1], "best_rmsd_A": round(best[2], 3),
                    "protocol": "chain A; ligand-centroid box 22 A; exhaustiveness 16; seed 2026",
                    "pass_le_2A": "yes" if best[2] <= 2.0 else "no"})

out_csv = f"{OUTDIR}/redocking_results.csv"
with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
    w.writeheader(); w.writerows(results)
print("\n===== redocking 汇总 =====")
for r_ in results:
    print(f"  {r_['case']}: top1 {r_['top1_rmsd_A']} A, best {r_['best_rmsd_A']} A, "
          f"pass={r_['pass_le_2A']}")
print(f"已写出: {out_csv}")
