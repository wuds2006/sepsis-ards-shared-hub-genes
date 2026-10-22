#!/usr/bin/env python3
"""对接姿态图：3D 结合口袋图 + 2D 相互作用分析（rdkit + matplotlib 替代 PyMOL）"""
import os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/09_docking_drug"
DOCK = f"{BASE}/docking_site"
OUT = f"{BASE}/poses"
os.makedirs(OUT, exist_ok=True)

# 最佳组合（结合能最低）
TOP = [
    ("VNN1", "luteolin"),
    ("VNN1", "quercetin"),
    ("UPP1", "curcumin"),
    ("GRB10", "EGCG"),
]

def parse_pdbqt(path, model=1):
    """解析 PDBQT，返回 (原子列表 [{name, x, y, z, elem, atype}])"""
    atoms = []
    cur_model = 0
    has_model = False
    with open(path) as f:
        for line in f:
            if line.startswith("MODEL"):
                cur_model = int(line.split()[1])
                has_model = True
                continue
            if line.startswith("ENDMDL"):
                if cur_model == model:
                    break
                continue
            if line.startswith(("ATOM", "HETATM")):
                if has_model and cur_model != model:
                    continue
                name = line[12:16].strip()
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                except: continue
                atype = line[77:79].strip() if len(line) >= 79 else ""
                elem = name[0] if name else line[76:78].strip()
                atoms.append({"name": name, "x": x, "y": y, "z": z, "elem": elem, "atype": atype})
    return atoms

def analyze_interactions(rec_atoms, lig_atoms, cutoff=4.0, hbond_cutoff=3.5):
    """分析配体-受体相互作用：接触残基 + 氢键"""
    contacts = []
    hbonds = []
    for la in lig_atoms:
        if la["elem"] == "H": continue
        for ra in rec_atoms:
            if ra["elem"] == "H": continue
            d = np.sqrt((la["x"]-ra["x"])**2 + (la["y"]-ra["y"])**2 + (la["z"]-ra["z"])**2)
            if d < cutoff:
                contacts.append((la, ra, d))
            # 氢键: 供体/受体为 O/N, 距离 < 3.5
            if d < hbond_cutoff and la["elem"] in ("O","N") and ra["elem"] in ("O","N"):
                hbonds.append((la, ra, d))
    return contacts, hbonds

def plot_pose(gene, lig_name, rec_atoms, lig_atoms, contacts, hbonds, aff):
    """画 3D 结合口袋图"""
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")
    # 受体（结合口袋，只画配体周围 8Å 内的）
    lig_center = np.mean([[a["x"],a["y"],a["z"]] for a in lig_atoms], axis=0)
    nearby = [a for a in rec_atoms if np.linalg.norm(np.array([a["x"],a["y"],a["z"]])-lig_center) < 10]
    if nearby:
        xs = [a["x"] for a in nearby]; ys = [a["y"] for a in nearby]; zs = [a["z"] for a in nearby]
        ax.scatter(xs, ys, zs, c="#8FA2D6", s=25, alpha=0.5, depthshade=True, label="Receptor pocket")
    # 配体（碳原子绿色高亮）
    for a in lig_atoms:
        if a["elem"] == "H": continue
        c = "#2E8B57" if a["elem"] == "C" else ("#D9534F" if a["elem"]=="O" else ("#5B7FA6" if a["elem"]=="N" else "grey"))
        ax.scatter(a["x"], a["y"], a["z"], c=c, s=60, depthshade=True)
    # 氢键（虚线）
    for la, ra, d in hbonds:
        ax.plot([la["x"],ra["x"]], [la["y"],ra["y"]], [la["z"],ra["z"]],
                color="#E7A889", linewidth=1.5, linestyle="--")
    ax.set_xlabel("X (Å)"); ax.set_ylabel("Y (Å)"); ax.set_zlabel("Z (Å)")
    ax.set_title(f"{gene} – {lig_name}\nBinding pose (affinity = {aff:.2f} kcal/mol, {len(hbonds)} H-bonds)")
    plt.tight_layout()
    plt.savefig(f"{OUT}/{gene}_{lig_name}_pose.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUT}/{gene}_{lig_name}_pose.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  3D 姿态图: {gene}_{lig_name}_pose.png")

def write_interactions(gene, lig_name, contacts, hbonds, aff):
    """写相互作用残基表"""
    # 接触残基（按距离排序，去重受体原子所属残基——简化：记录受体原子）
    contact_atoms = sorted(set((r["name"], r["elem"], round(d,2)) for l,r,d in contacts), key=lambda x: x[2])
    hb = [(l["name"], l["elem"], r["name"], r["elem"], round(d,2)) for l,r,d in hbonds]
    with open(f"{OUT}/{gene}_{lig_name}_interactions.txt", "w") as f:
        f.write(f"{gene} – {lig_name}  binding affinity = {aff:.2f} kcal/mol\n")
        f.write(f"氢键数: {len(hbonds)}  接触原子数: {len(contact_atoms)}\n\n")
        f.write("== 氢键 (配体原子 -> 受体原子, 距离Å) ==\n")
        for l_name, l_e, r_name, r_e, d in hb:
            f.write(f"  {l_name}({l_e}) -> {r_name}({r_e})  {d}Å\n")
        f.write("\n== 疏水/范德华接触 (受体原子, 距离Å) ==\n")
        for name, elem, d in contact_atoms[:30]:
            f.write(f"  {name}({elem})  {d}Å\n")

for gene, lig_name in TOP:
    rec_file = f"{DOCK}/{gene}_receptor.pdbqt"
    lig_file = f"{DOCK}/{gene}_{lig_name}_out.pdbqt"
    if not (os.path.exists(rec_file) and os.path.exists(lig_file)):
        print(f"跳过 {gene}-{lig_name}: 文件缺失"); continue
    rec_atoms = parse_pdbqt(rec_file)
    lig_atoms = parse_pdbqt(lig_file, model=1)
    # 从 out 文件读结合能
    aff = None
    with open(lig_file) as f:
        for line in f:
            if "VINA RESULT" in line:
                aff = float(line.split()[3]); break
    if aff is None:
        # 从 results csv 读
        import pandas as pd
        df = pd.read_csv(f"{DOCK}/{gene}_results.csv")
        row = df[(df.gene==gene) & (df.ligand==lig_name)]
        aff = row.affinity.values[0] if len(row) else None
    contacts, hbonds = analyze_interactions(rec_atoms, lig_atoms)
    print(f"\n{genes if False else gene} - {lig_name}: 结合能 {aff:.2f}, 接触 {len(contacts)}, 氢键 {len(hbonds)}")
    plot_pose(gene, lig_name, rec_atoms, lig_atoms, contacts, hbonds, aff)
    write_interactions(gene, lig_name, contacts, hbonds, aff)

print("\n===== 对接姿态图完成 =====")
