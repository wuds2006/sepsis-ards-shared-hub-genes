#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 70 —— MM/GBSA 单点自由能重评分（交叉学科专家 W2(b) 建议 3）
协议: rigid single-point MM/GBSA (bonded 项在 complex-protein-ligand 三体系中严格抵消)
  蛋白: OpenMM ff14SB + GBn2 (OBC-II, NoCutoff)
  配体: Gasteiger 电荷 + UFF LJ (元素表) + GBSAOBC2Force 标准半径 (getStandardParameters)
  非极性项: Shrake-Rupley SASA × gamma=0.005 kcal/mol/A^2 (KDTree 加速)
  塌缩防护: 极性 H 按 Amber 惯例零 LJ (ff14SB HO 类型 eps=0), 但配体无 bonded 项,
  eps=1 库仑吸引会在最小化中把带相反电荷的原子对推至 r->0 (实测 emodin/curcumin
  出现 -1e6 kcal/mol 量级伪能) → 最小化期间对含配体的原子对施加纯 r^-12 排斥墙
  (15 kcal/mol × (1.5 A/r)^12, 单调无势阱, 仅挡 <1.5 A 非物理接近), 评估系统不含墙。
背景: openforcefield 套件 2026 年从 PyPI/conda-forge 下架 (供应链事件), ambertools
  (600MB+) 本环境不可下载 → GAFF/AM1-BCC 不可用, 采用 UFF/Gasteiger 替代,
  Methods 与回应信中披露。
用法: python step70_mmgb_sasa.py VNN1 [UPP1 ...]   (按目标逐个跑, 每目标一个 CSV)
      环境变量 MMGB_LIGANDS=emodin,curcumin 可只跑指定配体 (快速验证用)
"""
import os, csv, sys, time
import numpy as np
import openmm
from openmm import app, unit
from openmm.app.internal import customgbforces
from pdbfixer import PDBFixer
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import rdPartialCharges
from scipy.spatial import cKDTree

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/09_docking_drug"
SITE = f"{BASE}/docking_site"
TMP  = f"{BASE}/_mmgb_tmp"
os.makedirs(TMP, exist_ok=True)

KCAL = 4.184
UFF_LJ = {"C": (3.851, 0.105), "H": (2.886, 0.044),
          "N": (3.660, 0.069), "O": (3.500, 0.060), "S": (4.035, 0.274)}
BONDI = {"C": 1.70, "N": 1.55, "O": 1.52, "H": 1.20, "S": 1.80}
GAMMA = 0.005
PROBE = 1.4
NPTS  = 960

TARGETS = {
    "UPP1":     (f"{BASE}/cocrystal/3EUF.pdb", "A"),
    "VNN1":     (f"{BASE}/cocrystal/7SLV.pdb", "A"),
    "CDK5RAP2": (f"{BASE}/receptors/CDK5RAP2_Q96SN8.pdb", None),
    "CD177":    (f"{BASE}/receptors/CD177_Q8N6Q3.pdb", None),
    "GRB10":    (f"{BASE}/receptors/GRB10_Q13322.pdb", None),
    "EXOSC4":   (f"{BASE}/receptors/EXOSC4_Q9NPD3.pdb", None),
}
LIGANDS = ["EGCG", "baicalein", "berberine", "curcumin", "emodin",
           "kaempferol", "luteolin", "quercetin", "resveratrol", "tetrahydroalstonine"]
if os.environ.get("MMGB_LIGANDS"):
    _sel = [x.strip() for x in os.environ["MMGB_LIGANDS"].split(",") if x.strip()]
    LIGANDS = [x for x in LIGANDS if x in _sel]
K_WALL_KJ = 62.76      # 15 kcal/mol
WALL_R0_NM = 0.15      # 1.5 A
FF = app.ForceField("amber14-all.xml")

def prep_protein(pdb_path, chain, out_path):
    lines = [l for l in open(pdb_path)
             if l.startswith("ATOM") and (chain is None or l[21] == chain)]
    raw = f"{TMP}/_raw_{os.path.basename(pdb_path)}"
    with open(raw, "w") as f:
        f.writelines(lines)
    fixer = PDBFixer(filename=raw)
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    mod = app.Modeller(fixer.topology, fixer.positions)
    mod.addHydrogens(FF, pH=7.4)
    with open(out_path, "w") as f:
        app.PDBFile.writeFile(mod.topology, mod.positions, f, keepIds=True)

def pose_mol(pdbqt_path):
    pq = PDBQTMolecule.from_file(pdbqt_path, skip_typing=True)
    return RDKitMolCreate.from_pdbqt_mol(pq)[0]

def lig_params(mol):
    rdPartialCharges.ComputeGasteigerCharges(mol)
    qs, sigs, epss = [], [], []
    for a in mol.GetAtoms():
        q = float(a.GetProp("_GasteigerCharge"))
        if np.isnan(q) or np.isinf(q):
            q = 0.0
        el = a.GetSymbol()
        x, D = UFF_LJ[el]
        # OpenMM LJ 形式 E=eps[(sig/r)^12-2(sig/r)^6], sig 即最小能量距离 = UFF 的 x
        sig_nm = x / 10.0
        eps_kJ = D * KCAL
        # 极性 H (连在 N/O/S 上) LJ 置零, 与 Amber ho/hn 类型惯例一致,
        # 否则 UFF 极性 H 的大 sigma 与蛋白 LJ 混合会在氢键距离处产生虚假排斥
        if el == "H" and any(nb.GetSymbol() in ("N", "O", "S") for nb in a.GetNeighbors()):
            sig_nm, eps_kJ = 0.0, 0.0
        qs.append(q)
        sigs.append(sig_nm)
        epss.append(eps_kJ)
    return qs, sigs, epss

def lig_exclusion_pairs(mol):
    """配体 1-2 与 1-3 原子对 (键连近邻), 非键相互作用须排除, 与常规力场一致"""
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    pairs = set(tuple(sorted(p)) for p in bonds)
    adj = {}
    for i, j in bonds:
        adj.setdefault(i, set()).add(j)
        adj.setdefault(j, set()).add(i)
    for i, j in bonds:
        for k in adj.get(j, ()):
            if k != i:
                pairs.add(tuple(sorted((i, k))))
        for k in adj.get(i, ()):
            if k != j:
                pairs.add(tuple(sorted((j, k))))
    return sorted(pairs)

def lig_topology(mol):
    top = app.Topology()
    ch = top.addChain()
    res = top.addResidue("LIG", ch)
    atoms = [top.addAtom(f"{a.GetSymbol()}{a.GetIdx()+1}",
                         app.Element.getBySymbol(a.GetSymbol()), res)
             for a in mol.GetAtoms()]
    for b in mol.GetBonds():
        top.addBond(atoms[b.GetBeginAtomIdx()], atoms[b.GetEndAtomIdx()])
    return top

def min_pl_dist(pos_p, pos_l):
    """蛋白-配体最小原子距 (nm), 用于塌缩诊断"""
    if len(pos_p) == 0 or len(pos_l) == 0:
        return float("nan")
    d, _ = cKDTree(pos_p).query(pos_l, k=1)
    return float(d.min())

def min_ll_dist(pos_l, excl):
    """配体内部最小非键 (非 1-2/1-3) 原子距 (nm)"""
    n = len(pos_l)
    dm = np.linalg.norm(pos_l[:, None, :] - pos_l[None, :, :], axis=-1)
    np.fill_diagonal(dm, np.inf)
    for i, j in excl:
        dm[i, j] = dm[j, i] = np.inf
    return float(dm.min())

def gb_std_params(top):
    arr = customgbforces.GBSAOBC2Force.getStandardParameters(top)
    return [(float(r), float(s)) for r, s in arr]

def find_forces(system):
    nbf = cgb = None
    for f in system.getForces():
        if isinstance(f, openmm.NonbondedForce):
            nbf = f
        elif isinstance(f, openmm.CustomGBForce):
            cgb = f
    return nbf, cgb

def eval_energy(system, positions):
    integ = openmm.VerletIntegrator(0.001)
    ctx = openmm.Context(system, integ)
    ctx.setPositions(positions)
    e = ctx.getState(getEnergy=True).getPotentialEnergy()
    del ctx, integ
    return e.value_in_unit(unit.kilojoule_per_mole) / KCAL

SPHERE = None
def fibonacci_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta),
                     np.cos(phi)], axis=1)

def sasa(coords, elements):
    """Shrake-Rupley, coords in Angstrom"""
    global SPHERE
    if SPHERE is None:
        SPHERE = fibonacci_sphere(NPTS)
    radii = np.array([BONDI.get(e, 1.70) for e in elements]) + PROBE
    tree = cKDTree(coords)
    pairs = tree.query_pairs(radii.max() * 2, output_type="ndarray")
    neigh = [[] for _ in range(len(coords))]
    for i, j in pairs:
        neigh[i].append(j); neigh[j].append(i)
    counts = np.zeros(len(coords))
    for i in range(len(coords)):
        pts = coords[i] + radii[i] * SPHERE
        nb = np.array(neigh[i], dtype=int)
        if len(nb):
            d2 = ((coords[nb][None, :, :] - pts[:, None, :]) ** 2).sum(-1)
            free = ~(d2 < radii[nb][None, :] ** 2).any(1)
        else:
            free = np.ones(NPTS, bool)
        counts[i] = free.sum()
    return (counts / NPTS * 4 * np.pi * radii ** 2).sum()

def run_target(gene):
    pdb_path, chain = TARGETS[gene]
    out_csv = f"{BASE}/mmgb_{gene}.csv"
    if os.path.exists(out_csv):
        print(gene, "already done"); return
    t0 = time.time()
    prot_pdb = f"{TMP}/{gene}_prot_H.pdb"
    if not os.path.exists(prot_pdb):
        prep_protein(pdb_path, chain, prot_pdb)
    pdb = app.PDBFile(prot_pdb)
    # OpenMM 8.5 起 createSystem 不再接受 implicitSolvent → 真空系统 + 手动 GBn2
    prot_system = FF.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff,
                                  constraints=None)
    nbf0, _ = find_forces(prot_system)
    gbc0 = customgbforces.GBSAOBC2Force(solventDielectric=78.5, soluteDielectric=1)
    for i, (r, s) in enumerate(gb_std_params(pdb.topology)):
        q = nbf0.getParticleParameters(i)[0]
        gbc0.addParticle([q.value_in_unit(unit.elementary_charge), r, s])
    gbc0.finalize()
    prot_system.addForce(gbc0)
    prot_pos = np.array(pdb.positions.value_in_unit(unit.nanometer))
    prot_elements = [a.element.symbol for a in pdb.topology.atoms()]
    E_prot = eval_energy(prot_system, prot_pos)
    SASA_prot = sasa(prot_pos * 10.0, prot_elements)
    xml = openmm.XmlSerializer.serialize(prot_system)
    print(f"{gene}: {len(prot_pos)} atoms, E_prot={E_prot:.1f} kcal/mol, "
          f"SASA={SASA_prot:.0f} A^2 [{time.time()-t0:.0f}s]", flush=True)

    rows = []
    for lig_name in LIGANDS:
        pose_path = f"{SITE}/{gene}_{lig_name}_out.pdbqt"
        if not os.path.exists(pose_path):
            continue
        mol = pose_mol(pose_path)
        lig_pos = np.array(mol.GetConformer().GetPositions()) / 10.0
        qs, sigs, epss = lig_params(mol)
        gbp = gb_std_params(lig_topology(mol))
        lig_elements = [a.GetSymbol() for a in mol.GetAtoms()]

        GB_OFFSET = 0.009   # customgbforces 落盘变换: r' = r - offset; s' = scale * r'
        excl = lig_exclusion_pairs(mol)
        nprot0 = len(prot_pos)
        def build_complex():
            sc_ = openmm.XmlSerializer.deserialize(xml)
            nbf_, cgb_ = find_forces(sc_)
            for q, sg, ep, (r, scc) in zip(qs, sigs, epss, gbp):
                sc_.addParticle(1.0)
                nbf_.addParticle(q, sg, ep)
                rt = r - GB_OFFSET
                cgb_.addParticle([q, rt, scc * rt])
            for i, j in excl:   # 配体键连对排除 (q=0, eps=0)
                nbf_.addException(nprot0 + i, nprot0 + j, 0.0, 0.5, 0.0)
            return sc_
        pos0 = np.vstack([prot_pos, lig_pos])

        # --- 复合物局部最小化 (配体位置约束 k=5000 kJ/mol/nm^2; 无配体 bonded FF,
        #     约束防止配体畸变; 蛋白侧链自由弛豫以解除对接姿态的 steric clash) ---
        # 塌缩防护: 纯 r^-12 排斥墙, 仅最小化期间生效。
        # 实现用 CustomBondForce 显式原子对 (配体 x 10A 内蛋白原子 + 配体内非 1-2/1-3 对):
        # 完全绕开 CustomNonbondedForce/GB 力的 exclusion 一致性限制,
        # 且 (1.5A/r)^12 在远处趋零, 有限对列举即完备。
        # 背景: 极性 H 按 Amber 惯例零 LJ (ff14SB HO eps=0), 无 bonded 防护的配体
        # 原子对会在 eps=1 库仑吸引下被最小化器推至 r->0 (库仑塌缩, 实测 -1e6 kcal/mol)。
        sys_c = build_complex()
        nlig = len(lig_pos)
        nprot = len(prot_pos)
        wall = openmm.CustomBondForce(f"Kwall*({WALL_R0_NM}/r)^12")
        wall.addGlobalParameter("Kwall", K_WALL_KJ)
        prot_tree = cKDTree(prot_pos)
        near_prot = prot_tree.query_ball_point(lig_pos, r=1.0)   # nm
        n_bonds = 0
        for li, plist in enumerate(near_prot):
            for pj in plist:
                wall.addBond(nprot + li, pj, [])
                n_bonds += 1
        excl_set = set(excl)
        for i in range(nlig):
            for j in range(i + 1, nlig):
                if (i, j) not in excl_set:
                    wall.addBond(nprot + i, nprot + j, [])
                    n_bonds += 1
        sys_c.addForce(wall)
        restraint = openmm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        restraint.addGlobalParameter("k", 5000.0)
        restraint.addPerParticleParameter("x0")
        restraint.addPerParticleParameter("y0")
        restraint.addPerParticleParameter("z0")
        nprot = len(prot_pos)
        for i in range(len(lig_pos)):
            restraint.addParticle(nprot + i, lig_pos[i])
        sys_c.addForce(restraint)
        integ = openmm.VerletIntegrator(0.001 * unit.picosecond)
        ctx = openmm.Context(sys_c, integ)
        ctx.setPositions(pos0)
        d0_pl = min_pl_dist(prot_pos, lig_pos)
        d0_ll = min_ll_dist(lig_pos, excl)
        openmm.LocalEnergyMinimizer.minimize(ctx, 10.0, 500)
        state = ctx.getState(getPositions=True)
        pos_c = np.array(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
        del ctx, integ, sys_c
        d1_pl = min_pl_dist(pos_c[:nprot], pos_c[nprot:])
        d1_ll = min_ll_dist(pos_c[nprot:], excl)
        print(f"    dmin(prot-lig) {d0_pl*10:.2f}->{d1_pl*10:.2f} A, "
              f"dmin(lig-lig) {d0_ll*10:.2f}->{d1_ll*10:.2f} A", flush=True)

        # 无 restraint 的新系统评估复合物能量 (单轨迹约定: 三体系同坐标)
        sys_e = build_complex()
        E_cplx = eval_energy(sys_e, pos_c)
        del sys_e

        sys_l = openmm.System()
        for _ in range(len(lig_pos)):
            sys_l.addParticle(1.0)
        nbl = openmm.NonbondedForce()
        nbl.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
        for q, sg, ep in zip(qs, sigs, epss):
            nbl.addParticle(q, sg, ep)
        for i, j in excl:
            nbl.addException(i, j, 0.0, 0.5, 0.0)
        sys_l.addForce(nbl)
        gbl = customgbforces.GBSAOBC2Force(solventDielectric=78.5, soluteDielectric=1)
        for q, (r, sc) in zip(qs, gbp):
            gbl.addParticle([q, r, sc])
        gbl.finalize()
        sys_l.addForce(gbl)
        E_lig = eval_energy(sys_l, pos_c[nprot:])
        del sys_l

        # 单轨迹约定: 蛋白/配体/SASA 均取最小化后复合物坐标
        E_prot_min = eval_energy(prot_system, pos_c[:nprot])
        SASA_cplx = sasa(pos_c * 10.0, prot_elements + lig_elements)
        SASA_prot_min = sasa(pos_c[:nprot] * 10.0, prot_elements)
        SASA_lig = sasa(pos_c[nprot:] * 10.0, lig_elements)
        dSASA = SASA_cplx - SASA_prot_min - SASA_lig
        dG = E_cplx - E_prot_min - E_lig + GAMMA * dSASA
        rows.append(dict(gene=gene, ligand=lig_name, E_complex=round(E_cplx, 2),
                         E_protein=round(E_prot_min, 2), E_ligand=round(E_lig, 2),
                         dSASA=round(dSASA, 1), deltaG=round(dG, 2)))
        print(f"  {gene}-{lig_name}: dG={dG:.2f} [{time.time()-t0:.0f}s]"
              + ("  [SUSPECT]" if dG < -150 else ""), flush=True)

    dock = {r["gene"] + "|" + r["ligand"]: float(r["affinity"])
            for r in csv.DictReader(open(f"{BASE}/docking_site_results.csv"))}
    for r in rows:
        r["vina"] = dock.get(r["gene"] + "|" + r["ligand"])
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print("saved", out_csv, flush=True)

if __name__ == "__main__":
    genes = sys.argv[1:] or list(TARGETS)
    for g in genes:
        run_target(g)
