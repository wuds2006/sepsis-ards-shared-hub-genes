#!/usr/bin/env python3
"""生成 Figure S25B 用的 3Dmol 页面：VNN1(7SLV chain A) 卡通 + 口袋残基棒模型
+ luteolin 对接姿态 + 氢键虚线 + 接触残基标签。随后由 Chrome 无头渲染截图。
"""
import os, json
import numpy as np

BASE = os.environ.get("DOCK_ROOT", "09_docking_drug")
POSE = os.path.join(BASE, "docking_site", "VNN1_luteolin_out.pdbqt")
REC = os.path.join(BASE, "cocrystal", "7SLV.pdb")
WORK = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_figS25_work")

AD2EL = {"A": "C", "C": "C", "N": "N", "NA": "N", "OA": "O", "OS": "O", "O": "O",
         "SA": "S", "S": "S", "H": "H", "HD": "H", "HS": "H"}

# ---------- 受体 A 链（只要标准氨基酸，去掉共晶配体 9TC 与结晶水） ----------
AA3 = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
       "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
prot, rec_atoms = [], []
for line in open(REC):
    if not line.startswith("ATOM"):
        continue
    if line[21] != "A":
        continue
    if line[17:20].strip() not in AA3:
        continue
    el = line[76:78].strip() or line[12:16].strip()[0]
    if el == "H":
        continue
    prot.append(line.rstrip("\n"))
    rec_atoms.append(dict(resi=int(line[22:26]), resn=line[17:20].strip(),
                          aname=line[12:16].strip(), el=el,
                          xyz=np.array([float(line[30:38]), float(line[38:46]),
                                        float(line[46:54])])))

# ---------- 配体姿态 ----------
lig_lines, lig_atoms = [], []
cur, inm = None, False
for line in open(POSE):
    if line.startswith("MODEL"):
        cur = int(line.split()[1]); inm = (cur == 1); continue
    if line.startswith("ENDMDL"):
        if inm: break
        continue
    if inm and line.startswith(("ATOM", "HETATM")):
        ad = line[77:79].strip()
        el = AD2EL.get(ad, line[12:16].strip()[0])
        xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        lig_atoms.append(dict(el=el, xyz=np.array(xyz)))
        lig_lines.append((el, xyz))

# 配体坐标（3Dmol 用 HETATM；列宽必须严格对齐，否则解析出 0 个原子）
#   1-6 记录 7-11 序号 13-16 原子名 18-20 残基名 22 链 23-26 残基号 31-54 坐标 77-78 元素
lig_pdb = []
for i, (el, xyz) in enumerate(lig_lines, 1):
    lig_pdb.append(f"HETATM{i:5d} {el:>4s} LIG A 501    "
                   f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
                   f"  1.00  0.00          {el:>2s}")

# ---------- 口袋残基 ----------
lig_xyz = np.array([a["xyz"] for a in lig_atoms])
pocket = {}
for a in rec_atoms:
    d = np.min(np.linalg.norm(lig_xyz - a["xyz"], axis=1))
    if d <= 5.0:
        pocket.setdefault(a["resi"], a["resn"])
print("5 Å 口袋残基:", ", ".join(f"{v}{k}" for k, v in sorted(pocket.items())))

# ---------- 氢键端点 ----------
aa3 = {"ALA": "Ala", "ARG": "Arg", "ASN": "Asn", "ASP": "Asp", "CYS": "Cys",
       "GLN": "Gln", "GLU": "Glu", "GLY": "Gly", "HIS": "His", "ILE": "Ile",
       "LEU": "Leu", "LYS": "Lys", "MET": "Met", "PHE": "Phe", "PRO": "Pro",
       "SER": "Ser", "THR": "Thr", "TRP": "Trp", "TYR": "Tyr", "VAL": "Val"}

HB = [("GLU164", "OE2", 1), ("LYS157", "NZ", 2), ("MET218", "N", 4)]
# 配体原子序号（pdbqt 序列）→ 坐标；上面的 1/2/4 是 pdbqt 序号-1
hb_shapes = []
for resn, aname, lig_i in HB:
    r = next(a for a in rec_atoms if f"{a['resn']}{a['resi']}" == resn and a["aname"] == aname)
    lx = lig_atoms[lig_i]["xyz"]
    hb_shapes.append(dict(res=aa3[r["resn"]] + str(r["resi"]),
                          start=[round(v, 3) for v in lx],
                          end=[round(v, 3) for v in r["xyz"]]))

# ---------- 标签：仅标 7 个接触残基 ----------
CONTACT_RES = ["GLU164", "LYS157", "MET218", "TRP64", "PHE161", "PHE191", "TRP217"]
rows = []
for tag in CONTACT_RES:
    resi = int(tag[3:])
    resn = tag[:3]
    if resi not in pocket:
        print(f"[警告] {tag} 不在 5 Å 口袋集合内")
    # 默认锚定到离配体最近的侧链原子；个别残基改用 CA 以免标签贴边或叠压
    USE_CA = {"TRP217", "PHE191"}
    cand = [a for a in rec_atoms if a["resi"] == resi]
    if not cand:
        continue
    if tag in USE_CA:
        tgt = next((a for a in cand if a["aname"] == "CA"), cand[0])
    else:
        tgt = min(cand, key=lambda a: float(np.min(np.linalg.norm(lig_xyz - a["xyz"], axis=1))))
    rows.append(dict(label=aa3[resn] + str(resi),
                     pos=[round(v, 3) for v in tgt["xyz"]]))

payload = dict(protein="\n".join(prot) + "\nEND\n",
               ligcenter=[round(float(v), 3) for v in lig_xyz.mean(0)],
               ligand="\n".join(lig_pdb) + "\nEND\n",
               pocket=sorted(pocket),
               hbonds=hb_shapes,
               labels=rows)
with open(f"{WORK}/scene3d.json", "w") as f:
    json.dump(payload, f)

html = """<!DOCTYPE html><html><head><meta charset="utf-8">
<style>html,body{margin:0;padding:0;background:#fff;overflow:hidden}
#v{width:2400px;height:1800px;position:relative}</style>
<script src="3Dmol-min.js"></script></head><body>
<div id="v"></div>
<script>
const S = __PAYLOAD__;
const viewer = $3Dmol.createViewer(document.getElementById("v"),
             {backgroundColor: "#FFFFFF"});
viewer.addModel(S.protein, "pdb");
viewer.setStyle({model: 0}, {cartoon: {color: "#EAEFF4", opacity: 0.55, thickness: 0.35}});
viewer.addStyle({model: 0, resi: S.pocket},
                {stick: {radius: 0.10, colorscheme: "whiteCarbon", opacity: 1.0}});
viewer.addModel(S.ligand, "pdb");
viewer.setStyle({model: 1},
                {stick: {radius: 0.20, colorscheme: "greenCarbon"},
                 sphere: {scale: 0.20, colorscheme: "greenCarbon"}});
// 氢键：分段虚线
S.hbonds.forEach(function (h) {
  const a = h.start, b = h.end, n = 9;
  for (let k = 0; k < n; k++) {
    if (k % 2) continue;
    const t0 = k / n, t1 = (k + 1) / n;
    viewer.addCylinder({
      start: {x: a[0] + (b[0]-a[0])*t0, y: a[1] + (b[1]-a[1])*t0, z: a[2] + (b[2]-a[2])*t0},
      end:   {x: a[0] + (b[0]-a[0])*t1, y: a[1] + (b[1]-a[1])*t1, z: a[2] + (b[2]-a[2])*t1},
      radius: 0.045, color: "#B05B5B", fromCap: 0, toCap: 0
    });
  }
});
const ligC = S.ligcenter;
S.labels.forEach(function (l) {
  const dx = l.pos[0]-ligC[0], dy = l.pos[1]-ligC[1], dz = l.pos[2]-ligC[2];
  const n = Math.sqrt(dx*dx+dy*dy+dz*dz) || 1, k = 3.2/n;
  viewer.addLabel(l.label, {position: {x: l.pos[0]+dx*k, y: l.pos[1]+dy*k, z: l.pos[2]+dz*k},
    backgroundColor: "#FFFFFF", backgroundOpacity: 0.82,
    fontColor: "#2B2B2B", fontSize: 52, font: "Helvetica",
    showBackground: true, inFront: true,
    borderThickness: 0, padding: 4});
});
viewer.zoomTo({model: 1});
viewer.zoom(0.70);
viewer.render();
viewer.render();
setTimeout(function(){ document.title = "READY"; }, 300);
</script></body></html>"""
html = html.replace("__PAYLOAD__", json.dumps(payload))
with open(f"{WORK}/scene3d_template.html", "w") as f:
    f.write(html)
print("写出 scene3d.json 与 scene3d_template.html")
print("口袋残基:", len(pocket), "氢键:", len(hb_shapes), "标签:", len(rows))
