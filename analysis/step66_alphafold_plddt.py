#!/usr/bin/env python3
"""
Step 66 -- AlphaFold pLDDT structural-quality assessment for the four AF-modeled
docking targets (CDK5RAP2, GRB10, CD177, EXOSC4).

pLDDT is stored in the B-factor column of AlphaFold PDB files.
For each target we report:
  - overall mean/median pLDDT, % residues >= 70 (confident), % >= 90 (very high)
  - the same statistics restricted to residues whose CA lies inside the
    30 A docking box centered on the whole-protein centroid (as used in
    step4_docking_site.py)
Output: 09_docking_drug/alphafold_plddt_assessment.csv
"""
import csv, os

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/09_docking_drug"
TARGETS = {
    "CDK5RAP2": "receptors/CDK5RAP2_Q96SN8.pdb",
    "CD177": "receptors/CD177_Q8N6Q3.pdb",
    "GRB10": "receptors/GRB10_Q13322.pdb",
    "EXOSC4": "receptors/EXOSC4_Q9NPD3.pdb",
}
BOX = 30.0  # Angstrom, cubic box used for AF centroid docking

def read_ca(path):
    """Return list of (resseq, resname, x, y, z, plddt) for CA atoms."""
    cas = []
    with open(path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                cas.append((
                    int(line[22:26]), line[17:20].strip(),
                    float(line[30:38]), float(line[38:46]), float(line[46:54]),
                    float(line[60:66]),
                ))
    return cas

def stats(plddts):
    n = len(plddts)
    if n == 0:
        return 0, 0.0, 0.0, 0.0, 0.0
    s = sorted(plddts)
    mean = sum(s) / n
    median = s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2
    ge70 = 100.0 * sum(1 for v in s if v >= 70) / n
    ge90 = 100.0 * sum(1 for v in s if v >= 90) / n
    return n, mean, median, ge70, ge90

rows = []
for gene, rel in TARGETS.items():
    cas = read_ca(os.path.join(BASE, rel))
    cx = sum(c[2] for c in cas) / len(cas)
    cy = sum(c[3] for c in cas) / len(cas)
    cz = sum(c[4] for c in cas) / len(cas)
    half = BOX / 2.0
    in_box = [c for c in cas
              if abs(c[2] - cx) <= half and abs(c[3] - cy) <= half and abs(c[4] - cz) <= half]
    n_all, m_all, md_all, g70_all, g90_all = stats([c[5] for c in cas])
    n_box, m_box, md_box, g70_box, g90_box = stats([c[5] for c in in_box])
    rows.append({
        "gene": gene, "structure": os.path.basename(rel),
        "n_residues": n_all,
        "mean_pLDDT": round(m_all, 1), "median_pLDDT": round(md_all, 1),
        "pct_pLDDT_ge70": round(g70_all, 1), "pct_pLDDT_ge90": round(g90_all, 1),
        "n_residues_in_box": n_box,
        "box_mean_pLDDT": round(m_box, 1), "box_median_pLDDT": round(md_box, 1),
        "box_pct_pLDDT_ge70": round(g70_box, 1), "box_pct_pLDDT_ge90": round(g90_box, 1),
    })
    print(f"{gene}: overall mean {m_all:.1f} (>=70: {g70_all:.1f}%), "
          f"box mean {m_box:.1f} (>=70: {g70_box:.1f}%, n={n_box})")

out = os.path.join(BASE, "alphafold_plddt_assessment.csv")
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print("Saved", out)
