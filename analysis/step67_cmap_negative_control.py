#!/usr/bin/env python3
"""
Step 67 -- CMap (Enrichr Old CMAP) negative-control specificity analysis.

Question (reviewer W1): is the glucocorticoid enrichment of the six-hub-gene
signature specific, or would any gene set score glucocorticoids?

Design:
  Real query     : the six hub genes (all up-regulated in both diseases) against
                   the Enrichr libraries Old_CMAP_up and Old_CMAP_down.
  Negative set A : 100 random 6-gene sets drawn from the 173 non-hub shared
                   candidate genes (same biological background).
  Negative set B : 100 random 6-gene sets drawn from the measured gene universe
                   of the discovery cohorts (genome-wide background).
Metric per set/library: presence of a glucocorticoid term in the top 10,
best (minimum) glucocorticoid P value, best overall P value.
Glucocorticoid terms are matched against a documented name list.
Outputs:
  09_docking_drug/cmap_negative_control_results.csv
  09_docking_drug/cmap_hub6_realtime_query.csv
  09_docking_drug/cmap_negative_control_summary.csv
"""
import csv, json, random, time, urllib.parse, urllib.request

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
ENRICHR = "https://maayanlab.cloud/Enrichr"
LIBS = ["Old_CMAP_up", "Old_CMAP_down"]
HUB6 = ["CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1"]
GLUCO = ["dexamethasone", "prednisolone", "prednisone", "hydrocortisone",
         "cortisone", "betamethasone", "triamcinolone", "methylprednisolone",
         "budesonide", "fludrocortisone", "corticosterone"]
random.seed(20260916)

def add_list(genes):
    boundary = "----wb67"
    body = ("\r\n".join([
        f"--{boundary}",
        'Content-Disposition: form-data; name="list"',
        "", "\n".join(genes),
        f"--{boundary}",
        'Content-Disposition: form-data; name="description"',
        "", "cmap_neg_control",
        f"--{boundary}--", ""])).encode()
    req = urllib.request.Request(f"{ENRICHR}/addList", data=body,
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)["userListId"]

def enrich(uid, lib):
    url = f"{ENRICHR}/enrich?userListId={uid}&backgroundType={urllib.parse.quote(lib)}"
    with urllib.request.urlopen(url, timeout=30) as r:
        return json.load(r).get(lib, [])

def evaluate(terms):
    """terms: list of [rank, term, pvalue, zscore, combined, overlap, adjp, ...]"""
    top10 = terms[:10]
    gluco_rows = [t for t in terms if any(g in t[1].lower() for g in GLUCO)]
    g_top10 = any(any(g in t[1].lower() for g in GLUCO) for t in top10)
    best_gluco_p = min((t[2] for t in gluco_rows), default=None)
    best_gluco_term = min(gluco_rows, key=lambda t: t[2])[1] if gluco_rows else ""
    best_p = min((t[2] for t in terms), default=None)
    return g_top10, best_gluco_p, best_gluco_term, best_p

# gene universes
candidates = [l.strip() for l in open(f"{BASE}/03_intersection/Shared_CGs.txt") if l.strip()]
bg_a = [g for g in candidates if g not in HUB6]
universe = set()
for f in ["01_DEG/GSE76293_DEG_all.csv", "01_DEG/GSE69528_DEG_all.csv"]:
    with open(f"{BASE}/{f}") as fh:
        r = csv.DictReader(fh)
        col = r.fieldnames[0]
        for row in r:
            universe.add(row[col])
bg_b = sorted(universe - set(HUB6))
print(f"background A: {len(bg_a)} non-hub candidates; background B: {len(bg_b)} measured genes")

results = []
# real hub6 query (fresh)
uid = add_list(HUB6); time.sleep(0.3)
real_rows = []
for lib in LIBS:
    terms = enrich(uid, lib); time.sleep(0.3)
    for t in terms:
        real_rows.append({"library": lib, "rank": t[0], "term": t[1],
                          "pvalue": t[2], "zscore": t[3], "combined_score": t[4],
                          "overlapping_genes": ";".join(t[5]) if isinstance(t[5], list) else t[5]})
    g10, gp, gt, bp = evaluate(terms)
    results.append({"set_id": "HUB6_REAL", "background": "real", "library": lib,
                    "gluco_in_top10": g10, "best_gluco_p": gp,
                    "best_gluco_term": gt, "best_p": bp})
    print(f"REAL {lib}: gluco top10={g10}, best gluco P={gp}, best P={bp:.2e}")

with open(f"{BASE}/09_docking_drug/cmap_hub6_realtime_query.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(real_rows[0].keys())); w.writeheader(); w.writerows(real_rows)

# negative controls
for bg_name, pool, n_sets in [("A_nonhub_candidates", bg_a, 100),
                              ("B_measured_universe", bg_b, 100)]:
    for i in range(n_sets):
        genes = random.sample(pool, 6)
        try:
            uid = add_list(genes); time.sleep(0.15)
            for lib in LIBS:
                terms = enrich(uid, lib); time.sleep(0.15)
                g10, gp, gt, bp = evaluate(terms)
                results.append({"set_id": f"{bg_name}_{i+1:03d}", "background": bg_name,
                                "library": lib, "gluco_in_top10": g10,
                                "best_gluco_p": gp, "best_gluco_term": gt, "best_p": bp})
        except Exception as e:
            print(f"  {bg_name} set {i+1} failed: {e}")
        if (i + 1) % 20 == 0:
            print(f"  {bg_name}: {i+1}/{n_sets} done")

out = f"{BASE}/09_docking_drug/cmap_negative_control_results.csv"
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(results[0].keys())); w.writeheader(); w.writerows(results)
print("Saved", out)

# summary
summary = []
for bg in ["real", "A_nonhub_candidates", "B_measured_universe"]:
    for lib in LIBS:
        sub = [r for r in results if r["background"] == bg and r["library"] == lib]
        n = len(sub)
        freq = 100.0 * sum(1 for r in sub if r["gluco_in_top10"]) / n if n else 0
        ps = [r["best_gluco_p"] for r in sub if r["best_gluco_p"] is not None]
        summary.append({"background": bg, "library": lib, "n_sets": n,
                        "gluco_top10_freq_pct": round(freq, 1),
                        "n_with_any_gluco_hit": len(ps)})
        print(f"{bg:22s} {lib:14s} n={n:3d} gluco-in-top10={freq:.1f}%  any-hit={len(ps)}")
with open(f"{BASE}/09_docking_drug/cmap_negative_control_summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0].keys())); w.writeheader(); w.writerows(summary)
print("Saved summary")
