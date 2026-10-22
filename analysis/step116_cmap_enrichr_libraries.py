#!/usr/bin/env python3
"""
Step 116 -- CMap / DSigDB / Drug-Perturbations primary drug-repositioning query
(provenance closure, mirrors step115 for ChEMBL).

Background: the manuscript's drug-repositioning results (Supplementary Methods S12;
Supplementary Table S18; RESULTS, final section; Figure 10B) query the six hub genes
against the Enrichr API. The original query was run ad hoc; the per-library outputs
were archived as

    09_docking_drug/Old_CMAP_up_hub_prediction.csv
    09_docking_drug/Old_CMAP_down_hub_prediction.csv
    09_docking_drug/DSigDB_hub_prediction.csv
    09_docking_drug/Drug_Perturbations_from_GEO_2014_hub_prediction.csv
    (exploratory only, not reported in the manuscript:
     DisGeNET_hub_prediction.csv, Jensen_DISEASES_hub_prediction.csv)

This script re-issues the identical query protocol (same gene list, same libraries,
same API) so the archived CSVs are reproducible. Enrichr libraries are versioned
snapshots, so a re-run should return the same terms; if the live library has been
updated since 2026-09, compare against the archived files rather than overwriting.

Protocol (identical to step67_cmap_negative_control.py):
  1. POST the six gene symbols to https://maayanlab.cloud/Enrichr/addList
  2. GET  /enrich?userListId=<id>&backgroundType=<LIB>
  3. keep the top 30 terms per library, columns:
     rank, term, pvalue, zscore, combined_score, overlapping_genes (comma-joined)

Usage:
  python3 step116_cmap_enrichr_libraries.py [--out-dir 09_docking_drug]

Outputs (one CSV per library; suffix _rerun to avoid clobbering archives):
  <out>/Old_CMAP_up_hub_prediction_rerun.csv
  <out>/Old_CMAP_down_hub_prediction_rerun.csv
  <out>/DSigDB_hub_prediction_rerun.csv
  <out>/Drug_Perturbations_from_GEO_2014_hub_prediction_rerun.csv
Plus a comparison report against the archived files, printed to stdout and written to
  <out>/step116_cmap_reproduction_report.csv
"""
import argparse
import csv
import json
import os
import sys
import urllib.parse
import urllib.request

ENRICHR = "https://maayanlab.cloud/Enrichr"
HUB6 = ["CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1"]
# Libraries reported in Supplementary Table S18 (S12 text):
#   Old CMap up/down, DSigDB, Drug Perturbations from GEO.
LIBS = ["Old_CMAP_up", "Old_CMAP_down", "DSigDB", "Drug_Perturbations_from_GEO_2014"]
TOP_N = 30  # archived files hold 30 data rows


def add_list(genes):
    boundary = "----wb116"
    body = ("\r\n".join([
        f"--{boundary}",
        'Content-Disposition: form-data; name="list"',
        "", "\n".join(genes),
        f"--{boundary}",
        'Content-Disposition: form-data; name="description"',
        "", "step116 primary repositioning query",
        f"--{boundary}--", ""])).encode()
    req = urllib.request.Request(
        f"{ENRICHR}/addList", data=body,
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)["userListId"]


def enrich(uid, lib, retries=3):
    url = f"{ENRICHR}/enrich?userListId={uid}&backgroundType={urllib.parse.quote(lib)}"
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r).get(lib, [])
        except Exception as e:  # transient API failures; step115 used the same pattern
            if attempt == retries:
                raise
            print(f"  [retry {attempt}] {lib}: {e}", file=sys.stderr)


def write_library(terms, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rank", "term", "pvalue", "zscore", "combined_score", "overlapping_genes"])
        for t in terms[:TOP_N]:
            # Enrichr row: [rank, term, pvalue, zscore, combined, overlap, adjp, old_p, old_adjp]
            w.writerow([t[0], t[1], t[2], t[3], t[4], ",".join(t[5].split(";"))])


def compare(archived, rerun):
    """Term-level comparison of archived vs re-run top-10 (what S18 reports)."""
    def top10(p):
        with open(p) as fh:
            rows = list(csv.DictReader(fh))
        return [(r["term"], float(r["pvalue"])) for r in rows[:10]]
    a, b = top10(archived), top10(rerun)
    a_terms = {t for t, _ in a}
    b_terms = {t for t, _ in b}
    rank1_a = a[0][0] if a else ""
    rank1_b = b[0][0] if b else ""
    return {
        "n_top10_archived": len(a_terms),
        "n_top10_rerun": len(b_terms),
        "n_top10_shared": len(a_terms & b_terms),
        "rank1_archived": rank1_a,
        "rank1_rerun": rank1_b,
        "rank1_match": rank1_a == rank1_b,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=os.path.join(
        "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis",
        "09_docking_drug"))
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("Uploading the six hub genes to Enrichr ...")
    uid = add_list(HUB6)
    print(f"  userListId = {uid}")

    report = []
    for lib in LIBS:
        print(f"Querying {lib} ...")
        terms = enrich(uid, lib)
        rerun_path = os.path.join(args.out_dir, f"{lib}_hub_prediction_rerun.csv")
        write_library(terms, rerun_path)
        print(f"  wrote {rerun_path} ({min(len(terms), TOP_N)} rows)")

        archived_path = os.path.join(args.out_dir, f"{lib}_hub_prediction.csv")
        row = {"library": lib, "rerun_file": os.path.basename(rerun_path)}
        if os.path.exists(archived_path):
            row.update(compare(archived_path, rerun_path))
        else:
            row["note"] = "archived file not found; nothing to compare"
        report.append(row)

    out_report = os.path.join(args.out_dir, "step116_cmap_reproduction_report.csv")
    keys = ["library", "n_top10_archived", "n_top10_rerun", "n_top10_shared",
            "rank1_archived", "rank1_rerun", "rank1_match", "rerun_file", "note"]
    with open(out_report, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for row in report:
            w.writerow({k: row.get(k, "") for k in keys})
    print(f"\nReproduction report: {out_report}")
    for row in report:
        print("  {library}: top-10 shared {n_top10_shared}/{n_top10_archived}, "
              "rank-1 {rank1_archived} -> {rank1_rerun} (match={rank1_match})".format(**row))


if __name__ == "__main__":
    main()
