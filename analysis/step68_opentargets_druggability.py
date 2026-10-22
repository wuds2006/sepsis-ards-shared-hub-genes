#!/usr/bin/env python3
"""
Step 68 -- Open Targets druggability / tractability query for the six hub genes.
Uses curl for HTTP (sandbox-friendly). Outputs:
  09_docking_drug/opentargets_druggability.csv
  09_docking_drug/opentargets_known_drugs.csv
"""
import csv, json, subprocess, time

API = "https://api.platform.opentargets.org/api/v4/graphql"
GENES = ["CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1"]

Q_SEARCH = 'query($sym:String!){search(queryString:$sym,entityNames:["target"],page:{index:0,size:1}){hits{id name}}}'
Q_TARGET = '''query($ensg:String!){target(ensemblId:$ensg){
  approvedSymbol approvedName
  subcellularLocations{location source}
  tractability{modality label value}
  drugAndClinicalCandidates{count rows{maxClinicalStage drug{name drugType maximumClinicalStage mechanismsOfAction{rows{mechanismOfAction}}}}}
}}'''

def gql(query, variables):
    payload = json.dumps({"query": query, "variables": variables})
    r = subprocess.run(["curl", "-s", "-m", "60", "-X", "POST", API,
                        "-H", "Content-Type: application/json", "-d", payload],
                       capture_output=True, text=True)
    return json.loads(r.stdout)["data"]

PHASE = {0: "Preclinical", 1: "Phase I", 2: "Phase II", 3: "Phase III", 4: "Approved"}

tract_rows, drug_rows = [], []
for sym in GENES:
    hit = gql(Q_SEARCH, {"sym": sym})["search"]["hits"][0]
    ensg = hit["id"]; time.sleep(0.5)
    t = gql(Q_TARGET, {"ensg": ensg})["target"]; time.sleep(0.5)
    tract = {x["label"]: x["value"] for x in (t["tractability"] or [])}
    sm = [k for k, v in tract.items() if v]
    locs = sorted({l["location"] for l in (t["subcellularLocations"] or [])})
    kd = t["drugAndClinicalCandidates"]
    n_drugs = kd["count"] if kd else 0
    max_phase = -1; drug_names = []
    if kd:
        for row in kd["rows"]:
            drug = row["drug"]
            p = drug.get("maximumClinicalStage")
            p = p if p is not None else 0
            max_phase = max(max_phase, p)
            moa_rows = (drug.get("mechanismsOfAction") or {}).get("rows") or []
            moa = "; ".join(m.get("mechanismOfAction", "") for m in moa_rows[:2])
            drug_names.append((drug["name"], p, moa))
    top = sorted(drug_names, key=lambda x: -x[1])[:5]
    tract_rows.append({
        "gene": sym, "ensembl": ensg, "approved_name": t["approvedName"],
        "subcellular_locations": "; ".join(locs[:6]),
        "tractable_buckets": "; ".join(sm) if sm else "none",
        "n_tractable_buckets": len(sm),
        "n_known_drugs": n_drugs,
        "max_clinical_phase": PHASE.get(max_phase, "none") if n_drugs else "none",
        "top_known_drugs": "; ".join(f"{n} ({PHASE.get(p, p)})" for n, p, _ in top),
    })
    for n, p, moa in drug_names:
        drug_rows.append({"gene": sym, "drug": n, "phase": PHASE.get(p, p), "mechanism": moa})
    print(f"{sym}: {len(sm)} tractable buckets, {n_drugs} known drugs, "
          f"max phase {PHASE.get(max_phase, 'none') if n_drugs else 'none'}", flush=True)

with open("09_docking_drug/opentargets_druggability.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(tract_rows[0].keys()))
    w.writeheader(); w.writerows(tract_rows)
with open("09_docking_drug/opentargets_known_drugs.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["gene", "drug", "phase", "mechanism"])
    w.writeheader(); w.writerows(drug_rows)
print("Saved opentargets_druggability.csv and opentargets_known_drugs.csv")
