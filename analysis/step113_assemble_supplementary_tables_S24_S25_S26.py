#!/usr/bin/env python3
"""step113: assemble the three supplementary tables that the manuscript already cites
but that had no physical file (S24, S25, S26), from existing analysis outputs.

  S24 - prognostic statistical checks (three facts added in response to reviewer items):
        A. Cox proportional-hazards assumption (global and per gene), GSE65682
        B. Kaplan-Meier cut-point sensitivity
        C. ssGSEA immune-comparison correction coverage (BH)
  S25 - CellChat pathway-level communication strength, recomputed for BOTH single-cell
        cohorts (netP$strength was NULL in the stored objects; see step110 / step111)
  S26 - cohort metadata (ancestry and sampling time point flagged for author input)

Sources are read, never recomputed, except the ssGSEA summary counts.
Date: 2026-10-19
"""
import csv
import os

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "AJRCCM_submission/09_review_deAI_2026-09-25/figures/supplementary_tables")
os.makedirs(OUT, exist_ok=True)


def read_csv(path):
    with open(os.path.join(ROOT, path), newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


# ---------------- S24 ----------------
s24 = []
s24.append(["A", "Cox PH assumption, GLOBAL (six-gene joint model)",
            "chi-square; df; P", "13.87; 6; 0.031",
            "Global test rejects proportionality; the joint estimate is a time-averaged effect."])
for r in read_csv("13_survival_prognosis/tables/step109_coxzph_hub6.csv"):
    term = r["term"]
    if term == "GLOBAL":
        continue
    s24.append(["A", f"Cox PH assumption, {term.replace('z_', '')}",
                "chi-square; df; P",
                f"{float(r['chisq']):.2f}; {r['df']}; {float(r['p']):.3f}",
                "No single gene violates proportionality at P<0.05."])
for r in read_csv("13_survival_prognosis/tables/step109_cutpoint_sensitivity.csv"):
    hr = r["HR"]
    hr_s = "NA" if hr in ("NA", "") else f"{float(hr):.2f}"
    lo = "NA" if r["lo"] in ("NA", "") else f"{float(r['lo']):.2f}"
    hi = "NA" if r["hi"] in ("NA", "") else f"{float(r['hi']):.2f}"
    s24.append(["B", f"KM cut-point scheme: {r['scheme']}",
                "HR (95% CI); log-rank P",
                f"{hr_s} ({lo} to {hi}); {float(r['P']):.5f}",
                f"n groups = {r['n_groups']}"])
for path, lab in [("07_hub_immune/GSE76293_immune_diff_BH.csv", "GSE76293 (ARDS)"),
                  ("07_hub_immune/GSE69528_immune_diff_BH.csv", "GSE69528 (sepsis)")]:
    rows = read_csv(path)
    n = len(rows)
    nom = sum(1 for r in rows if float(r["p"]) < 0.05)
    bh = sum(1 for r in rows if float(r["padj"]) < 0.05)
    s24.append(["C", f"ssGSEA lineage-marker comparison, {lab}",
                "cell types tested; nominal P<0.05; BH-adjusted P<0.05",
                f"{n}; {nom}; {bh}",
                "Benjamini-Hochberg applied across the tested cell types within each cohort."])

with open(os.path.join(OUT, "TableS24_prognostic_statistical_checks.csv"), "w",
          newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["block", "item", "statistic", "value", "note"])
    w.writerows(s24)

# ---------------- S25 ----------------
s25 = []
for r in read_csv("08_singlecell/cellchat/pathway_diff_network_recalc_2026-10-19.csv"):
    s25.append(["GSE151263 (PBMC)", "sepsis versus ARDS", r["pathway"],
                r["strength_ARDS"], r["strength_Sepsis"], r["diff_Sepsis_minus_ARDS"],
                r["ratio_Sepsis_over_ARDS"], r["dominant_side"],
                r["weight_LR_pairs_ARDS"], r["weight_LR_pairs_Sepsis"], r["rank_by_abs_diff"]])
for r in read_csv("08_singlecell/GSE167363/cellchat/pathway_diff_network_GSE167363_recalc_2026-10-19.csv"):
    s25.append(["GSE167363 (CD14+ monocytes)", "sepsis versus healthy control", r["pathway"],
                r["strength_HC"], r["strength_Sepsis"], r["diff_Sepsis_minus_HC"],
                r["ratio_Sepsis_over_HC"], r["dominant_side"],
                r["weight_LR_pairs_HC"], r["weight_LR_pairs_Sepsis"], r["rank_by_abs_diff"]])

with open(os.path.join(OUT, "TableS25_CellChat_pathway_strength_recalc.csv"), "w",
          newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["cohort", "comparison", "pathway",
                "strength_group1", "strength_group2",
                "diff_group2_minus_group1", "ratio_group2_over_group1",
                "dominant_side", "weight_LR_pairs_group1", "weight_LR_pairs_group2",
                "rank_by_abs_diff"])
    w.writerows(s25)

# ---------------- S26 ----------------
NEED = "[AUTHOR_INPUT_NEEDED]"
s26 = [
    ["GSE76293", "Discovery (ARDS)", "GPL570, Affymetrix HG-U133 Plus 2.0",
     "Peripheral-blood polymorphonuclear neutrophils", "ARDS", "12", "12", "24",
     NEED + ": ancestry", NEED + ": sampling time point",
     "70 in-vitro cultured samples removed before analysis"],
    ["GSE69528", "Discovery (sepsis)", "GPL10558, Illumina HumanHT-12 v4",
     "Whole blood", "Sepsis", "83", "28", "111",
     NEED + ": ancestry", NEED + ": sampling time point",
     "27 non-infected type 2 diabetes samples removed before analysis"],
    ["GSE100159", "External validation", "GPL6884, Illumina HumanWG-6 v3.0",
     "Whole blood", "Sepsis", "33", "12", "45",
     NEED + ": ancestry", NEED + ": sampling time point", ""],
    ["GSE243066", "External validation", "GPL30209, Illumina NovaSeq 6000 (RNA-seq)",
     "Peripheral-blood mononuclear cells", "ARDS", "34", "15", "49",
     NEED + ": ancestry", NEED + ": sampling time point", ""],
    ["GSE65682", "Prognostic (septic shock)", "GPL570, Affymetrix HG-U133 Plus 2.0",
     "Whole blood", "Septic shock", "479", "0", "479",
     "European (as reported by the original study; " + NEED + ": confirm)",
     "At intensive care unit admission",
     "28-day mortality, 114 deaths; Mars1-Mars4 transcriptomic endotypes"],
    ["GSE151263", "Single cell", "10x Genomics", "Peripheral-blood mononuclear cells",
     "ARDS and sepsis", "3 ARDS + 4 sepsis donors", "0", "13,996 cells",
     NEED + ": ancestry", NEED + ": sampling time point",
     "Disease-versus-disease contrast at the patient level (n=7)"],
    ["GSE167363", "Single cell", "10x Genomics", "Peripheral-blood CD14+ monocytes",
     "Sepsis", "2 sepsis donors", "2 healthy donors", "12,359 cells",
     NEED + ": ancestry", NEED + ": sampling time point",
     "One donor sampled twice"],
    ["ieu-b-4980 (IEU OpenGWAS; GWAS Catalog GCST90281174)", "Mendelian randomization",
     "GWAS summary statistics", "Germline", "Sepsis", "NA", "NA", "NA",
     "European", "NA", "Exposure for the forward MR analysis"],
    ["FinnGen R11", "Mendelian randomization", "GWAS summary statistics", "Germline",
     "ARDS", "NA", "NA", "NA", "European", "NA", "Outcome for the forward MR analysis"],
    ["eQTLGen", "SMR-HEIDI reference", "Blood cis-eQTL summary statistics", "Whole blood",
     "NA", "NA", "NA", "NA", "European", "NA",
     "*EXOSC4* has no usable cis-eQTL in this reference"],
]

with open(os.path.join(OUT, "TableS26_cohort_metadata.csv"), "w",
          newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["dataset", "role", "platform", "tissue_or_cell_type", "disease",
                "cases", "controls", "total", "ancestry", "sampling_time_point", "notes"])
    w.writerows(s26)

print("S24 rows:", len(s24))
print("S25 rows:", len(s25))
print("S26 rows:", len(s26))
for f in ["TableS24_prognostic_statistical_checks.csv",
          "TableS25_CellChat_pathway_strength_recalc.csv",
          "TableS26_cohort_metadata.csv"]:
    print("wrote", os.path.join(OUT, f), os.path.getsize(os.path.join(OUT, f)), "bytes")
