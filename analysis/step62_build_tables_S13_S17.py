#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 62 —— 生成补充表 TableS13-S17 xlsx（审稿回应新增）"""
import pandas as pd, os

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT = f"{BASE}/AJRCCM_submission/05_supplementary/Supplementary_Figures_Tables/Tables"

def build(path_sheets, out_name, notes=None):
    path = os.path.join(OUT, out_name)
    with pd.ExcelWriter(path, engine="openpyxl") as w:
        for sheet, src in path_sheets.items():
            df = pd.read_csv(src)
            df.to_excel(w, sheet_name=sheet, index=False)
        if notes:
            pd.DataFrame({"Notes": notes}).to_excel(w, sheet_name="README", index=False)
    print("saved:", path)

# ---- Table S13: permutation / apparent vs CV / exact LB ----
build({
    "permutation_500x": f"{BASE}/05_ML_SHAP/stability/permutation_500x_results.csv",
    "apparent_vs_cv_ARDS": f"{BASE}/05_ML_SHAP/stability/apparent_vs_cv_AUC_ARDS.csv",
    "apparent_vs_cv_Sepsis": f"{BASE}/05_ML_SHAP/stability/apparent_vs_cv_AUC_Sepsis.csv",
    "auc_exact_lower_bounds": f"{BASE}/05_ML_SHAP/stability/auc_exact_lower_bounds.csv",
}, "TableS13_ML_overfitting_diagnostics.xlsx", notes=[
    "Permutation test: 500 label shuffles per classifier (ARDS training cohort GSE76293); "
    "perm_P = (1 + #{null AUC >= real}) / (1 + 500). GBM could not be refit under permutation (NA).",
    "Exact one-sided lower confidence bound for observed AUC = 1.000: LB = alpha^(1/(n1*n0)).",
])

# ---- Table S14: feature-selection stability ----
build({
    "lasso_bootstrap_ARDS": f"{BASE}/05_ML_SHAP/stability/lasso_bootstrap_stability_ARDS.csv",
    "lasso_bootstrap_Sepsis": f"{BASE}/05_ML_SHAP/stability/lasso_bootstrap_stability_Sepsis.csv",
    "rf_seed_ARDS": f"{BASE}/05_ML_SHAP/stability/rf_seed_stability_ARDS.csv",
    "rf_seed_Sepsis": f"{BASE}/05_ML_SHAP/stability/rf_seed_stability_Sepsis.csv",
    "threshold_jaccard_ARDS": f"{BASE}/05_ML_SHAP/stability/threshold_perturbation_jaccard_ARDS.csv",
    "threshold_jaccard_Sepsis": f"{BASE}/05_ML_SHAP/stability/threshold_perturbation_jaccard_Sepsis.csv",
}, "TableS14_ML_feature_selection_stability.xlsx", notes=[
    "LASSO bootstrap: 200 resamples; selection frequency per gene.",
    "RF seed stability: 50 random seeds; frequency of top-30 and final-candidate inclusion.",
    "Threshold perturbation: Jaccard similarity of candidate sets under +/-10% ranking-threshold shifts.",
])

# ---- Table S15: MR power analysis ----
build({
    "power_summary": f"{BASE}/04_MR_SMR/MR_power_summary.csv",
    "power_by_OR": f"{BASE}/04_MR_SMR/MR_power_by_OR.csv",
    "IVs_forward": (f"{BASE}/04_MR_SMR/power_forward_sepsis_to_ARDS_IVs.csv"),
    "IVs_reverse": (f"{BASE}/04_MR_SMR/power_reverse_ARDS_to_sepsis_IVs.csv"),
}, "TableS15_MR_power_analysis.xlsx", notes=[
    "Analytic power (Brion 2013 / mRnd equivalent) with IVW SE under the multiplicative "
    "random-effects model as reported by TwoSampleMR; simulated power uses 5,000 draws.",
])

# ---- Table S16: validation diagnostics (W5) ----
build({
    "single_gene_AUC_CI": f"{BASE}/06_external_validation/reviewer_w5/single_gene_AUC_CI_validation.csv",
    "distribution_overlap": f"{BASE}/06_external_validation/reviewer_w5/distribution_overlap_train_vs_validation.csv",
    "metaDE_single_gene": f"{BASE}/06_external_validation/reviewer_w5/metaDE_single_gene_validation.csv",
    "combat_sensitivity": f"{BASE}/06_external_validation/reviewer_w5/combat_sensitivity_model_transfer.csv",
}, "TableS16_validation_diagnostics.xlsx", notes=[
    "OVL: within-cohort z-standardized overlap coefficient; Cohen d: within-cohort "
    "disease-vs-control effect; Hedges g: fixed- and random-effects meta-analysis across "
    "the two validation cohorts.",
    "ComBat sensitivity pools discovery and validation data and therefore breaks the "
    "validation firewall; exploratory only.",
])

# ---- Table S17: donor-level single-cell ----
build({
    "celloracle_per_donor": f"{BASE}/08_singlecell/celloracle_per_donor/CellOracle_CEBPB_KO_per_donor.csv",
    "UPP1_pseudotime_CEBPB_pb": f"{BASE}/08_singlecell/donor_level/GSE167363_donor_level_UPP1_CEBPB.csv",
    "pseudobulk_DESeq2_GSE167363": f"{BASE}/08_singlecell/donor_level/pseudobulk_DESeq2_GSE167363_Sepsis_vs_HC.csv",
    "pseudobulk_DESeq2_GSE151263": f"{BASE}/08_singlecell/donor_level/pseudobulk_DESeq2_GSE151263_ARDS_vs_Sepsis.csv",
}, "TableS17_donor_level_singlecell.xlsx", notes=[
    "CellOracle per donor: independent gene regulatory network built for each of the five "
    "GSE167363 donors (same parameters as the pooled analysis); delta = simulated minus "
    "observed mean UPP1 expression upon CEBPB in silico knockout.",
    "Pseudobulk DESeq2: raw counts aggregated per donor x cell type; GSE167363 contrast "
    "sepsis vs control (2 vs 2 donors); GSE151263 contrast ARDS vs sepsis (3 vs 4 donors).",
])
print("All done.")
