#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# step108b_supp_tables.py
# 生成补充表 TableS22（模型比较）与 TableS23（LASSO 系数）CSV，放入投稿包。

import csv, io, os

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
TAB = f"{BASE}/13_survival_prognosis/tables"
DST = f"{BASE}/AJRCCM_submission/09_review_deAI_2026-09-25/figures/supplementary_tables"

def rd(fn):
    with io.open(os.path.join(TAB, fn), encoding="utf-8") as f:
        return list(csv.DictReader(f))

perf = rd("ext_model_performance.csv")
boot = rd("ext_bootstrap_optimism.csv")
lrt = rd("ext_lrt.csv")
nri = rd("ext_nri_idi.csv")

boot_d = {r["model"]: r for r in boot}
lrt_d = {"Clinical + Hub6": lrt[0], "Clinical + LASSO-179": lrt[1]}
nri_d = {"Clinical + Hub6": nri[0], "Clinical + LASSO-179": nri[1]}

rows = []
for r in perf:
    m = r["model"]
    b = boot_d.get(m, {})
    l = lrt_d.get(m, {})
    n = nri_d.get(m, {})
    rows.append({
        "model": m,
        "n_parameters": r["n_par"],
        "AIC": r["AIC"],
        "C_apparent": r["C_apparent"],
        "bootstrap_optimism": b.get("optimism", ""),
        "C_corrected": b.get("C_corrected", ""),
        "AUC_7d": r["AUC_7d"],
        "AUC_14d": r["AUC_14d"],
        "AUC_28d": r["AUC_28d"],
        "AUC28_95CI": r["AUC28_CI"],
        "LRT_chisq_vs_clinical": l.get("chisq", ""),
        "LRT_df": l.get("df", ""),
        "LRT_P": l.get("P", ""),
        "NRI_28d_vs_clinical": n.get("NRI", ""),
        "NRI_P": n.get("NRI_P", ""),
        "IDI_28d_vs_clinical": n.get("IDI", ""),
    })
with io.open(os.path.join(DST, "TableS22_prognostic_model_comparison.csv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

# Table S23：LASSO 系数 + 惩罚对照
lasso = rd("ext_lasso_coefficients.csv")
with io.open(os.path.join(DST, "TableS23_lasso_prognostic_signature.csv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["gene", "beta_standardised", "HR", "selected_hub"])
    w.writeheader(); w.writerows(lasso)

# 惩罚对照作为附加表（同 S23 说明）
pen = rd("ext_penalty_comparison.csv")
with io.open(os.path.join(DST, "TableS23b_penalty_comparison.csv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(pen[0].keys()))
    w.writeheader(); w.writerows(pen)

print("TableS22/S23 写入", DST)
