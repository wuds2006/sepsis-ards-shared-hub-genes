#!/usr/bin/env bash
# =============================================================================
# run_all.sh -- master analysis pipeline for
#   "Shared Hub Genes and Immune Signatures of Sepsis and ARDS:
#    A Multi-Omics Integration" (AJRCCM submission).
#
# Reproduces every RESULTS number from raw public data. Stage order follows the
# manuscript's Supplementary Methods (S1-S14); within a stage, scripts run in
# numeric order. A full run takes several days on a single workstation
# (single-cell CellOracle and docking stages dominate).
#
# Prerequisites
#   1. bash analysis/download_data.sh && bash analysis/download_data2.sh
#      plus the manual downloads listed in section 0 below.
#   2. R 4.6.1 with the packages in environment/R_packages.csv
#      (Rscript analysis/install_R_packages.R installs missing CRAN packages).
#   3. Python 3.13 with: scanpy, harmonypy, anndata, pandas, numpy, scipy,
#      rdkit, meeko, pymupdf, PIL, openpyxl, python-docx.
#   4. External binaries NOT run here (see the noted scripts):
#      SMR v1.03, PLINK, LDSC -- see step117_run_smr_heidi.sh (a documented
#      runbook; SMR/HEIDI outputs are archived under results/04_MR_SMR/),
#      AutoDock Vina 1.2.5 (invoked by the docking scripts below).
#
# Usage
#   bash analysis/run_all.sh [--stage N[,N,...]]   run selected stages
#   bash analysis/run_all.sh --list                list stages and scripts
#   bash analysis/run_all.sh                       run everything
#
# Stage scripts run from the project root (the parent of analysis/), because
# the archived scripts resolve inputs/outputs relative to the project root or
# use absolute BASE paths baked in at their original run paths. BEFORE running,
# set the five path variables in section 0 to your local layout.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ------------------------------------------------------------- 0. local paths
# The archived scripts were written for the original workstation and contain
# absolute paths. Adjust these to match your machine (see PATHS.md section 2
# for the full inventory and the recommended sed commands).
# export BASE="/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
# export PROJ="$BASE"
# export SUBM="$BASE/AJRCCM_submission"
# export REVIEW="$SUBM/09_review_deAI_2026-09-25"
BASE="${BASE:-placeholder-edit-me}"

R() { echo "  [R]    $*"; Rscript "$@"; }
P() { echo "  [py]   $*"; python3  "$@"; }

# ------------------------------------------------------------- stage registry
# name|scripts (semicolon-separated, relative to $BASE/scripts/)
STAGES=(
"S1_data_acquisition|download_data.sh;download_data2.sh;parse_geo_matrix.py"
"S8_microarray_preprocessing|step1a_DEG.R"
"S9_deg_wgcna|step1b_WGCNA.R;step1c_intersection.R;step55_wgcna_preservation.R;step27_reviewer_fix_R.R"
"S2_ml_shap|step2a_ML.R;step2a2_ML_sepsis.R;step50_ml_stability_overfit.R;step27_reviewer_fix_py.py"
"S6_external_validation|step2b_GSE100159_validate.R;step2b_validation_ARDS.R;step2b2_GSE100159_candidates.R;step2b3_GSE243066_sepsis_candidates.R;step2b4_GSE100159_ROC.R;step2b_AUC_CI.R;step49_cd177_ablation.R;step75_muc1_sensitivity.R;step52_validation_diagnostics.R;step59_cohort_pca.R;step2d_expression_ROC.R;step2d2_boxplot_ROC.R"
"S2_prediction_models|step2e_nomogram.R"
"S3_deconvolution|step41_cibersort_ssgsea.R;step45_hub6_gsva.R;step2c_hub_immune.R;step47_ssgsea_immune_4cohorts.R"
"S10_gsva_immune|step2c2_GSVA_immune.R;step2c3_GSVA_immune_top20.R"
"S4_singlecell|step3a_scanpy.py;step3a_singlecell.R;step3a_annotate.R;step3b_prepare_cellchat.py;step3b_cellchat.R;step3b_cellchat_compare.R;step3c_prepare_slingshot.py;step3c_slingshot.R;step3e_GSE167363_annotate.R;step3e_fix_annotate.R;step3e_fix2_annotate.R;step3e_fix3_annotate.R;step3f_GSE167363_cellchat.R;step3f_fix_cellchat_viz.R;step3f_fix2_cellchat_viz.R;step3g_aucell.R;step3h_GSE167363_slingshot.R;step3h_finish_slingshot.R;step3o_monocle3.R;step3n_seurat_integration.R;step3n_metabolism.py;step3n_clonality.py;step54_donor_level.R;step53_celloracle_per_donor.py"
"S4_cellchat_quantitative|step110_cellchat_netP_recalc.R;step111_cellchat_GSE167363_netP_recalc.R;step3m_cellchat_supp.R"
"S5_virtual_knockout|step3d_scTenifoldKnk.R;step3d_EXOSC4_ko.R;step3i_GSE167363_scTenifoldKnk.R;step3i_extra_knockout.R;step3d_CellOracle_local.py;step3d_CellOracle.py;step3j_GSE167363_CellOracle.py;step3n_celloracle_stratified.py;step3n_ko_stratified.R;step3n_ko_stratified_patch.R;step3n_ko_stratified_compare.R;step22_KO_top10_enrichment.R;step3l_cross_dataset_comparison.py"
"S7_mendelian_randomization|extract_MR_4980.py;match_ards_4980.py;extract_MR_py.py;step1d_MR.R;step1d_MR2.R;step1d_MR_4980.R;step48a_extract_reverse_MR.py;step48b_reverse_MR.R;step51_mr_power.R;prepare_ldsc.py;convert_gwas_ma.py;align_coloc.py;step1d_coloc.R;step1d_heidi.R"
# S7 SMR/HEIDI proper: NOT automated -- external SMR binary; see
# step117_run_smr_heidi.sh (runbook) and results/04_MR_SMR/SMR.txt (archived).
"S11_protein_validation|step_hpa_validation.R;step_hpa_immune_cell.R"
"S12_repositioning_docking|step67_cmap_negative_control.py;step76_cmap_deg_null.py;step116_cmap_enrichr_libraries.py;step68_opentargets_druggability.py;step115_chembl_known_drugs.py;step4_docking_site.py;step4_dock_one.py;step57_upp1_redock_fixed_box.py;step56_redocking.py;step66_alphafold_plddt.py;step70_mmgb_sasa.py;step69_plip_vnn1_luteolin.py;step5_network_ppi.R"
"S13_single_gene_gsea|step4a_single_gene_gsea.R;step46_single_gene_gsea_4cohorts.R"
"S14_prognosis_GSE65682|step105_survival_prognosis_GSE65682.R;step105b_pergene_KM_GSE65682.R;step106_extended_modeling_GSE65682.R;step109_three_statistical_facts.R"
"X_verify_results|verify_manuscript.py"
)

# ------------------------------------------------------------- CLI
run_stage() {
  local spec="$1" name="${1%%|*}" scripts="${1#*|}"
  echo "=================================================================="
  echo "STAGE $name"
  echo "=================================================================="
  IFS=';' read -ra arr <<< "$scripts"
  for s in "${arr[@]}"; do
    local f="scripts/$s"
    [ -f "$f" ] || { echo "  [skip] $s (not present)"; continue; }
    echo "--- $s"
    case "$s" in
      *.R)  R "$f" ;;
      *.sh) bash "$f" ;;
      *)    P "$f" ;;
    esac || { echo "  [FAILED] $s -- fix and re-run this stage"; return 1; }
  done
}

if [ "${1:-}" = "--list" ]; then
  for spec in "${STAGES[@]}"; do
    name="${spec%%|*}"
    n=$(grep -o ';' <<< "${spec#*|}" | wc -l | awk '{print $1+1}')
    printf "%-28s %s script(s)\n" "$name" "$n"
  done
  exit 0
fi

# Everything below executes analysis scripts: the project root must exist.
cd "$BASE" 2>/dev/null || { echo "Set BASE in run_all.sh (see section 0)."; exit 1; }

if [ "${1:-}" = "--stage" ]; then
  IFS=',' read -ra WANT <<< "${2:-}"
  for spec in "${STAGES[@]}"; do
    name="${spec%%|*}"
    for w in "${WANT[@]}"; do
      [ "$w" = "$name" ] && run_stage "$spec"
    done
  done
  exit 0
fi

# default: run all, in registry order
FAILED=()
for spec in "${STAGES[@]}"; do
  run_stage "$spec" || FAILED+=("${spec%%|*}")
done

echo "=================================================================="
if [ ${#FAILED[@]} -eq 0 ]; then
  echo "ALL STAGES COMPLETED."
  echo "Numeric outputs land in $BASE/01_DEG .. 13_survival_prognosis;"
  echo "manuscript-grade copies are under $BASE/AJRCCM_submission/."
else
  echo "COMPLETED WITH FAILURES in: ${FAILED[*]}"
  exit 1
fi
