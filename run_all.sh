#!/usr/bin/env bash
# =============================================================================
# run_all.sh -- 全流程驱动脚本
#   Shared hub genes and immune signatures of sepsis and ARDS (AJRCCM)
#
# 说明
#   脚本按阶段编号顺序执行，阶段与 Supplementary Methods S1-S14 的对应关系
#   见 docs/METHODS_TO_CODE.md。所有脚本的工作目录都是仓库根目录，
#   路径由 R/00_common.R 的 ROOT 解析，不需要改脚本。
#
# 用法
#   bash run_all.sh                 # 全流程
#   bash run_all.sh --stage 1,3     # 只跑指定阶段
#   bash run_all.sh --list          # 列出阶段与脚本
#
# 依赖
#   R 4.6.1 + environment/R_packages_used.txt 中的包
#   Python 3.13 + environment/python_packages_used.txt 中的包
#   SMR 1.03 / PLINK 1.9 / AutoDock Vina 1.2.5 / OpenMM（仅第 6、8 阶段用到）
# =============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE" || exit 1

R()  { echo "  [R]  $*"; Rscript "$@"; }
PY() { echo "  [py] $*"; python3 "$@"; }
SH() { echo "  [sh] $*"; bash "$@"; }

stage_0() {   # S1 数据获取与样本清洗
  SH shell/01_fetch_geo_series.sh
  SH shell/02_fetch_annotation_and_gwas.sh
  PY python/01_parse_geo_matrix.py
}

stage_1() {   # S8/S9 预处理、差异表达、WGCNA、交集与富集
  R R/10_deg.R
  R R/11_wgcna.R
  R R/12_wgcna_preservation.R
  R R/13_intersection.R
  R R/14_enrichment_go_kegg.R
  R R/15_ppi_network.R
}

stage_2() {   # S2 机器学习与 SHAP
  R R/20_ml_ards.R
  R R/21_ml_sepsis.R
  R R/22_ml_stability.R
}

stage_3() {   # S6 外部验证与 meta 分析
  R R/23_validation_ards.R
  R R/24_validation_sepsis_candidates.R
  R R/25_validation_ards_candidates.R
  R R/26_validation_sepsis_roc.R
  R R/27_validation_auc_ci.R
  R R/28_validation_diagnostics.R
  R R/29_muc1_sensitivity.R
  R R/30_cd177_ablation.R
}

stage_4() {   # S3/S10/S13 免疫浸润、GSVA、单基因 GSEA
  R R/35_cibersort_ssgsea.R
  R R/36_gsva_immune.R
  R R/37_hub6_gsva.R
  R R/38_single_gene_gsea.R
  R R/39_ssgsea_4cohorts.R
  R R/40_cohort_pca.R
}

stage_5() {   # S14 GSE65682 预后建模
  R R/45_prognosis_cox.R
  R R/46_prognosis_pergene_km.R
  R R/47_prognosis_extended_models.R
  R R/48_prognosis_figure.R
  R R/49_prognosis_statistics.R
}

stage_6() {   # S7 孟德尔随机化、coloc、SMR-HEIDI
  SH shell/03_prepare_mr_inputs.sh
  PY python/02_mr_extract_forward.py
  PY python/03_convert_gwas_sumstats.py
  PY python/04_match_ards_against_sepsis.py
  PY python/05_align_coloc_inputs.py
  PY python/06_prepare_ldsc_inputs.py
  R R/50_mr_forward.R
  PY python/07_mr_extract_reverse.py
  R R/51_mr_reverse.R
  R R/52_mr_power.R
  R R/53_colocalization.R
  R R/54_smr_heidi.R
  echo "  ! SMR/LDSC 需外部二进制，见 shell/04_run_smr_heidi.sh"
}

stage_7() {   # S4/S5 单细胞、通讯、拟时序、虚拟敲除
  PY python/10_scanpy_cluster.py
  R R/60_annotate_gse167363.R
  PY python/11_cellchat_prepare.py
  R R/61_cellchat_pbmc.R
  R R/62_cellchat_netp_pbmc.R
  R R/63_cellchat_monocytes.R
  R R/64_cellchat_netp_monocytes.R
  R R/65_monocle3.R
  R R/66_aucell.R
  R R/67_sctenifold_six_hub.R
  R R/68_sctenifold_monocytes.R
  R R/69_sctenifold_stratified.R
  R R/70_donor_level.R
  PY python/12_celloracle_pbmc.py
  PY python/13_celloracle_monocytes.py
  PY python/14_celloracle_stratified.py
  PY python/15_celloracle_per_donor.py
  PY python/16_metabolism.py
  PY python/17_cross_dataset.py
}

stage_8() {   # S12 分子对接与重打分
  PY python/21_docking_box_definition.py
  PY python/20_docking_screen.py
  PY python/22_redock_upp1.py
  PY python/23_redocking_validation.py
  PY python/24_mmgbsa_rescoring.py
  PY python/25_alphafold_plddt.py
  PY python/26_plip_contacts.py
}

stage_9() {   # S12 药物重定位与可成药性
  PY python/30_cmap_random_null.py
  PY python/31_cmap_deg_matched_null.py
  PY python/32_cmap_enrichr.py
  PY python/33_chembl_known_drugs.py
  PY python/34_opentargets_druggability.py
}

stage_10() {  # S11 蛋白水平验证
  R R/75_hpa_plasma_protein.R
  R R/76_hpa_immune_cells.R
}

stage_11() {  # 出图：主图 1-10 + 补充图 S1-S27
  R R/figures/01_prep_tables.R
  for f in R/figures/0*_fig*.R R/figures/10_fig10.R; do R "$f"; done
  for f in R/figures/*_supp_S*.R; do R "$f"; done
  R R/figures/99_check_figures.R
}

list_stages() {
  sed -n 's/^stage_\([0-9]*\)() *{ *# *\(.*\)/  \1  \2/p' "$0"
}

if [ "${1:-}" = "--list" ]; then
  echo "阶段 / stage:"; list_stages; exit 0
fi

WANT="${2:-all}"
if [ "${1:-}" = "--stage" ]; then WANT="$2"; else WANT="all"; fi

for s in 0 1 2 3 4 5 6 7 8 9 10 11; do
  if [ "$WANT" != "all" ]; then
    case ",$WANT," in *",$s,"*) ;; *) continue ;; esac
  fi
  echo "== stage $s =="
  "stage_$s" || echo "  ! stage $s 中断，请检查上面的报错"
done
echo "done"
