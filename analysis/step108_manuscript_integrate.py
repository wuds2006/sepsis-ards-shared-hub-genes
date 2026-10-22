#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# step108_manuscript_integrate.py
# 将 GSE65682 预后分析与扩展建模整合入手稿：
#   02_Manuscript_Revised.md（英文）、03_Manuscript_Translation_CN.md（中文）、
#   Supplementary_Materials.md、Supplementary_Methods.md、04_Modification_Log.md
# 数字来自 13_survival_prognosis/（step105/106/107）。

import io, sys

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
PKG = f"{BASE}/AJRCCM_submission/09_review_deAI_2026-09-25"

def sub_once(text, old, new, path):
    n = text.count(old)
    assert n == 1, f"锚点计数异常({n}) in {path}: {old[:60]!r}"
    return text.replace(old, new, 1)

# ============================================================
# 英文修订稿 02_Manuscript_Revised.md
# ============================================================
p2 = f"{PKG}/02_Manuscript_Revised.md"
t = io.open(p2, encoding="utf-8").read()

# 1. Abstract Methods 句
t = sub_once(t,
    "Findings were extended by immune deconvolution, single-cell transcriptomics with dual virtual knockout, drug repurposing prediction and molecular docking, and Mendelian randomization.",
    "Findings were extended by immune deconvolution, single-cell transcriptomics with dual virtual knockout, drug repurposing prediction and molecular docking, Mendelian randomization, and prognostic evaluation in an independent septic shock cohort.", p2)

# 2. Abstract M&M 段末
t = sub_once(t,
    "whereas the reverse analysis was informative (minimum detectable odds ratio 1.04).",
    "whereas the reverse analysis was informative (minimum detectable odds ratio 1.04). In the prognostic cohort (479 patients with septic shock, 114 deaths within 28 days), the six-gene score stratified mortality (log-rank P=0.002) but added nothing to age, sex, and endotype (likelihood-ratio P=0.16), whereas a penalized signature drawn from the 179 shared candidates retained prognostic information (optimism-corrected C index 0.622) and contained none of the six hub genes.", p2)

# 3. Abstract Conclusions 句
t = sub_once(t,
    "they are candidate disease-state markers requiring experimental validation.",
    "they are candidate disease-state markers rather than prognostic markers and require experimental validation.", p2)

# 4. Data Sources 段末
t = sub_once(t,
    "Single-cell analysis used GSE151263 (PBMCs; three ARDS and four sepsis donors; approximately 14,000 cells) and GSE167363 (CD14+ monocytes; two healthy controls and two sepsis, one sampled at two time points; 12,359 cells).",
    "Single-cell analysis used GSE151263 (PBMCs; three ARDS and four sepsis donors; approximately 14,000 cells) and GSE167363 (CD14+ monocytes; two healthy controls and two sepsis, one sampled at two time points; 12,359 cells). Prognostic evaluation used GSE65682 (GPL13667; whole blood from 479 patients with septic shock of the Molecular Diagnosis and Risk Stratification of Sepsis [MARS] study, sampled within 24 hours of intensive-care-unit admission), the only available cohort with a survival endpoint (28-day mortality; 114 deaths, 23.8%), with age, sex, and Mars1-to-Mars4 endotype assignments (39).", p2)

# 5. Methods 新小节（Mendelian Randomization 之前）
methods_new = """Prognostic Modeling in an Independent Septic Shock Cohort

Prognostic analysis was confined to GSE65682, the only cohort with a survival endpoint. Expression was taken from the processed series matrix; probes were collapsed to gene symbols by the probe with the highest mean expression, and 171 of the 179 shared candidate genes were represented on the platform. The endpoint was 28-day mortality; patients with a recorded follow-up time of zero were assigned 0.5 days. Six-hub-gene Cox models were fitted with expression standardized to z scores, separately for each gene and jointly for all six; the multivariable model's linear predictor was used as a risk score, and patients were grouped at its median for Kaplan-Meier curves (log-rank test). Time-dependent discrimination was assessed by inverse-probability-of-censoring-weighted receiver-operating-characteristic curves at 7, 14, and 28 days (40), with the 28-day evaluation placed at 27.99 days to avoid the administrative censoring boundary. To test whether prognostic information exists beyond the six hub genes, a clinical model (age, sex, and Mars endotype) was compared with a model adding the six-gene score and with a penalized Cox model over the 171 shared candidate genes fitted by LASSO (glmnet) with the penalty selected by 10-fold cross-validation (41); because the one-standard-error rule retained no gene, the penalty at minimum cross-validated deviance was used, which retained 15 genes, and elastic-net and ridge penalties were fitted for comparison. A random survival forest (ranger, 500 trees) over the same genes plus clinical variables served as a non-linear comparator, with individual-level risk taken from out-of-fold predictions of 10-fold cross-validation because in-bag predictions were strongly optimistic at this sample size (42). Discrimination was summarized by the Harrell C index with optimism correction from 300 bootstrap resamples; within each resample the LASSO penalty was held at the full-cohort value, so the optimism estimate for the penalized models is, if anything, conservative. Calibration was assessed in five risk groups as mean predicted versus Kaplan-Meier observed 28-day mortality. Incremental value over the clinical model was tested by likelihood-ratio tests and quantified by the net reclassification index and integrated discrimination improvement at 28 days (43), and clinical usefulness was assessed by decision-curve analysis over threshold probabilities of 0.05 to 0.60.

"""
t = sub_once(t, "Mendelian Randomization\n\nBidirectional two-sample MR",
             methods_new + "Mendelian Randomization\n\nBidirectional two-sample MR", p2)

# 6. Results 新小节（MR Results 之前）
results_new = """Prognostic Performance in an Independent Septic Shock Cohort

None of the cohorts analyzed above reports survival, so prognostic association was tested in GSE65682 (479 patients with septic shock; 114 deaths within 28 days, 23.8%; median age 63 years, interquartile range 53 to 71; Mars1 to Mars4 endotypes in 132, 176, 118, and 53 patients) (39). No hub gene was significantly associated with 28-day mortality on its own (P=0.136 to 0.640; Supplementary Figure S20B). In the joint six-gene model, *UPP1* was the only significant term (hazard ratio per SD 0.64; 95% CI 0.44 to 0.95; P=0.026), *GRB10* was marginal (P=0.078), and the model as a whole was not significant (likelihood-ratio P=0.067). The six-gene risk score nonetheless separated outcomes: 28-day survival was 0.70 in the high-risk half versus 0.83 in the low-risk half (log-rank P=0.002; Supplementary Figure S20A), with time-dependent AUCs of 0.638, 0.634, and 0.613 at 7, 14, and 28 days (Supplementary Figure S20C), a bootstrap-corrected C index of 0.566 (apparent 0.602; 95% CI 0.547 to 0.653), and adequate 28-day calibration (Supplementary Figure S20D). The score remained associated with mortality after adjustment for age and sex (hazard ratio per unit 2.77; 95% CI 1.57 to 4.90; P=0.0004), but added nothing to a model of age, sex, and Mars endotype (likelihood-ratio P=0.16; optimism-corrected C index 0.603 versus 0.606; net reclassification index 0.044, P=0.37) (Supplementary Figure S20E; Supplementary Table S22).

A LASSO Cox model over the 171 shared candidate genes represented on the platform told a different story. The one-standard-error rule retained no gene, making explicit that prognostic signal in this panel is weak and diffuse, and cross-validated deviance differed little between LASSO (2.085), elastic net (2.099), and ridge (2.053) penalties (Supplementary Figure S21A). At the penalty minimizing cross-validated deviance, 15 genes were retained (SNHG1, CD244, RAB13, MUC1, DEFA4, BTN3A2, TLR5, F5, NOV, LCK, ITFG2, MEGF6, EVL, FCER1A, and SGK1), none of them a hub gene (Supplementary Figure S21B; Supplementary Table S23). The penalized signature outperformed the six-gene score: optimism-corrected C index 0.622 alone and 0.638 with clinical covariates, 28-day AUC 0.698 (95% CI 0.644 to 0.753), likelihood-ratio P=9.7×10⁻¹² against the clinical model, net reclassification index 0.132 (P=0.026), and integrated discrimination improvement 0.094, with 28-day calibration close to the diagonal and positive decision-curve net benefit across the tested thresholds (Supplementary Figures S21C to S21E; Supplementary Table S22). A random survival forest over the same predictors gave no advantage (out-of-bag C index 0.610; cross-validated 28-day AUC 0.627; Supplementary Figures S21G and S21H), so the gain came from variable selection rather than from non-linear modeling. These estimates are internal to a single cohort and need confirmation in an independent cohort with a survival endpoint.

"""
t = sub_once(t, "Mendelian Randomization Does Not Detect Causal Signal (Under Limited Power)\n\nInverse-variance-weighted MR",
             results_new + "Mendelian Randomization Does Not Detect Causal Signal (Under Limited Power)\n\nInverse-variance-weighted MR", p2)

# 7. Discussion 新段（Target prioritization 之前）
disc_new = """Prognostic information is separable from the diagnostic signal

The prognostic test sharpens the marker interpretation. The six hub genes separated disease from health across cohorts but stratified 28-day mortality only modestly, and what stratification they achieved was absorbed by age, sex, and transcriptomic endotype. The prognostic content of the shared candidate set resided in 15 other genes, none of them a hub gene, and even there it was weak enough that the conservative one-standard-error selection rule retained nothing. Diagnostic and prognostic information are therefore largely separable in this gene set, consistent with markers that report the intensity of a shared host response rather than the trajectory of an individual patient. Two bounds apply. The prognostic estimates are internal to one cohort (10-fold cross-validation and bootstrap optimism correction), and the penalized signature's advantage over the six hub genes is an internally validated estimate that still requires an independent septic shock cohort with a survival endpoint; and the weak-signal reading is the safer one, since models of this kind trained and tested in the same cohort without optimism correction have been reported with considerably higher performance elsewhere (6). The six genes should therefore be read as disease-state markers; if prognostic use is pursued, it should start from the broader candidate set rather than from the hub genes.

"""
t = sub_once(t, "Target prioritization: druggability is concentrated in VNN1 and UPP1",
             disc_new + "Target prioritization: druggability is concentrated in VNN1 and UPP1", p2)

# 8. Limitations 段首句后加预后局限
t = sub_once(t,
    "These findings should be interpreted against several constraints. The discovery ARDS cohort is small",
    "These findings should be interpreted against several constraints. The prognostic analysis rests on a single cohort with an internal validation design (10-fold cross-validation and bootstrap optimism correction) and no independent prognostic cohort, so the penalized 15-gene signature requires external confirmation and its selection stability was not assessed. The discovery ARDS cohort is small", p2)

# 9. DATA AVAILABILITY
t = sub_once(t, "GSE151263, and GSE167363.", "GSE151263, GSE167363, and GSE65682.", p2)

# 10. REFERENCES 追加 39-43
t = sub_once(t,
    "38. Liu C, Yang Y, Song Y, Wu W, He P, Han Y, et al. Vascular non-inflammatory molecule-1 aggravates acute kidney injury induced acute lung injury via promoting neutrophil extracellular traps formation. *Inflammation*. 2026;49(1):60. doi:10.1007/s10753-026-02457-y.",
    "38. Liu C, Yang Y, Song Y, Wu W, He P, Han Y, et al. Vascular non-inflammatory molecule-1 aggravates acute kidney injury induced acute lung injury via promoting neutrophil extracellular traps formation. *Inflammation*. 2026;49(1):60. doi:10.1007/s10753-026-02457-y.\n39. Davenport EE, Burnham KL, Radhakrishnan J, Humburg P, Hutton P, Mills TC, et al. Genomic landscape of the individual host response and outcomes in sepsis: a prospective cohort study. *Lancet Respir Med*. 2016;4(4):259–271. doi:10.1016/S2213-2600(16)00046-1.\n40. Blanche P, Dartigues JF, Jacqmin-Gadda H. Estimating and comparing time-dependent areas under receiver operating characteristic curves for censored event times with competing risks. *Stat Med*. 2013;32(30):5381–5397. doi:10.1002/sim.5958.\n41. Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized linear models via coordinate descent. *J Stat Softw*. 2010;33(1):1–22. doi:10.18637/jss.v033.i01.\n42. Wright MN, Ziegler A. ranger: A fast implementation of random forests for high dimensional data in C++ and R. *J Stat Softw*. 2017;77(1):1–17. doi:10.18637/jss.v077.i01.\n43. Pencina MJ, D'Agostino RB Sr, D'Agostino RB Jr, Vasan RS. Evaluating the added predictive ability of a new marker: from area under the ROC curve to reclassification and beyond. *Stat Med*. 2008;27(2):157–172. doi:10.1002/sim.2929.", p2)

io.open(p2, "w", encoding="utf-8").write(t)
print("02_Manuscript_Revised.md 更新完成")

# ============================================================
# 中文直译稿 03_Manuscript_Translation_CN.md
# ============================================================
p3 = f"{PKG}/03_Manuscript_Translation_CN.md"
t = io.open(p3, encoding="utf-8").read()

t = sub_once(t,
    "随后通过免疫去卷积、单细胞转录组学（含双重虚拟敲除）、药物重定位预测与分子对接、以及孟德尔随机化对结果进行拓展。",
    "随后通过免疫去卷积、单细胞转录组学（含双重虚拟敲除）、药物重定位预测与分子对接、孟德尔随机化、以及在一个独立脓毒性休克队列中的预后评价对结果进行拓展。", p3)

t = sub_once(t,
    "而反向分析具有信息量（最小可检测比值比为1.04）。",
    "而反向分析具有信息量（最小可检测比值比为1.04）。在预后队列（479例脓毒性休克患者，28天内死亡114例）中，6基因评分对死亡有分层能力（log-rank P=0.002），但在年龄、性别与内型之上没有增量（似然比P=0.16）；而从179个共享候选基因中抽取的惩罚回归标签保留了预后信息（乐观校正C指数为0.622），且不含6个hub基因中的任何一个。", p3)

t = sub_once(t,
    "它们是尚需实验验证的候选疾病状态标志物。",
    "它们是候选疾病状态标志物而非预后标志物，且尚需实验验证。", p3)

t = sub_once(t,
    "其中1例有两个时间点样本；12359个细胞）。",
    "其中1例有两个时间点样本；12359个细胞）。预后评价使用GSE65682（GPL13667；脓毒症分子诊断与风险分层研究[MARS]479例脓毒性休克患者的全血，于重症监护病房入院24小时内采样），这是唯一具有生存终点的队列（28天死亡；114例死亡，23.8%），含年龄、性别与Mars1至Mars4内型分组（39）。", p3)

methods_new_cn = """### 独立脓毒性休克队列的预后建模

预后分析仅限于GSE65682——唯一具有生存终点的队列。表达值取自处理后的系列矩阵；探针按每基因平均表达最高者折叠为基因符号，179个共享候选基因中有171个在该平台上有代表。终点为28天死亡；随访时间为零的患者记为0.5天。6个hub基因的Cox模型以z分数标化表达分别（单因素）与联合（多因素）拟合；多因素模型的线性预测值作为风险评分，患者按其中位数分组绘制Kaplan-Meier曲线（log-rank检验）。时间依赖判别在7、14与28天用删失逆概率加权受试者工作特征曲线评估（40），28天评价置于27.99天以规避行政删失边界。为检验6个hub基因之外是否存在预后信息，将临床模型（年龄、性别与Mars内型）分别与加入6基因评分的模型、以及基于171个共享候选基因的惩罚Cox模型（LASSO，glmnet，惩罚系数经10折交叉验证选择）（41）比较；由于一倍标准误规则未保留任何基因，采用交叉验证偏似然偏差最小对应的惩罚系数（保留15个基因），并拟合弹性网与岭回归作为对照。随机生存森林（ranger，500棵树）以相同基因加临床变量作为非线性对照，个体水平风险取10折交叉验证的折外预测——因为该样本量下in-bag预测严重乐观（42）。判别以Harrell C指数汇总，并用300次bootstrap重采样做乐观校正；每次重采样内LASSO惩罚系数固定于全队列选定值，因此惩罚模型的乐观估计若有偏倚也是偏保守的。校准以5个风险组中平均预测与Kaplan-Meier观察的28天死亡率比较评估。相对临床模型的增量价值用似然比检验检验，并以28天净重分类指数与综合判别改善量化（43）；临床实用性用决策曲线分析（阈值概率0.05至0.60）评估。

"""
t = sub_once(t, "### 孟德尔随机化\n\n", methods_new_cn + "### 孟德尔随机化\n\n", p3)

results_new_cn = """### 独立脓毒性休克队列的预后表现

上述队列均不报告生存，故预后关联在GSE65682中检验（479例脓毒性休克患者；28天内死亡114例，23.8%；中位年龄63岁，四分位距53至71；Mars1至Mars4内型分别为132、176、118与53例）（39）。没有任何hub基因单独与28天死亡显著相关（P=0.136至0.640；补充图S20B）。在6基因联合模型中，*UPP1*是唯一显著项（每SD的风险比0.64；95% CI 0.44至0.95；P=0.026），*GRB10*处于边缘（P=0.078），模型整体不显著（似然比P=0.067）。6基因风险评分仍能区分结局：高风险半组的28天生存率为0.70，低风险半组为0.83（log-rank P=0.002；补充图S20A），7、14与28天的时间依赖AUC分别为0.638、0.634与0.613（补充图S20C），bootstrap校正C指数0.566（表观0.602；95% CI 0.547至0.653），28天校准良好（补充图S20D）。该评分在调整年龄与性别后仍与死亡相关（每单位风险比2.77；95% CI 1.57至4.90；P=0.0004），但在年龄、性别与Mars内型的模型之上没有增量（似然比P=0.16；乐观校正C指数0.603对0.606；净重分类指数0.044，P=0.37）（补充图S20E；补充表S22）。

基于该平台上有代表的171个共享候选基因的LASSO Cox模型给出了不同结果。一倍标准误规则未保留任何基因，明确表明该面板的预后信号弱且弥散；LASSO（2.085）、弹性网（2.099）与岭回归（2.053）的交叉验证偏似然偏差差异很小（补充图S21A）。在交叉验证偏差最小的惩罚系数下保留了15个基因（SNHG1、CD244、RAB13、MUC1、DEFA4、BTN3A2、TLR5、F5、NOV、LCK、ITFG2、MEGF6、EVL、FCER1A与SGK1），其中没有一个是hub基因（补充图S21B；补充表S23）。该惩罚标签优于6基因评分：单独的乐观校正C指数为0.622，联合临床协变量后为0.638；28天AUC为0.698（95% CI 0.644至0.753）；相对临床模型的似然比检验P=9.7×10⁻¹²；净重分类指数0.132（P=0.026），综合判别改善0.094；28天校准贴合对角线，决策曲线在所检验阈值范围内净获益为正（补充图S21C至S21E；补充表S22）。基于相同预测变量的随机生存森林没有优势（out-of-bag C指数0.610；交叉验证28天AUC 0.627；补充图S21G与S21H），说明增益来自变量选择而非非线性建模。这些估计均为单队列内部估计，需在有生存终点的独立队列中确认。

"""
t = sub_once(t, "### 孟德尔随机化未检测到因果信号（效能有限）\n\n",
             results_new_cn + "### 孟德尔随机化未检测到因果信号（效能有限）\n\n", p3)

disc_new_cn = """预后信息与诊断信号可分离

预后检验进一步锐化了标志物定位。6个hub基因在各队列中能区分疾病与健康，但对28天死亡仅有中等的分层能力，且其分层作用可被年龄、性别与转录组内型吸收。共享候选基因集的预后内容存在于另外15个基因中，其中没有一个是hub基因；即便如此，其信号弱到保守的一倍标准误选择规则未保留任何基因。因此，在该基因集中诊断信息与预后信息大体可分离，这与作为共享宿主反应强度而非个体患者病程轨迹的标志物定位一致。有两条边界。预后估计为单队列内部估计（10折交叉验证与bootstrap乐观校正），惩罚标签相对6个hub基因的优势是内部验证的估计，仍需具有生存终点的独立脓毒性休克队列确认；且弱信号解读更为稳妥，因为这类模型在同一队列内训练与检验且不做乐观校正时，在别处曾被报告过明显更高的性能（6）。因此这6个基因应解读为疾病状态标志物；若开展预后用途，应从更宽的候选基因集出发而非从hub基因出发。

"""
t = sub_once(t, "靶点优先级：可成药性集中于VNN1与UPP1",
             disc_new_cn + "靶点优先级：可成药性集中于VNN1与UPP1", p3)

t = sub_once(t,
    "这些结果应在若干约束下解读。ARDS发现队列较小",
    "这些结果应在若干约束下解读。预后分析基于单队列的内部验证设计（10折交叉验证与bootstrap乐观校正）且无独立预后队列，因此惩罚15基因标签需要外部确认，其选择稳定性未评估。ARDS发现队列较小", p3)

t = sub_once(t,
    "登录号为GSE76293、GSE69528、GSE100159、GSE243066、GSE151263和GSE167363。",
    "登录号为GSE76293、GSE69528、GSE100159、GSE243066、GSE151263、GSE167363和GSE65682。", p3)

io.open(p3, "w", encoding="utf-8").write(t)
print("03_Manuscript_Translation_CN.md 更新完成")

# ============================================================
# Supplementary_Materials.md
# ============================================================
p4 = f"{PKG}/Supplementary_Materials.md"
t = io.open(p4, encoding="utf-8").read()

table_new = """### Table S22. Prognostic model comparison in GSE65682

Discrimination and incremental-value metrics for the prognostic models fitted in GSE65682 (479 patients with septic shock; 28-day mortality 23.8%): apparent and bootstrap-optimism-corrected C index (300 resamples), time-dependent AUC at 7, 14, and 28 days, likelihood-ratio tests against the clinical model (age, sex, and Mars endotype), and net reclassification index and integrated discrimination improvement at 28 days. Models: Clinical (age, sex, endotype); Hub6 (six hub genes, z-standardized); Clinical + Hub6; LASSO-179 (penalized Cox signature over the 171 shared candidate genes represented on the platform; penalty at minimum 10-fold cross-validated deviance, 15 genes retained); Clinical + LASSO-179; and a random survival forest over the same predictors (out-of-bag and 10-fold out-of-fold estimates).
*Source: `13_survival_prognosis/tables/ext_model_performance.csv`, `ext_bootstrap_optimism.csv`, `ext_lrt.csv`, `ext_nri_idi.csv`; scripts `scripts/step105_survival_prognosis_GSE65682.R`, `scripts/step106_extended_modeling_GSE65682.R`*

### Table S23. LASSO-179 prognostic signature coefficients

Coefficients of the 15 genes retained by the LASSO Cox model at the penalty minimizing 10-fold cross-validated deviance in GSE65682, with hazard ratios per unit of log2 expression; none of the six hub genes was retained, and the more conservative one-standard-error rule retained no gene at all.
*Source: `13_survival_prognosis/tables/ext_lasso_coefficients.csv`; script `scripts/step106_extended_modeling_GSE65682.R`*

"""
t = sub_once(t, "---\n\n## Supplementary Figure Legends",
             table_new + "---\n\n## Supplementary Figure Legends", p4)

fig_new = """### Figure S20. Six-hub-gene prognostic performance in the independent septic shock cohort GSE65682

479 patients with septic shock (MARS cohort; 114 deaths within 28 days). (A) Kaplan-Meier curves for patients grouped at the median of the six-hub-gene risk score (28-day survival 0.83 low-risk versus 0.70 high-risk; log-rank P=0.002; shaded bands, 95% confidence bands). (B) Kaplan-Meier curves for each hub gene individually (median split; log-rank P shown per panel); none of the six genes stratified mortality significantly. (C) Time-dependent receiver-operating-characteristic curves of the risk score at 7, 14, and 28 days (inverse-probability-of-censoring weighting; the 28-day evaluation was placed at 27.99 days to avoid the administrative censoring boundary). (D) 28-day calibration of the six-gene Cox model (rms calibrate, hare method, 400 bootstrap resamples); the bias-corrected curve follows the diagonal. (E) Multivariable Cox model adjusted for age and sex: risk-score hazard ratio 2.77 per unit (95% CI 1.57 to 4.90; P=0.0004), age hazard ratio 1.02 per year (P=0.004), and sex not significant (P=0.327).
*Source: `13_survival_prognosis/figures/FigureS20_sixgene_prognosis_GSE65682.{pdf,png,tiff}`; scripts `scripts/step105_survival_prognosis_GSE65682.R`, `scripts/step107_sixgene_prognosis_figure.R`*

### Figure S21. Extended prognostic modeling over the 179 shared candidate genes in GSE65682

(A) Ten-fold cross-validated partial-likelihood deviance of the LASSO Cox model over the 171 shared candidate genes represented on the platform; the penalty at minimum deviance (lambda.min, 15 genes retained) and the one-standard-error penalty (lambda.1se, no gene retained) are marked. (B) Coefficients of the 15 genes retained at lambda.min (none of the six hub genes was retained). (C) Time-dependent receiver-operating-characteristic curves at 28 days for the clinical model (age, sex, Mars endotype), the six-hub-gene model, their combination, the LASSO signature alone and with clinical covariates, and the random survival forest (out-of-fold predictions). (D) Five-group calibration of predicted versus Kaplan-Meier-observed 28-day mortality for the same models. (E) Decision-curve analysis (net benefit versus threshold probability, 0.05 to 0.60). (F) Incremental value over the clinical model: change in bootstrap-optimism-corrected C index, net reclassification index, and integrated discrimination improvement (x100) at 28 days. (G) Corrected C-index comparison including the random survival forest (out-of-bag). (H) Random-survival-forest variable importance (impurity) for the 12 most important predictors; age dominates, and no hub gene appears.
*Source: `13_survival_prognosis/figures/FigureS21_extended_prognostic_modeling.{pdf,png,tiff}`; script `scripts/step106_extended_modeling_GSE65682.R`*

"""
t = sub_once(t, "## File Manifest",
             fig_new + "## File Manifest", p4)

t = sub_once(t, "| `Supplementary_Materials.md` | This index (Tables S1-S21 and Figure legends S1-S19) |",
             "| `Supplementary_Materials.md` | This index (Tables S1-S23 and Figure legends S1-S21) |", p4)
t = sub_once(t, "| `figures/supplementary/FigureS1-S19.*` | Supplementary figures (renumbered), PDF/PNG/TIFF triple export |",
             "| `figures/supplementary/FigureS1-S21.*` | Supplementary figures (renumbered), PDF/PNG/TIFF triple export |", p4)
t = sub_once(t, "| `figures/supplementary_tables/TableS1-S21.xlsx` | Supplementary tables (renumbered, multi-sheet workbooks) |",
             "| `figures/supplementary_tables/TableS1-S23.*` | Supplementary tables (renumbered; S22 and S23 as CSV) |", p4)

io.open(p4, "w", encoding="utf-8").write(t)
print("Supplementary_Materials.md 更新完成")

# ============================================================
# Supplementary_Methods.md：新增 S14
# ============================================================
p5 = f"{PKG}/Supplementary_Methods.md"
t = io.open(p5, encoding="utf-8").read()
s14 = """
## S14. Prognostic modeling in GSE65682

Prognostic analysis was confined to GSE65682 (GPL13667, Affymetrix HG-U219), the only cohort among those analyzed with a survival endpoint. The processed series matrix was used as provided; probe-to-gene mapping used the platform annotation extracted from the family SOFT file, and probes were collapsed to gene symbols by the probe with the highest mean expression across samples. Of the 179 shared candidate genes, 171 were represented; four genes were not annotated on the platform (ZBTB42, PIK3R6, ZNF329, CCL23), and four additional mapped genes had no measurable expression and were dropped. The analysis set comprised the 479 patients with septic shock with a 28-day endpoint, after removing duplicate draws from the same patient (first sample retained). The endpoint was 28-day mortality (114 events, 23.8%); a recorded follow-up time of zero was recoded as 0.5 days. Time-dependent ROC used the inverse-probability-of-censoring weighting estimator (40), with the 28-day evaluation placed at 27.99 days to avoid the administrative censoring boundary. Bootstrap optimism correction used 300 resamples; the LASSO penalty was held fixed at the full-cohort value within each resample, so the optimism estimate for the penalized models is, if anything, slightly conservative. For the net reclassification index and integrated discrimination improvement at 28 days (43), events were deaths within 28 days (n=114), non-events were patients alive and still under follow-up at 28 days (n=363; the two patients whose last follow-up preceded 28 days without an event were excluded from this calculation), and the classification cutoff was the cohort event rate (0.239). The random survival forest used ranger (500 trees, minimum node size 15, mtry set to the square root of the number of predictors) (42); individual-level risks were taken from the out-of-fold predictions of a 10-fold cross-validation because in-bag predictions were strongly optimistic at this sample size, and the forest's out-of-bag C index was reported alongside.
"""
t = t.rstrip() + "\n" + s14 + "\n"
io.open(p5, "w", encoding="utf-8").write(t)
print("Supplementary_Methods.md 更新完成")

# ============================================================
# 04_Modification_Log.md：新增小节
# ============================================================
p6 = f"{PKG}/04_Modification_Log.md"
t = io.open(p6, encoding="utf-8").read()
log_new = """

## 2026-10-15：GSE65682 生存/预后分析与扩展建模（整合）

### 触发
在四个数据集上完善生存分析、预后模型构建与评价。核查后仅 GSE65682（MARS 脓毒性休克队列，479 例，28 天死亡 114 例）含生存终点；GSE76293（体外 PMN 实验）、GSE100159（类鼻疽）、GSE243066（ARDS vs 健康）均无生存字段。

### 分析（13_survival_prognosis/）
| 脚本 | 内容 |
|---|---|
| step105_survival_prognosis_GSE65682.R | 6 hub 基因 Cox、风险评分 KM、timeROC、C 指数、校准、nomogram、调整后森林图 |
| step105b_pergene_KM_GSE65682.R | 单基因 KM |
| step106_extended_modeling_GSE65682.R | 171 基因 LASSO-Cox（10 折 CV）、弹性网/岭回归对照、随机生存森林、bootstrap 乐观校正（B=300）、NRI/IDI、DCA、校准比较 |
| step107_sixgene_prognosis_figure.R | Figure S20 组装 |

### 主要结果
- 六基因：单基因均不显著；多因素 Cox 中 UPP1 显著（HR 0.64, P=0.026）；风险评分 KM P=0.002（28 天生存 0.70 vs 0.83）；校正 C 指数 0.566；在年龄+性别+endotype 模型上无增量（LRT P=0.161）。
- 扩展建模：LASSO 1-SE 规则未选基因（信号弱），lambda.min 保留 15 基因（无 hub 基因）；校正 C 指数 0.622（联合临床 0.638）；28 天 AUC 0.698；NRI 0.132（P=0.026）、IDI 0.094；RSF 无增益（OOB C 0.610）。
- 解读：诊断标志物与预后信息可分离；hub 基因为疾病状态标志物而非预后标志物。

### 手稿修改（英文+中文同步）
- Abstract：Methods 句加预后评价；M&M 加预后结果句；Conclusions 改为"疾病状态标志物而非预后标志物"。
- Methods：Data Sources 加 GSE65682；新增小节"Prognostic Modeling in an Independent Septic Shock Cohort"。
- Results：新增小节"Prognostic Performance in an Independent Septic Shock Cohort"（两段）。
- Discussion：新增段落"Prognostic information is separable from the diagnostic signal"；limitations 段加预后局限句。
- DATA AVAILABILITY：加 GSE65682。
- REFERENCES：新增 39-43（Davenport/Blanche/Friedman/Wright/Pencina）。
- Supplementary_Materials.md：新增 Table S22、S23 条目与 Figure S20、S21 图注；File Manifest 更新。
- Supplementary_Methods.md：新增 S14（GSE65682 预后建模）。
- 新增补充图：FigureS20_sixgene_prognosis_GSE65682、FigureS21_extended_prognostic_modeling（PDF/PNG/TIFF 300dpi）。
- 新增补充表：TableS22_prognostic_model_comparison.csv、TableS23_lasso_prognostic_signature.csv。
"""
t = t.rstrip() + "\n" + log_new + "\n"
io.open(p6, "w", encoding="utf-8").write(t)
print("04_Modification_Log.md 更新完成")

print("=== 全部整合完成 ===")
