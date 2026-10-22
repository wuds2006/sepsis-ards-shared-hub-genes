#!/usr/bin/env python3
# step100: rewrite Introduction per 7-step framework (anti-defensive), merge Q+approach paragraphs,
# sync CN translation, regenerate docx, append modification log.
import re, subprocess, sys

BASE = '/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission/09_review_deAI_2026-09-25'
EN = f'{BASE}/02_Manuscript_Revised.md'
CN = f'{BASE}/03_Manuscript_Translation_CN.md'
LOG = f'{BASE}/04_Modification_Log.md'

EN_INTRO = """Sepsis, life-threatening organ dysfunction caused by a dysregulated host response to infection (1), remains a leading cause of death in the critically ill. The lung is the organ it injures most often, and sepsis is the most common predisposing condition for the acute respiratory distress syndrome (ARDS) (2), a heterogeneous syndrome of acute hypoxemic respiratory failure marked by increased alveolar\u2013capillary permeability and neutrophil-rich inflammation (3). Despite lung-protective ventilation, mortality of sepsis-associated ARDS remains high and no targeted molecular therapy exists (4). The two syndromes share much of their inflammatory machinery, including complement activation, cytokine release, and engagement of neutrophils and monocytes/macrophages (5), yet they have been studied largely as separate entities.

Considerable groundwork exists. Shared-gene screens spanning the two diseases have been reported (6\u20138), and the clinical immunology of each condition is well characterized in isolation: sepsis induces a biphasic host response in which early cytokine storm is followed by protracted immunosuppression marked by lymphocyte apoptosis, expansion of myeloid-derived suppressor cells, and downregulation of monocyte HLA-DR (9\u201312), whereas ARDS carries neutrophil-rich, complement-activating lung inflammation (2, 5). Independent work has converged on individual sepsis-relevant molecules: the dynamic neutrophil CD177/CD10 ratio differentiates infection from non-infectious inflammation and stratifies mortality in critically ill patients (13), CD177+ neutrophils have been linked to necroptosis-driven inflammation and adverse outcomes in sepsis (14), CD177 expression defines functionally stable neutrophil subtypes (15), and macrophage UPP1 drives a mitochondrial-ROS\u2013cGAS\u2013NLRP3 inflammatory axis (16). Computational drug repurposing has likewise been attempted, with Connectivity Map (CMap) queries against sepsis and ARDS expression signatures nominating candidate compounds, including a recent septic-ARDS analysis that prioritized etoposide from 49 key genes (17).

Three limitations persist across this literature. First, shared-gene screens were built on processed matrices without systematic removal of off-target samples: GSE76293 mixes disease samples with in vitro-cultured neutrophils treated with PI3K inhibitor or GM-CSF, and GSE69528 includes uninfected patients with type 2 diabetes, contaminants that confound differential expression and machine-learning inference. Second, candidates were nominated and validated within the same data, so the shared molecular intersection has never been tested under discovery-strict external validation. Third, genetic causality and cell-level pseudoreplication remain unaddressed, and signature-matching repurposing nominations inherit well-documented reproducibility limits of reference compendia profiled in immortalized cancer cell lines (18, 19). A biological question is also unresolved: whether the immune changes of sepsis-associated ARDS are simply the pulmonary continuation of the septic response or diverge from it in organ-specific ways has rarely been examined at the level of coordinately regulated genes.

This study therefore asks a specific question: which genes and immune processes are coordinately dysregulated in both sepsis and ARDS and remain supported after raw-data re-curation, external validation in independent cohorts, and orthogonal causal and single-cell testing? To answer it, we re-analyzed both discovery cohorts from raw probe-level matrices with off-target samples removed. Differentially expressed genes intersected with WGCNA disease-associated modules yielded 179 shared candidate genes; 11 machine-learning classifiers with SHAP interpretation, trained separately in each disease, narrowed these to candidate hub genes, and cross-disease candidacy defined within the discovery data alone yielded six hub genes (*CDK5RAP2*, *CD177*, *GRB10*, *UPP1*, *EXOSC4*, *VNN1*), which were then tested in two independent validation cohorts used strictly for confirmation. The working premise was that genuinely shared hub genes should surface consistently across cohorts, cell types, and analytic angles. We tested this premise along four orthogonal lines: immune deconvolution and single-cell mapping in two cohorts (cell-cell communication, pseudotime, AUCell gene-set activity, and dual virtual knockout by scTenifoldKnk and CellOracle); pathway-level enrichment; genetic causal inference by Mendelian randomization with SMR-HEIDI, linkage-disequilibrium score regression, and colocalization; and in silico drug repurposing prediction around the six hub genes, with CMap specificity quantified against random gene-set nulls and every docking claim checked by redocking, predicted-structure confidence, and free-energy rescoring.

The result is a stringently curated molecular intersection of sepsis and ARDS: six hub genes whose shared dysregulation survives re-curation and independent validation; an immune and regulatory map showing where and how these genes act in defined myeloid populations; genetic evidence on their causal status; and a bounded set of candidate compounds for experimental follow-up. For CD177 and UPP1, the analysis extends two well-supported sepsis markers into the shared sepsis\u2013ARDS context and adds the single-cell regulatory and genetic-causal depth that marker-only screens cannot supply. The curation-first, discovery-strict design is directly reusable for other comorbid disease pairs in public transcriptomic resources."""

CN_INTRO = """脓毒症是感染所致宿主反应失调引起的危及生命的器官功能障碍（1），仍是危重症患者的主要死因。肺是脓毒症最常损伤的器官，脓毒症也是急性呼吸窘迫综合征（ARDS）最常见的易感条件（2）；ARDS是一种以肺泡-毛细血管通透性增高和中性粒细胞富集性炎症为特征的急性低氧性呼吸衰竭异质性综合征（3）。尽管采用肺保护性通气，脓毒症相关ARDS的病死率仍然很高，且尚无靶向分子治疗（4）。两种综合征共享大部分炎症机制，包括补体激活、细胞因子释放以及中性粒细胞和单核/巨噬细胞的参与（5），但长期以来基本被作为独立实体研究。

前人已打下相当的基础。跨越两种疾病的共享基因筛选已有报道（6-8），两种疾病各自的临床免疫学也已在孤立状态下被充分刻画：脓毒症诱导双相宿主反应，早期细胞因子风暴之后是迁延的免疫抑制期，其特征为淋巴细胞凋亡、髓源性抑制细胞扩增和单核细胞HLA-DR下调（9-12）；ARDS则表现为中性粒细胞富集、补体激活的肺部炎症（2，5）。独立研究已汇聚到若干脓毒症相关分子上：动态中性粒细胞CD177/CD10比值可区分感染与非感染性炎症，并对危重症患者的死亡风险进行分层（13）；CD177+中性粒细胞与脓毒症中坏死性凋亡驱动的炎症及不良结局相关（14）；CD177表达定义了功能稳定的中性粒细胞亚型（15）；巨噬细胞UPP1驱动线粒体ROS-cGAS-NLRP3炎症轴（16）。计算药物重定位同样有尝试：针对脓毒症与ARDS表达标签的Connectivity Map（CMap）查询提名了候选化合物，包括近期一项从49个关键基因出发优先提名依托泊苷的脓毒性ARDS分析（17）。

这一文献体系仍存在三方面局限。第一，共享基因筛选建立在处理后矩阵之上，未系统剔除脱靶样本：GSE76293将疾病样本与PI3K抑制剂或GM-CSF处理的体外培养中性粒细胞混在一起，GSE69528包含未感染的2型糖尿病患者；这些污染样本会混淆差异表达与机器学习推断。第二，候选基因在同一数据内提名并验证，共享分子交集从未在严格限定发现集的外部验证下接受检验。第三，遗传因果性与细胞水平假重复问题仍未处理，标签匹配式重定位提名还继承了参考数据集在永生化癌细胞系中测定所带来的有据可查的可重复性局限（18，19）。一个生物学问题同样悬而未决：脓毒症相关ARDS的免疫改变究竟只是脓毒症反应在肺部的延续，还是以器官特异性方式发生偏离，这一问题很少在协同调控基因的层面上被检验。

因此，本研究提出一个具体问题：哪些基因与免疫过程在脓毒症与ARDS中协同失调，并在原始数据再整编、独立队列外部验证以及正交的因果与单细胞检验之后仍获支持？为回答这一问题，我们从原始探针水平矩阵重新分析两个发现队列，并剔除脱靶样本。差异表达基因与WGCNA疾病相关模块取交集得到179个共享候选基因；分别在每种疾病中训练的11种机器学习分类器结合SHAP解释将其缩减为候选hub基因；跨疾病候选资格仅在发现数据内定义，最终得到6个hub基因（*CDK5RAP2*、*CD177*、*GRB10*、*UPP1*、*EXOSC4*、*VNN1*），随后在两个严格用于确认的独立验证队列中检验。我们的工作前提是：真正共享的hub基因应当在不同队列、细胞类型和分析角度下一致地浮现。我们沿四条正交证据线检验这一前提：在两个队列中进行免疫去卷积与单细胞图谱绘制（细胞间通讯、拟时序、AUCell基因集活性，以及scTenifoldKnk与CellOracle双重虚拟敲除）；通路水平富集分析；通过孟德尔随机化结合SMR-HEIDI、连锁不平衡评分回归与共定位进行遗传因果推断；围绕6个hub基因开展计算机模拟药物重定位预测，其中CMap特异性以随机基因集零分布量化，每一项对接论断均经回对接、预测结构置信度和自由能重评分核查。

最终得到的是一个经严格整编的脓毒症与ARDS分子交集：6个共享失调通过再整编与独立验证留存的hub基因；一幅展示这些基因在特定髓系细胞群中作用位置与方式的免疫与调控图谱；关于其因果地位的遗传学证据；以及一组边界限定、可供实验随访的候选化合物。对CD177与UPP1而言，本分析把两个证据充分的脓毒症标志物拓展到脓毒症-ARDS共享背景中，并补充了仅靠标志物筛选无法提供的单细胞调控与遗传因果深度。整编优先、严格限定发现集的设计可直接复用于公共转录组资源中的其他共病疾病对。"""

# --- 1. English manuscript ---
txt = open(EN, encoding='utf-8').read()
m = re.search(r'INTRODUCTION\n\n(.*?)\n\n---\n\nMETHODS', txt, re.S)
assert m, 'EN Introduction block not found'
old = m.group(1)
txt = txt[:m.start(1)] + EN_INTRO + txt[m.end(1):]
open(EN, 'w', encoding='utf-8').write(txt)
print(f'EN Introduction replaced: {old.count(chr(10)+chr(10))+1} -> 5 paragraphs')

# --- 2. CN translation ---
txt = open(CN, encoding='utf-8').read()
m = re.search(r'## 引言\n\n(.*?)\n\n---\n\n## 方法', txt, re.S)
assert m, 'CN Introduction block not found'
txt = txt[:m.start(1)] + CN_INTRO + txt[m.end(1):]
open(CN, 'w', encoding='utf-8').write(txt)
print('CN Introduction replaced -> 5 paragraphs')

# --- 3. Modification log entry ---
entry = """

---

## 2026-10-07：Introduction 按 7 步框架重写（去防御性表达）

按用户指定的 7 步框架（研究背景 → 研究现状 → research gap → 研究问题 → commitment/论证思路合并段 → 贡献）将 Introduction 由 7 段重写为 5 段；应用户要求将"研究问题"与"commitment 与论证思路"合并为第 4 段。引用 1–19 全部保留且首引顺序连续，无需重排；无新增事实。

防御性表达处理：
- 删除 "the contribution lies not in nominating these genes de novo, but in extending..." 自贬结构，改为正向陈述（末段 "the analysis extends ... and adds ... depth"）。
- 删除 "framed, from the outset, as in silico drug repurposing prediction rather than drug repositioning in the validated sense" 预防性划界，改为正向描述所做工作及对照。
- "we bounded every downstream claim ..." 总结式对冲融入四条正交证据线的具体描述。
- "Part of the reason may be simple" 弱断言改为直接陈述（两病共享炎症机制却分开研究）。
- "to our knowledge without systematic off-target sample removal" 不确定措辞改为肯定句。
- "has rarely been asked" 明确为 "has rarely been examined at the level of coordinately regulated genes"。

同步更新：02_Manuscript_Revised.md、03_Manuscript_Translation_CN.md、manuscript_revised.docx（step95 重生成）。脚本：scripts/step100_update_introduction.py。
"""
open(LOG, 'a', encoding='utf-8').write(entry)
print('Modification log appended')

# --- 4. Regenerate docx ---
r = subprocess.run([sys.executable, '/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/step95_make_revised_docx.py'],
                   capture_output=True, text=True)
print(r.stdout[-500:] if r.stdout else '')
if r.returncode != 0:
    print('STEP95 FAILED:', r.stderr[-1000:])
    sys.exit(1)
print('docx regenerated')
