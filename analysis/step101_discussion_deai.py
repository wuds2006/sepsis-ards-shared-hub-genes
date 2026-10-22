#!/usr/bin/env python3
# step101: apply 14 confirmed anti-defensive edits to DISCUSSION in EN revised manuscript,
# mirror them in CN translation, regenerate docx, append modification log.
import re, subprocess, sys

BASE = '/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission/09_review_deAI_2026-09-25'
EN = f'{BASE}/02_Manuscript_Revised.md'
CN = f'{BASE}/03_Manuscript_Translation_CN.md'
LOG = f'{BASE}/04_Modification_Log.md'

# (old, new) literal replacements for EN; each must occur exactly once
EN_EDITS = [
    # #1 P1 stacked caveat -> single bounded statement
    ("supports a mixed reading: partial extension, with an organ-level imprint, although that imprint rests on a single-cohort cell-cell-communication ranking and two-to-four-donor-per-group comparisons and cannot adjudicate the question.",
     "supports a mixed reading of partial extension with an organ-level imprint, pending confirmation in larger multi-cohort data."),
    # #11 P1 preemptive terminology disclaimer -> delete
    (" Throughout this report, \"hub gene\" is used in its operational, selection-rule sense (concordant nomination by LASSO selection, random-forest ranking, and SHAP importance within discovery data), not to imply degree centrality in a protein-interaction network.",
     ""),
    # #2a P2 apology clause -> delete
    ("The six genes are biologically connected rather than a random list, an inference drawn from existing literature rather than from new experimental evidence generated here.",
     "The six genes are biologically connected rather than a random list."),
    # #2b/#6 P2 double disclaimer -> single plain statement
    ("The GO-derived groupings of these genes (pyrimidine metabolism, RNA exosome, pantetheinase signaling) should be read as annotation-based categories rather than experimentally demonstrated co-regulated modules; nothing in the present data establishes coordinate regulation across those categories, and the classification is a hypothesis-generating framework rather than a confirmed functional architecture.",
     "These GO-derived groupings (pyrimidine metabolism, RNA exosome, pantetheinase signaling) are annotation-based categories; whether they function as co-regulated modules requires experimental testing."),
    # #3 P3 donor-count caveat repeated -> delete (kept in limitations)
    (", a hypothesis that rests on two to four donors per group and requires confirmation in larger cohorts.",
     "."),
    # #7 P3 stacked hedge -> single calibrated statement
    ("the inverted-U shape should be read with its cross-cohort provenance in mind (see above), and its apparent alignment with the hyperinflammation-to-immunosuppression transition of clinical sepsis (9, 10) is consistent with that transition although it cannot be taken as evidence for it, whereas",
     "the inverted-U shape mirrors the hyperinflammation-to-immunosuppression transition of clinical sepsis (9, 10), with the cross-cohort provenance limiting it to a concordant pattern, whereas"),
    # #12 P3 weak hedge
    ("is presumably carried by its mononuclear-cell components:",
     "is carried by its mononuclear-cell components:"),
    # #8 P3 directive framing -> direct statement
    ("The transferred model should therefore be read as a myeloid-compartment signature rather than a neutrophil readout in PBMCs.",
     "The transferred model is therefore a myeloid-compartment signature, not a neutrophil readout, in PBMCs."),
    # #4b P4 forward-null caveat second restatement -> compress
    ("under such weak genetic architecture the forward null cannot exclude small or moderate causal effects (power analysis in Results), so any modest true effect would have gone undetected.",
     "under such weak genetic architecture any modest true effect would have gone undetected (power analysis in Results)."),
    # #4c/#9 P4 third restatement -> delete
    ("low power. The null should therefore be read as absence of detectable causal signal rather than proof of absence of causality. Protein-level",
     "low power. Protein-level"),
    # #13a P4 defensive tail
    ("a hypothesis that requires experimental validation rather than a prioritized target claim.",
     "a hypothesis requiring experimental validation."),
    # #10 P5 directive framing
    ("so the CDK5RAP2 docking energies should be read as exploratory (Supplementary Table S5)",
     "so the CDK5RAP2 docking energies are exploratory (Supplementary Table S5)"),
    # #5 P7 circularity third time -> delete
    ("claims. The CMap glucocorticoid signal remains partly circular in expectation, although the random gene-set nulls show that its top-ten concentration exceeds that of same-background random sets (a glucocorticoid entered the top ten for the real signature versus 17% of same-background and 0% of transcriptome-wide nulls; see Target prioritization). The sepsis nomogram",
     "claims. The sepsis nomogram"),
    # #14 P7 generic closing disclaimer -> delete
    ("result. Finally, everything here is computational and awaits prospective or experimental validation.",
     "result."),
    # #13b P8 defensive tail -> positive close
    ("and computational screening nominates VNN1 as a candidate for experimental validation rather than an established target.",
     "and computational screening nominates VNN1 as the lead candidate for that validation."),
]

# CN edits: (pattern, new, is_regex)
CN_EDITS = [
    ("支持一种混合解读：部分延伸，叠加器官水平印记；但该印记仅基于单队列的细胞间通讯排序和每组2至4个供体的比较，无法裁决这一问题。",
     "支持一种部分延伸、叠加器官水平印记的混合解读，但有待更大多队列数据确认。", False),
    (r'在本报告中，[""]hub基因[""]按其操作性、选择规则的含义使用（在发现数据内由LASSO选择、随机森林排序与SHAP重要性一致提名），而非指蛋白互作网络中的度中心性。',
     "", True),
    ("这6个基因在生物学上相互关联，而非随机清单；这一推断来自已有文献，而非本研究产生的新实验证据。",
     "这6个基因在生物学上相互关联，而非随机清单。", False),
    ("这些基因的GO衍生分组（嘧啶代谢、RNA外切体、pantetheinase信号）应视为基于注释的类别，而非实验证明的共调控模块；现有数据并未确立这些类别间的协调调控，该分类是提出假说的框架，而非已确认的功能架构。",
     "这些基因的GO衍生分组（嘧啶代谢、RNA外切体、pantetheinase信号）是基于注释的类别；它们是否作为共调控模块发挥作用，需要实验检验。", False),
    ("而非转录上独立的实体；该假说建立在每组2至4个供体之上，需在更大队列中确认。",
     "而非转录上独立的实体。", False),
    ("倒U形应结合其跨队列来源解读（见上），其表面上与临床脓毒症高炎症向免疫抑制转变（9，10）的对齐与该转变相容，但不能作为该转变的证据；而",
     "倒U形反映了临床脓毒症高炎症向免疫抑制的转变（9，10），但受其跨队列来源所限，仅为与之相符的模式；而", False),
    ("ARDS模型向PBMC的可迁移性推测由其单核细胞成分承载：",
     "ARDS模型向PBMC的可迁移性由其单核细胞成分承载：", False),
    ("因此，迁移模型应解读为髓系区室标签，而非PBMC中的中性粒细胞读数。",
     "因此，在PBMC中，迁移模型是髓系区室标签，而非中性粒细胞读数。", False),
    ("在如此弱的遗传架构下，正向零结果不能排除中小因果效应（效能分析见结果），任何中等的真实效应都不会被检出。",
     "在如此弱的遗传架构下，任何中等的真实效应都不会被检出（效能分析见结果）。", False),
    ("但无法挽回低效能。因此，该零结果应解读为未检测到因果信号，而非无因果的证明。对大多数基因而言",
     "但无法挽回低效能。对大多数基因而言", False),
    ("该假说需要实验验证，而非优先靶点论断。",
     "该假说有待实验验证。", False),
    ("故CDK5RAP2的对接能量应作探索性解读（补充表S5）",
     "故CDK5RAP2的对接能量为探索性结果（补充表S5）", False),
    (r'论断。CMap糖皮质激素信号在期望上仍有部分循环性，尽管随机基因集零分布显示其前10集中度超过同背景随机集（真实标签前10出现糖皮质激素，同背景零基因集为17%，全转录组零基因集为0%；见[""]靶点优先级[""]）。脓毒症列线图',
     "论断。脓毒症列线图", True),
    ("也可能影响CMap糖皮质激素结果。最后，本研究的一切均为计算结果，有待前瞻性或实验验证。",
     "也可能影响CMap糖皮质激素结果。", False),
    ("计算筛选提名VNN1作为实验验证的候选，而非已确立的靶点。",
     "计算筛选提名VNN1作为该验证的首要候选。", False),
]

def apply_literal(path, edits):
    txt = open(path, encoding='utf-8').read()
    for i, (old, new) in enumerate(edits, 1):
        n = txt.count(old)
        assert n == 1, f'{path} edit #{i}: found {n} occurrences (expected 1): {old[:60]}...'
        txt = txt.replace(old, new)
    open(path, 'w', encoding='utf-8').write(txt)
    print(f'{path.split("/")[-1]}: {len(edits)} literal edits applied')

def apply_cn(path, edits):
    txt = open(path, encoding='utf-8').read()
    for i, (old, new, is_rx) in enumerate(edits, 1):
        if is_rx:
            hits = re.findall(old, txt)
            assert len(hits) == 1, f'CN edit #{i}: {len(hits)} regex hits: {old[:50]}'
            txt = re.sub(old, new, txt)
        else:
            n = txt.count(old)
            assert n == 1, f'CN edit #{i}: found {n} occurrences: {old[:50]}'
            txt = txt.replace(old, new)
    open(path, 'w', encoding='utf-8').write(txt)
    print(f'{path.split("/")[-1]}: {len(edits)} CN edits applied')

def disc_words(path, start, end):
    txt = open(path, encoding='utf-8').read()
    m = re.search(start + r'(.*?)' + end, txt, re.S)
    return len(m.group(1).split())

w0 = disc_words(EN, 'DISCUSSION\n', '\n---\n\nDATA AVAILABILITY')
apply_literal(EN, EN_EDITS)
w1 = disc_words(EN, 'DISCUSSION\n', '\n---\n\nDATA AVAILABILITY')
print(f'Discussion word count: {w0} -> {w1}')

apply_cn(CN, CN_EDITS)

entry = """

---

## 2026-10-07（二）：DISCUSSION 去防御性表达（14 处，用户确认后执行）

按已确认清单修改 Discussion：①重复堆叠免责各保留一处（混合解读 caveat、GO 分组非实验证据、donor 数、MR 正向零结果、CMap 循环性）；②6 处 "should be read as" 指导式框架保留 limitations 引导句 1 处、其余改直接陈述；③删预防性术语声明（"hub gene" 操作性含义整句）、弱断言（presumably）、防御性收尾（P4/P8 两遍 "rather than" 改正向）、泛化免责声明（"everything here is computational"）。所有具体 limitation（效能、LDSC、VNN2 选择性、MM/GBSA、n=24、非感染对照未测、nomogram 0.708、无临床协变量）全部保留。不新增事实、不改数据、不动引用编号。Discussion 字数变化见脚本输出。同步更新 02_Manuscript_Revised.md、03_Manuscript_Translation_CN.md、manuscript_revised.docx（step95 重生成）。脚本：scripts/step101_discussion_deai.py。
"""
open(LOG, 'a', encoding='utf-8').write(entry)
print('log appended')

r = subprocess.run([sys.executable, '/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/step95_make_revised_docx.py'],
                   capture_output=True, text=True)
print(r.stdout[-300:] if r.stdout else '')
if r.returncode != 0:
    print('STEP95 FAILED:', r.stderr[-1000:]); sys.exit(1)
print('docx regenerated')
