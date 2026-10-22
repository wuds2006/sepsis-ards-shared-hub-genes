#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
step104_add_context_dependence_para.py
2026-10-15: 在 Discussion 第一段与第二段之间新增"共享基因名单分析情境依赖"段，
内容为 CX3CR1/PID1/PTGDS（Jiang et al., Front Immunol 2024;15:1480542，本文参考文献 6）
是否进入本文 179 候选集的核对结论。三稿同步：02 EN md、03 CN md、04 修改日志；
docx 由 step95_make_revised_docx.py 重生成（本脚本调用）。
"""
import subprocess, sys, os

O9 = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission/09_review_deAI_2026-09-25"
EN = os.path.join(O9, "02_Manuscript_Revised.md")
CN = os.path.join(O9, "03_Manuscript_Translation_CN.md")
LOG = os.path.join(O9, "04_Modification_Log.md")

EN_PARA = (
    "The identity of shared genes depends on the analytic context. "
    "A recent screen of peripheral-blood cohorts with substantially more "
    "permissive differential expression criteria nominated *CX3CR1*, *PID1*, "
    "and *PTGDS* as shared sepsis-ARDS markers (6). In the present curated "
    "discovery data, all three genes are downregulated in sepsis blood, but "
    "none would have entered the intersection from the ARDS side: in "
    "GSE76293, *PID1* falls short of the pre-set fold-change threshold "
    "(log2 fold change \u22120.79, adjusted P=0.006), and *CX3CR1* and *PTGDS* "
    "show no significant differential expression (log2 fold change \u22120.25 "
    "and \u22120.27; adjusted P=0.691 and P=0.309). The divergence is not a "
    "refutation of that screen; it quantifies how strongly shared-gene lists "
    "depend on cohort composition and threshold choice, and it is the reason "
    "the present analysis fixed thresholds in advance, removed off-target "
    "samples before testing, and required candidates to survive independent "
    "validation."
)

CN_PARA = (
    "共享基因的名单取决于分析情境。近期一项基于外周血队列、差异表达标准宽松得多"
    "的筛查提名*CX3CR1*、*PID1*和*PTGDS*为脓毒症与ARDS的共享标志物（6）。"
    "在本研究整编后的发现数据中，这三个基因在脓毒症血液中均下调，但没有任何一个"
    "能从ARDS一侧进入交集：在GSE76293中，*PID1*的效应量未达到预设倍数变化阈值"
    "（log2倍数变化\u22120.79，校正P=0.006），*CX3CR1*和*PTGDS*则无显著差异表达"
    "（log2倍数变化分别为\u22120.25和\u22120.27；校正P分别为0.691和0.309）。"
    "这一分歧并非对该筛查的否定；它量化了共享基因名单对队列构成与阈值选择的依赖"
    "程度，也正是本研究预先固定阈值、在检验前剔除脱靶样本、并要求候选基因通过"
    "独立验证的原因。"
)

LOG_SECTION = """

---

## 2026-10-15：Discussion 补入"共享基因名单分析情境依赖"段（CX3CR1/PID1/PTGDS 核对）

背景：用户要求核对 Jiang et al.（Front Immunol 2024;15:1480542，本文参考文献 6）提名的 *CX3CR1*、*PID1*、*PTGDS* 是否在本文 179 个共享候选基因集中，并将核对结果补入 Discussion。

核对数据（01_DEG/GSE69528_DEG_all.csv、GSE76293_DEG_all.csv；候选集 03_intersection/Shared_CGs.txt）：

| 基因 | GSE69528（脓毒症） | GSE76293（ARDS） | 是否入 179 候选集 |
|---|---|---|---|
| *CX3CR1* | log2FC=−1.13，校正P=0.000195（Down DEG） | log2FC=−0.25，校正P=0.691（NS） | 否 |
| *PID1* | log2FC=−1.10，校正P=2.97×10⁻⁸（Down DEG） | log2FC=−0.79，校正P=0.006（显著但未达 \|log2FC\|>1，NS） | 否 |
| *PTGDS* | log2FC=−1.26，校正P=1.48×10⁻⁵（Down DEG） | log2FC=−0.27，校正P=0.309（NS） | 否 |

三基因均死于"两病交集"的 ARDS 侧 DEG 门（阈值 |log2FC|>1 且校正P<0.05）；其中 *CX3CR1*、*PTGDS* 在 GSE76293 连名义显著都不满足。FIMMU 文（文献 6）的检出依赖其不同的 ARDS 数据集（GSE32707、GSE10474）与宽松阈值（|fold change|>0.3、P<0.05）。

修改：在 Discussion 第一段（六基因总述）与第二段（生物学关联）之间新增一段，以文献 6 为外部对照，说明共享基因名单对队列构成与阈值选择的依赖，并衔接本文 curation-first 设计（预先固定阈值、检验前剔除脱靶样本、独立验证）的合理性。引用编号不变（复用文献 6，无新增、无重排）。

同步更新：02_Manuscript_Revised.md、03_Manuscript_Translation_CN.md、manuscript_revised.docx（step95 重生成）。脚本：scripts/step104_add_context_dependence_para.py。
"""


def insert_para(path, anchor, para):
    with open(path, encoding="utf-8") as f:
        txt = f.read()
    assert txt.count(anchor) == 1, f"anchor not unique in {path}: {txt.count(anchor)}"
    assert para not in txt, f"paragraph already present in {path}"
    txt = txt.replace(anchor, para + "\n\n" + anchor)
    with open(path, "w", encoding="utf-8") as f:
        f.write(txt)
    print(f"inserted -> {os.path.basename(path)}")


def append_log(path, section):
    with open(path, encoding="utf-8") as f:
        txt = f.read()
    assert "2026-10-15：Discussion 补入" not in txt, "log section already present"
    txt = txt.rstrip("\n") + section
    with open(path, "w", encoding="utf-8") as f:
        f.write(txt)
    print("appended -> 04_Modification_Log.md")


def main():
    insert_para(EN, "The six genes are biologically connected rather than a random list.", EN_PARA)
    insert_para(CN, "这6个基因在生物学上相互关联，而非随机清单。", CN_PARA)
    append_log(LOG, LOG_SECTION)

    # regenerate docx via step95
    r = subprocess.run(
        [sys.executable,
         "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/scripts/step95_make_revised_docx.py"],
        capture_output=True, text=True)
    print("step95 stdout:", r.stdout[-400:] if r.stdout else "")
    print("step95 stderr:", r.stderr[-400:] if r.stderr else "")
    if r.returncode != 0:
        sys.exit("docx regeneration FAILED")
    print("docx regenerated OK")


if __name__ == "__main__":
    main()
