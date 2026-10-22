#!/usr/bin/env python3
"""step96: Rewire main-figure citations from the old Figures_final set to the
03_figures_RESULTS_NPG set (Plan B), in both the EN revised manuscript and the
CN translation. Every replacement asserts its expected occurrence count.
"""
import sys

BASE = '/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission/09_review_deAI_2026-09-25'
EN = f'{BASE}/02_Manuscript_Revised.md'
CN = f'{BASE}/03_Manuscript_Translation_CN.md'

EN_REPL = [
    # Results 1: shared candidate genes
    ('(Figure 2; Supplementary Table S7)', '(Figure 1; Supplementary Table S7)', 1),
    # Results 2: enrichment
    ('Figure 3B', 'Figure 2A', 1),
    ('Figure 3C', 'Figure 2B', 2),
    ('Figure 3D', 'Figure 2C', 2),
    # Results 3: ML + PPI
    ('(Figure 4F; Supplementary Figure S5)', '(Figure 3A, B; Supplementary Figure S5)', 1),
    ('marker (Figure 3A)', 'marker (Figure 3D)', 1),
    # Results 4: boxplots unchanged; GSEA gains main-panel citation
    ('differed by tissue context (Supplementary Figure S6).',
     'differed by tissue context (Figure 4E; Supplementary Figure S6).', 1),
    # Results 5: Figure 5 unchanged
    # Results 6: external confirmation
    ('(Figure 4E, Table 2; case-to-control', '(Figure 6; Table 2; case-to-control', 1),
    ('differ by disease (Supplementary Table S11).',
     'differ by disease (Figure 6B; Supplementary Table S11).', 1),
    # Results 7: MR gains main-figure citations
    ('was likewise null (IVW P=0.329).', 'was likewise null (IVW P=0.329) (Figure 7A, B).', 1),
    ('(minimum detectable odds ratio 1.04) (Supplementary Table S13).',
     '(minimum detectable odds ratio 1.04) (Figure 7C; Supplementary Table S13).', 1),
    ('(PP.H4<0.05 for all) (Supplementary Figure S8, Supplementary Table S14).',
     '(PP.H4<0.05 for all) (Figure 7D; Supplementary Figure S8; Supplementary Table S14).', 1),
    # Results 8: immune microenvironment
    ('(Figure 6A).', '(Figure 8A).', 1),
    ('(Figures 6B–E and 7B, Supplementary Table S2)', '(Figure 8B–E; Supplementary Table S2)', 1),
    ('(Figure 7A, Supplementary Figure S9)', '(Figure 8F; Supplementary Figure S9)', 1),
    ('(Figure 8; Supplementary Table S15)', '(Figure 8G, H; Supplementary Table S15)', 1),
    # Results 9: single cell
    ('complementary windows on the disease spectrum: PBMCs',
     'complementary windows on the disease spectrum (Figure 9): PBMCs', 1),
    ('one early-sepsis donor excepted; Figure 9C)', 'one early-sepsis donor excepted; Figure 9C, G)', 1),
    # Results 10: protein validation and therapeutics
    ('VNN2 was upregulated in all 10 diseases (Supplementary Figure S19).',
     'VNN2 was upregulated in all 10 diseases (Figure 10A; Supplementary Figure S19).', 1),
    ('(Figure 10C; Supplementary Figure S3, Supplementary Table S4)',
     '(Figure 10B; Supplementary Figure S3; Supplementary Table S4)', 1),
    ('(Figure 10, Table 3, Supplementary Table S19)',
     '(Figure 10C, D; Table 3; Supplementary Table S19)', 1),
    ('(Figure 10D).', '(Figure 10E).', 1),
]

CN_REPL = [
    ('（图2；补充表S7）', '（图1；补充表S7）', 1),
    ('图3B', '图2A', 1),
    ('图3C', '图2B', 2),
    ('图3D', '图2C', 2),
    ('（图4F；补充图S5）', '（图3A、B；补充图S5）', 1),
    ('（图3A）', '（图3D）', 1),
    ('炎症富集（补充图S6）。', '炎症富集（图4E；补充图S6）。', 1),
    ('（图4E，表2；', '（图6；表2；', 1),
    ('标志物相符（补充表S11）。', '标志物相符（图6B；补充表S11）。', 1),
    ('（IVW P=0.329）。', '（IVW P=0.329）（图7A、B）。', 1),
    ('（补充表S13）。', '（图7C；补充表S13）。', 1),
    ('（均PP.H4<0.05）（补充图S8，补充表S14）。', '（均PP.H4<0.05）（图7D；补充图S8；补充表S14）。', 1),
    ('（图6A）。', '（图8A）。', 1),
    ('（图6B-E与7B，补充表S2）', '（图8B-E；补充表S2）', 1),
    ('（图7A，补充图S9）', '（图8F；补充图S9）', 1),
    ('（图8；补充表S15）', '（图8G、H；补充表S15）', 1),
    ('两个单细胞队列提供了疾病谱的互补窗口：', '两个单细胞队列提供了疾病谱的互补窗口（图9）：', 1),
    ('1例早期脓毒症供体例外；图9C）', '1例早期脓毒症供体例外；图9C、G）', 1),
    ('（补充图S19）。治疗学分析', '（图10A；补充图S19）。治疗学分析', 1),
    ('（图10C；补充图S3，补充表S4）', '（图10B；补充图S3；补充表S4）', 1),
    ('（图10，表3，补充表S19）', '（图10C、D；表3；补充表S19）', 1),
    ('（图10D）。', '（图10E）。', 1),
]


def apply(path, repls):
    t = open(path, encoding='utf-8').read()
    for old, new, n in repls:
        c = t.count(old)
        assert c == n, f'{path}: expected {n} of {old!r}, found {c}'
        t = t.replace(old, new)
    open(path, 'w', encoding='utf-8').write(t)
    print(f'{path}: {len(repls)} replacements OK')


apply(EN, EN_REPL)
apply(CN, CN_REPL)
print('DONE')
