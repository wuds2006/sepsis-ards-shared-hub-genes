#!/usr/bin/env python3
# step94: De-AI revision of manuscript.md -> 02_Manuscript_Revised.md
# - reference renumbering by first appearance
# - supplementary Figure/Table/Methods renumbering by first appearance
# - statistical format normalization (P=, AUC=, n=, etc.)
# - targeted content edits from review report (must-fix items)
import re

SRC = '/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission/01_manuscript/manuscript.md'
DST = '/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission/09_review_deAI_2026-09-25/02_Manuscript_Revised.md'

text = open(SRC).read()
body, refs = text.split('\nREFERENCES\n', 1)

REF_MAP = {1:1,2:2,3:3,4:4,5:5,6:9,7:10,8:11,9:12,10:13,11:14,12:15,13:16,
           14:6,15:7,16:8,17:20,18:21,19:22,20:23,21:24,22:25,23:26,24:27,
           25:28,26:29,27:30,28:31,29:32,30:33,31:34,32:35,33:36,34:17,35:18,
           36:19,37:37,38:38}
SF_MAP = {1:1,2:2,18:3,19:4,3:5,4:6,5:7,6:8,7:9,8:10,9:11,10:12,11:13,12:14,
          13:15,14:16,15:17,16:18,17:19}
ST_MAP = {1:1,2:2,3:3,18:4,19:5,4:6,5:7,6:8,7:9,8:10,9:11,10:12,11:13,12:14,
          13:15,14:16,15:17,16:18,17:19,20:21,21:20}
SM_MAP = {1:1,4:2,5:6,6:7,7:3,8:4,9:5}

# ---- 1. reference citation renumbering (numeric-only parentheticals) ----
def renum_cite(m):
    inner = m.group(1)
    parts = []
    for tok in re.split(r'(,\s*)', inner):
        tok_s = tok.strip()
        if re.fullmatch(r'\d+\s*[–-]\s*\d+', tok_s):
            a, b = re.split(r'[–-]', tok_s)
            a, b = REF_MAP[int(a)], REF_MAP[int(b)]
            parts.append(f'{a}\u2013{b}')
        elif tok_s.isdigit():
            parts.append(str(REF_MAP[int(tok_s)]))
        else:
            parts.append(tok)
    return '(' + ''.join(parts) + ')'

body = re.sub(r'\(([\d,–\-\s]+)\)',
              lambda m: renum_cite(m) if re.fullmatch(r'[\d,–\-\s]+', m.group(1)) else m.group(0),
              body)

# ---- 2. supplementary series renumbering ----
def renum_supp(m):
    label, plural, grp = m.group(1), m.group(2), m.group(3)
    mp = {'Figure': SF_MAP, 'Table': ST_MAP, 'Methods': SM_MAP}[label]
    grp2 = re.sub(r'S(\d+)', lambda x: 'S' + str(mp[int(x.group(1))]), grp)
    return f'Supplementary {label}{plural} {grp2}'

body = re.sub(r'Supplementary\s+(Figure|Table|Methods)(s?)\s+((?:S\d+[,\s]*(?:and\s*)?)+)',
              renum_supp, body)

# ---- 3. statistical format normalization ----
for a, b in [
    ('P = ', 'P='), ('P < ', 'P<'), ('P > ', 'P>'), ('P ≤ ', 'P≤'),
    ('AUC = ', 'AUC='), ('AUC > ', 'AUC>'), ('AUC < ', 'AUC<'),
    ('n = ', 'n='), ('β = ', 'β='), ('rg = ', 'rg='), ('rho = ', 'rho='),
    ('NES = ', 'NES='), ('r = ', 'r='), ('Zsummary = ', 'Zsummary='),
    ('P_HEIDI > ', 'P_HEIDI>'), ('PP.H4 < ', 'PP.H4<'),
    ('|r| > ', '|r|>'), ('F > ', 'F>'), ('h2 ≤ ', 'h2≤'), ('α = ', 'α='),
    ('γ = ', 'γ='), (' × 10', '×10'), ('95% CI, ', '95% CI '),
    ('I² ', 'I²='), ('Δ = ', 'Δ='),
]:
    body = body.replace(a, b)

# ---- 4. targeted content edits (review must-fix / should-fix) ----
E = [
# Abstract: MR asymmetry
("Mendelian randomization was null but uninformative given near-zero single-nucleotide polymorphism heritability.",
 "Mendelian randomization detected no causal signal; the forward analysis was underpowered at plausible effect sizes (minimum detectable odds ratio 2.26), whereas the reverse analysis was informative (minimum detectable odds ratio 1.04)."),
# Abstract: bound complement claim
("complement activation was the only shared pathway.",
 "complement activation was the only immune pathway shared among the top-ranked dysregulated processes of both diseases."),
# Intro: soften blanket claim on prior screens
("but they were built on processed matrices without off-target sample removal, nominated and validated candidates within the same data, and did not test genetic causality or account for cell-level pseudoreplication.",
 "but they were built on processed matrices, to our knowledge without systematic off-target sample removal, nominated and validated candidates within the same data, and did not test genetic causality or account for cell-level pseudoreplication."),
# Intro: compress closing results preview
("The working premise was simple: if sepsis and ARDS share a core of coordinately regulated hub genes, those genes should surface consistently across cohorts, cell types, and analytic angles; they did. The six genes were nominated entirely within discovery data and confirmed as single-gene markers in independent validation cohorts, and the shared immune axis centered on complement activation, with myeloid expansion and cytotoxic-lymphocyte depletion in blood. Two in silico findings carry hypothesis-level therapeutic interest: glucocorticoid repurposing predictions, and VNN1 as the most structurally supported of the docking-derived target hypotheses for flavonoid natural products (predicted binding energy as low as −10.01 kcal/mol), an inference that requires experimental validation rather than a prioritized drug-target claim. The same work demonstrates a reproducible, curation-first pipeline for cross-disease transcriptomic discovery.",
 "The working premise was simple: if sepsis and ARDS share a core of coordinately regulated hub genes, those genes should surface consistently across cohorts, cell types, and analytic angles. We tested that premise under a curation-first, discovery-strict design, and we bounded every downstream claim, from immune deconvolution to drug prediction, by orthogonal or null-based controls."),
# Methods: add GSE100159 exclusion
("External validation used GSE100159 (sepsis; 33 cases, 12 controls) and GSE243066",
 "External validation used GSE100159 (sepsis; 33 cases and 12 controls, after exclusion of two recovery-phase samples) and GSE243066"),
# Methods heading
("Immune infiltrations\n", "Immune Infiltration\n"),
# Methods ML: remove outcome leak
("External validation cohorts were used to confirm the pre-specified cross-disease candidates and, in the same step, to screen the seven ARDS-nominated candidates down to the final set: they confirmed six candidates in both cohorts and refuted *MUC1* (AUC=0.575 in the ARDS validation cohort), after which the six confirmed genes (*CDK5RAP2*, *CD177*, *GRB10*, *UPP1*, *EXOSC4*, *VNN1*) were defined as the core hub genes. Because this screening step consumed validation information, performance estimates for the six genes in these cohorts are mildly optimistic; a sensitivity analysis refitting the prediction models with all seven nominated genes (including *MUC1*) left the ARDS validation AUC unchanged (0.984) and slightly reduced the sepsis transfer AUC (0.708 to 0.699), indicating that the screening step did not inflate the reported validation performance.",
 "External validation cohorts were then used to confirm the pre-specified cross-disease candidates and, in the same step, to screen the nominated candidates down to the final set (outcomes in Results). Because this screening step consumed validation information, performance estimates in these cohorts are mildly optimistic and were bounded by a seven-gene sensitivity analysis (see Results)."),
# Results: fix exact lower bound
("exact one-sided 95% lower confidence bounds of 0.994 (GSE243066) and 0.992 (GSE100159)",
 "exact one-sided 95% lower confidence bounds of 0.994 (GSE243066) and 0.993 (GSE100159)"),
# Discussion: beyond-sepsis-state rephrase
("in which ARDS monocytes sit beyond the sepsis state along a shared differentiation axis while MIF signaling dominates the ARDS interactome",
 "in which disease-group monocytes occupy later pseudotime positions than their respective controls along a shared differentiation axis while MIF signaling dominates the ARDS interactome"),
# Discussion: compress MR power repetition
("under such weak genetic architecture, a null MR cannot exclude small causal effects: the forward analysis had 28% power at a causal odds ratio of 1.5 and could detect, at 80% power, only effects of odds ratio 2.26 or larger, so any modest true effect would have gone undetected, and our orthogonal checks",
 "under such weak genetic architecture the forward null cannot exclude small or moderate causal effects (power analysis in Results), so any modest true effect would have gone undetected. Our orthogonal checks"),
# Discussion: compress CMap enumeration
("Our random gene-set nulls quantify how much information remains after this expectation is accounted for (see Results): in the up-regulated query the six-hub-gene signature returned dexamethasone as the rank-one term (P=3.7×10⁻⁴), an enrichment level attained or exceeded by 13 of 100 size-matched random signatures drawn from the non-hub candidate genes, with a glucocorticoid in the top ten for 17 of those nulls versus the real signature's every ranked position, and no transcriptome-wide random set produced a top-ten glucocorticoid or any glucocorticoid hit at all; 87 of the 100 candidate-background nulls contained at least one glucocorticoid term somewhere in their full results, confirming that glucocorticoid matching is the expected behavior of inflammatory signatures. The candidate-gene background shares the disease co-expression structure of the real signature and therefore quantifies how much of the signal any similarly structured inflammatory signature would reproduce (13 of 100), whereas the transcriptome-wide background matches neither disease relevance nor coordinated up-regulation and therefore overstates specificity; a stringency-matched background of random six-gene sets drawn from the genes up-regulated in both discovery cohorts matched or exceeded the real glucocorticoid P value in only 2 of 100 and placed a glucocorticoid in the top ten in 18 of 100, indicating that the enrichment level is uncommon even among disease-relevant, coordinated inflammatory signatures of the same size.",
 "Our random gene-set nulls quantify how much information remains after this expectation is accounted for (full distributions in Results). Against the candidate-gene background, which shares the disease co-expression structure of the real signature, 13 of 100 random signatures matched or exceeded the observed enrichment; against the transcriptome-wide background, none did, and most candidate-background nulls (87 of 100) returned at least one glucocorticoid term somewhere in their full results. These two backgrounds bracket rather than match the real signature's provenance, so a third, stringency-matched background of random six-gene sets drawn from genes up-regulated in both discovery cohorts was tested; it matched or exceeded the real enrichment in only 2 of 100 cases, indicating that the enrichment level is uncommon even among disease-relevant, coordinated inflammatory signatures of the same size."),
# Discussion: un-bold the two run-in headings
("**Shared hub genes: reliable markers of disease state, but causal roles remain unverified**",
 "Shared hub genes: reliable markers of disease state, but causal roles remain unverified"),
("**Target prioritization: druggability is concentrated in VNN1 and UPP1**",
 "Target prioritization: druggability is concentrated in VNN1 and UPP1"),
]
for a, b in E:
    if a not in body:
        print('MISS:', a[:70])
    body = body.replace(a, b)

# ---- 5. reference list rebuild ----
entries = {}
for m in re.finditer(r'(?m)^(\d+)\.\s+(.*?)(?=\n\d+\.\s|\Z)', refs, re.S):
    entries[int(m.group(1))] = ' '.join(m.group(2).split())
new_refs = {}
for old, new in REF_MAP.items():
    new_refs[new] = entries[old]
out_refs = '\n'.join(f'{n}. {new_refs[n]}' for n in sorted(new_refs))

open(DST, 'w').write(body + '\nREFERENCES\n' + out_refs + '\n')
print('written', DST, len(body + out_refs), 'chars')
