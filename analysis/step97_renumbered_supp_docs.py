#!/usr/bin/env python3
"""step97: Build renumbered Supplementary_Materials.md and Supplementary_Methods.md
for the 09_review_deAI_2026-09-25 package.

Renumbering (identical to step94 maps used for the revised manuscript):
  SF_MAP / ST_MAP / SM_MAP below.
Methods sections cited in the main text are renumbered S1-S7 by first-citation
order; uncited sections keep their relative order as S8-S13.
Internal cross-references, bare Table/Figure S refs, and stats formatting are
normalized to match the revised manuscript.
"""
import re

SRC = '/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission/05_supplementary'
DST = '/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission/09_review_deAI_2026-09-25'

SF_MAP = {1:1, 2:2, 18:3, 19:4, 3:5, 4:6, 5:7, 6:8, 7:9, 8:10, 9:11, 10:12,
          11:13, 12:14, 13:15, 14:16, 15:17, 16:18, 17:19}
ST_MAP = {1:1, 2:2, 3:3, 18:4, 19:5, 4:6, 5:7, 6:8, 7:9, 8:10, 9:11, 10:12,
          11:13, 12:14, 13:15, 14:16, 15:17, 16:18, 17:19, 20:21, 21:20}
# methods: cited sections by first-citation order, then uncited in original order
SM_MAP = {1:1, 4:2, 7:3, 8:4, 9:5, 5:6, 6:7, 2:8, 3:9, 10:10, 11:11, 12:12, 13:13}

STATS = [
    ('P = ', 'P='), ('P < ', 'P<'), ('P > ', 'P>'), ('P ≤ ', 'P≤'),
    ('AUC = ', 'AUC='), ('AUC > ', 'AUC>'), ('AUC < ', 'AUC<'),
    ('n = ', 'n='), ('r = ', 'r='), ('rho = ', 'rho='), ('F = ', 'F='),
    ('rg = ', 'rg='), ('h² = ', 'h²='), ('R² = ', 'R²='), ('β = ', 'β='),
    ('delta = ', 'delta='), ('power = ', 'power='), ('Q = ', 'Q='),
]


def norm_stats(t):
    for a, b in STATS:
        t = t.replace(a, b)
    return t


def renum_crossrefs(t):
    """Renumber Supplementary Table/Figure and bare Table/Figure S refs in body text."""
    t = re.sub(r'Supplementary Table[s]? S(\d+)',
               lambda m: f'Supplementary Table S{ST_MAP[int(m.group(1))]}', t)
    t = re.sub(r'Supplementary Figure[s]? S(\d+)',
               lambda m: f'Supplementary Figure S{SF_MAP[int(m.group(1))]}', t)
    t = re.sub(r'(?<!Supplementary )(?<!Supplementary)\bTable S(\d+)',
               lambda m: f'Table S{ST_MAP[int(m.group(1))]}', t)
    t = re.sub(r'(?<!Supplementary )(?<!Supplementary)\bFigure S(\d+)',
               lambda m: f'Figure S{SF_MAP[int(m.group(1))]}', t)
    return t


def split_sections(text, kind):
    """Split into (preamble, [(old_num, heading_line, body)], tail).
    kind = 'Table' or 'Figure' matching '### Table S1.' style headings.
    Sections are delimited by headings; optional '---' separators attach to
    the following section's body start (kept with preceding body as-is)."""
    lines = text.split('\n')
    pat = re.compile(rf'^### {kind} S(\d+)\.\s*(.*)$')
    idx = [i for i, l in enumerate(lines) if pat.match(l)]
    if not idx:
        return text, [], ''
    preamble = '\n'.join(lines[:idx[0]])
    secs = []
    for k, i in enumerate(idx):
        j = idx[k + 1] if k + 1 < len(idx) else len(lines)
        m = pat.match(lines[i])
        secs.append((int(m.group(1)), lines[i], '\n'.join(lines[i + 1:j]).strip('\n')))
    return preamble, secs, ''


def rebuild(prefix, secs, mp, kind):
    out = [prefix.rstrip('\n')]
    for old, head, body in sorted(secs, key=lambda s: mp[s[0]]):
        new = mp[old]
        title = re.match(rf'^### {kind} S\d+\.\s*(.*)$', head).group(1)
        out.append(f'### {kind} S{new}. {title}')
        out.append(renum_crossrefs(norm_stats(body)))
    return '\n\n'.join(out) + '\n'


# ---------------- Materials ----------------
mat = open(f'{SRC}/Supplementary_Materials.md', encoding='utf-8').read()
head_end = mat.index('## Supplementary Tables')
fig_marker = mat.index('## Supplementary Figure Legends')
manifest_marker = mat.index('## File Manifest')

header = mat[:head_end]
tables_txt = mat[head_end:fig_marker]
figs_txt = mat[fig_marker:manifest_marker]

tab_pre, tab_secs, _ = split_sections(tables_txt, 'Table')
fig_pre, fig_secs, _ = split_sections(figs_txt, 'Figure')
assert sorted(n for n, _, _ in tab_secs) == list(range(1, 22)), [n for n, _, _ in tab_secs]
assert sorted(n for n, _, _ in fig_secs) == list(range(1, 20)), [n for n, _, _ in fig_secs]

new_header = (
    '# Supplementary Materials\n\n'
    'This file accompanies the manuscript "Integrated Transcriptomic and '
    'Machine-Learning Analysis Identifies Shared Hub Genes and Immune Signatures '
    'of Sepsis and Acute Respiratory Distress Syndrome." Tables and figures are '
    'numbered S1-S21 and S1-S19 by first-citation order in the revised main text '
    '(renumbering of 2026-09-25; see 05_Renumbering_Maps.md). Source paths refer '
    'to the analysis archive and retain their original directory names; the '
    'renumbered files accompanying this package are in figures/supplementary/ '
    'and figures/supplementary_tables/.\n\n'
    '## Contents\n\n'
    '- **Supplementary Methods** (see `Supplementary_Methods.md`)\n'
    '- **Supplementary Tables** (below)\n'
    '- **Supplementary Figure Legends** (below)\n\n---\n'
)

tables_out = rebuild('## Supplementary Tables', tab_secs, ST_MAP, 'Table')
figs_out = rebuild('## Supplementary Figure Legends', fig_secs, SF_MAP, 'Figure')

manifest = (
    '## File Manifest\n\n'
    '| File | Description |\n'
    '|------|-------------|\n'
    '| `Supplementary_Methods.md` | Detailed methods (S1-S13; S1-S7 follow first-citation order in the main text) |\n'
    '| `Supplementary_Materials.md` | This index (Tables S1-S21 and Figure legends S1-S19) |\n'
    '| `figures/supplementary/FigureS1-S19.*` | Supplementary figures (renumbered), PDF/PNG/TIFF triple export |\n'
    '| `figures/supplementary_tables/TableS1-S21.xlsx` | Supplementary tables (renumbered, multi-sheet workbooks) |\n'
    '| `figures/supplementary_tables/CellChat_LR_*.csv` | Complete differential ligand-receptor interaction lists |\n'
)

materials = new_header + '\n' + tables_out + '\n---\n\n' + figs_out + '\n' + manifest
open(f'{DST}/Supplementary_Materials.md', 'w', encoding='utf-8').write(materials)

# ---------------- Methods ----------------
meth = open(f'{SRC}/Supplementary_Methods.md', encoding='utf-8').read()
lines = meth.split('\n')
pat = re.compile(r'^## S(\d+)\.\s*(.*)$')
idx = [i for i, l in enumerate(lines) if pat.match(l)]
secs = []
for k, i in enumerate(idx):
    j = idx[k + 1] if k + 1 < len(idx) else len(lines)
    m = pat.match(lines[i])
    secs.append((int(m.group(1)), m.group(2), '\n'.join(lines[i + 1:j]).strip('\n')))
assert sorted(n for n, _, _ in secs) == list(range(1, 14)), [n for n, _, _ in secs]

out = ['# Supplementary Methods',
       '',
       'Sections S1-S7 are numbered by first-citation order in the revised main '
       'text; sections S8-S13 provide additional detail in their original order '
       '(renumbering of 2026-09-25; see 05_Renumbering_Maps.md).']
for old, title, body in sorted(secs, key=lambda s: SM_MAP[s[0]]):
    new = SM_MAP[old]
    body = renum_crossrefs(norm_stats(body))
    # method-internal section cross-reference: "(S4)" -> "(S2)" etc.
    body = re.sub(r'\(S(\d+)\)', lambda m: f'(S{SM_MAP[int(m.group(1))]})', body)
    out.append(f'## S{new}. {title}')
    out.append(body)
open(f'{DST}/Supplementary_Methods.md', 'w', encoding='utf-8').write(
    '\n\n'.join(out) + '\n')

print('Materials: tables', len(tab_secs), 'figures', len(fig_secs))
print('Methods: sections', len(secs))
print('DONE')
