#!/usr/bin/env python3
# Convert Table3_docking.md to .docx with a real table (italic gene targets, bold best value)
from docx import Document
from docx.shared import Pt

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission"

doc = Document()
style = doc.styles['Normal']
style.font.name = 'Times New Roman'
style.font.size = Pt(11)

doc.add_heading('Table 3. Molecular-docking binding energies (AutoDock Vina, kcal/mol) of ten '
                'natural-product ligands against the six hub-gene protein targets', level=1)

headers = ['Target', 'EGCG', 'Quercetin', 'Luteolin', 'Kaempferol', 'Baicalein',
           'Resveratrol', 'Curcumin', 'Berberine', 'Emodin', 'Tetrahydroalstonine']
# (gene, footnote_marker, {col_index: (value, bold)})
rows = [
    ('CDK5RAP2', '*', {}),
    ('CD177', '*', {0: ('-7.64', True)}),
    ('GRB10', '*', {0: ('-9.61', True)}),
    ('UPP1', '\u2020', {6: ('-8.85', True)}),
    ('EXOSC4', '*', {0: ('-7.59', True)}),
    ('VNN1', '\u2020', {0: ('-7.53', True), 2: ('-10.01', True)}),
]
data = [
    ['-6.09', '-6.37', '-6.42', '-6.18', '-6.63', '-6.21', '-6.34', '-6.68', '-6.91', '-6.62'],
    ['-7.64', '-7.54', '-7.24', '-6.95', '-6.98', '-5.82', '-5.76', '-6.68', '-7.54', '-7.62'],
    ['-9.61', '-8.41', '-8.30', '-8.18', '-8.03', '-7.13', '-8.29', '-8.41', '-8.41', '-9.52'],
    ['-8.29', '-8.00', '-8.12', '-7.41', '-7.88', '-6.61', '-8.85', '-7.10', '-7.57', '-8.28'],
    ['-7.59', '-7.04', '-7.37', '-6.87', '-6.62', '-6.59', '-6.27', '-7.11', '-6.68', '-7.39'],
    ['-7.53', '-9.84', '-10.01', '-9.57', '-9.83', '-9.40', '-9.18', '-6.62', '-7.70', '-8.67'],
]

table = doc.add_table(rows=1, cols=len(headers))
table.style = 'Light Grid Accent 1'
hdr = table.rows[0].cells
for i, h in enumerate(headers):
    hdr[i].text = h
    for p in hdr[i].paragraphs:
        for r in p.runs:
            r.bold = True

for (gene, marker, highlights), vals in zip(rows, data):
    cells = table.add_row().cells
    # gene name italic + footnote marker upright
    run = cells[0].paragraphs[0].add_run(gene)
    run.italic = True
    if marker:
        cells[0].paragraphs[0].add_run(marker)
    for i, v in enumerate(vals):
        cell = cells[i + 1]
        r = cell.paragraphs[0].add_run(v.replace('-', '\u2212'))
        if i in highlights and highlights[i][0] == v:
            r.bold = True

notes = [
    ('Values are AutoDock Vina binding energies (kcal/mol); more negative values indicate '
     'stronger predicted binding. The strongest predicted binding across all 60 protein\u2013ligand '
     'combinations was VNN1\u2013luteolin (\u221210.01 kcal/mol, shown in bold).'),
    ('* Docked at the AlphaFold centroid box (no co-crystallized ligand available).'),
    ('\u2020 Docked at the co-crystallized ligand site (UPP1, PDB 3EUF; VNN1, PDB 7SLV).'),
    ('Abbreviations: EGCG = epigallocatechin gallate.'),
    ('Note: The best ligand per target was EGCG for GRB10 (\u22129.61 kcal/mol), emodin for '
     'CDK5RAP2 (\u22126.91 kcal/mol), EGCG for CD177 (\u22127.64 kcal/mol), curcumin for UPP1 '
     '(\u22128.85 kcal/mol), EGCG for EXOSC4 (\u22127.59 kcal/mol), and luteolin for VNN1 '
     '(\u221210.01 kcal/mol). VNN1 showed the strongest predicted binding of all targets for the '
     'flavonoid ligands tested; these are computational predictions requiring experimental '
     'confirmation. The docking protocol was validated by redocking the co-crystallized ligands '
     '(BAU at 3EUF and 9TC at 7SLV, chain A): top-ranked poses reproduced the crystallographic '
     'poses with heavy-atom root-mean-square deviation of 1.54 and 0.69 \u00c5, respectively, '
     'below the 2.0 \u00c5 criterion. As an orthogonal affinity estimate, single-point MM/GBSA '
     'free energies were computed for the five targets with high-confidence structures '
     '(Supplementary Table S21); CDK5RAP2 was excluded from rescoring owing to moderate '
     'box-region model confidence (mean pLDDT 73.5; Supplementary Table S19), and its docking '
     'energies should be read as exploratory.'),
]
doc.add_paragraph('')
for n in notes:
    doc.add_paragraph(n)

doc.save(f"{BASE}/02_tables/Table3_docking.docx")
print("Saved Table3_docking.docx")
