#!/usr/bin/env python3
# Convert manuscript.md and Table1 to .docx with proper gene-name italics
import re
from docx import Document
from docx.shared import Pt, RGBColor, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission"

def add_inline(paragraph, text):
    """Parse **bold** and *italic* inline markdown."""
    tokens = re.split(r'(\*\*[^*]+\*\*|\*[^*]+\*)', text)
    for tok in tokens:
        if not tok:
            continue
        if tok.startswith('**') and tok.endswith('**') and len(tok) > 4:
            r = paragraph.add_run(tok[2:-2]); r.bold = True
        elif tok.startswith('*') and tok.endswith('*') and len(tok) > 2:
            r = paragraph.add_run(tok[1:-1]); r.italic = True
        else:
            paragraph.add_run(tok)

def md_to_docx(md_path, docx_path):
    doc = Document()
    style = doc.styles['Normal']
    style.font.name = 'Times New Roman'
    style.font.size = Pt(12)
    with open(md_path, encoding='utf-8') as f:
        lines = f.read().split('\n')
    for line in lines:
        s = line.rstrip()
        if not s.strip():
            continue
        if s.startswith('<!--'):
            continue  # skip HTML comments (author notes)
        if s.startswith('### '):
            h = doc.add_heading(level=3); add_inline(h, s[4:])
        elif s.startswith('## '):
            h = doc.add_heading(level=2); add_inline(h, s[3:])
        elif s.startswith('# '):
            h = doc.add_heading(level=1); add_inline(h, s[2:])
        elif s.startswith('- ') or s.startswith('* '):
            p = doc.add_paragraph(style='List Bullet'); add_inline(p, s[2:])
        elif re.match(r'^\d+\.\s', s):
            p = doc.add_paragraph(style='List Number'); add_inline(p, re.sub(r'^\d+\.\s', '', s))
        else:
            p = doc.add_paragraph(); add_inline(p, s)
    doc.save(docx_path)
    print("Saved", docx_path)

# Manuscript
md_to_docx(f"{BASE}/01_manuscript/manuscript.md", f"{BASE}/01_manuscript/manuscript.docx")

# Table 1 as a real docx table
def table1_to_docx():
    doc = Document()
    doc.add_heading('Table 1. Characteristics of the transcriptomic datasets used in this study', level=1)
    headers = ['Characteristic', 'GSE76293', 'GSE69528', 'GSE100159', 'GSE243066', 'GSE151263', 'GSE167363']
    rows = [
        ['Platform', 'Affymetrix HG-U133 Plus 2.0 (GPL570)', 'Illumina HumanHT-12 v4 (GPL10558)',
         'Illumina HumanWG-6 v3.0 (GPL6884)', 'Illumina NovaSeq 6000 RNA sequencing (GPL30209)',
         '10x Genomics', '10x Genomics'],
        ['Disease', 'ARDS', 'Sepsis', 'Sepsis', 'ARDS', 'ARDS + Sepsis', 'Sepsis'],
        ['Tissue / sample', 'Peripheral-blood polymorphonuclear neutrophils', 'Whole blood', 'Whole blood',
         'Peripheral-blood mononuclear cells', 'Peripheral-blood mononuclear cells',
         'Peripheral-blood CD14+ monocytes'],
        ['Cases, n', '12', '83', '33', '34', '3 ARDS, 4 sepsis',
         '2 sepsis (3 samples; one patient sampled at two time points)'],
        ['Controls, n', '12 (healthy)', '28 (uninfected)', '12 (healthy)', '15 (healthy)', 'N/A', '2 (healthy)'],
        ['Excluded samples, n', '70 (in vitro-cultured PMNs)', '27 (uninfected type 2 diabetes)', '2 (recovery)',
         'N/A', 'N/A', 'N/A'],
        ['Total samples, n', '94', '138', '47', '49', '7', '5 (12,359 cells)'],
        ['Role in study', 'Discovery', 'Discovery', 'External validation', 'External validation',
         'Single-cell analysis', 'Single-cell analysis'],
    ]
    table = doc.add_table(rows=1, cols=len(headers))
    table.style = 'Light Grid Accent 1'
    hdr = table.rows[0].cells
    for i, h in enumerate(headers):
        hdr[i].text = h
        for p in hdr[i].paragraphs:
            for r in p.runs:
                r.bold = True
    for row in rows:
        cells = table.add_row().cells
        for i, val in enumerate(row):
            cells[i].text = val
    doc.add_paragraph('')
    doc.add_paragraph('Abbreviations: ARDS = acute respiratory distress syndrome; PMN = polymorphonuclear neutrophil; T2DM = type 2 diabetes mellitus.')
    doc.save(f"{BASE}/02_tables/Table1_datasets.docx")
    print("Saved Table1_datasets.docx")

table1_to_docx()

# Table 2 as a real docx table with italic gene names
def table2_to_docx():
    doc = Document()
    doc.add_heading('Table 2. External confirmation of the six hub genes and combined signature', level=1)
    headers = ['Gene', 'GSE100159 (Sepsis) AUC', 'GSE100159 95% CI', 'GSE243066 (ARDS) AUC', 'GSE243066 95% CI']
    rows = [
        ('CD177', '0.990', '0.971–1.000', '0.998', '0.993–1.000'),
        ('UPP1', '0.753', '0.540–0.965', '0.935', '0.863–1.000'),
        ('GRB10', '0.982', '0.947–1.000', '0.971', '0.929–1.000'),
        ('VNN1', '0.995', '0.983–1.000', '0.982', '0.954–1.000'),
        ('EXOSC4', '0.995', '0.983–1.000', '0.831', '0.716–0.947'),
        ('CDK5RAP2', '0.934', '0.864–1.000', '0.782', '0.655–0.909'),
        ('Combined signature (in-cohort refit, descriptive)', '1.000', '\u22650.993 (exact)', '1.000', '\u22650.994 (exact)'),
    ]
    table = doc.add_table(rows=1, cols=len(headers))
    table.style = 'Light Grid Accent 1'
    hdr = table.rows[0].cells
    for i, h in enumerate(headers):
        hdr[i].text = h
        for p in hdr[i].paragraphs:
            for r in p.runs:
                r.bold = True
    for row in rows:
        cells = table.add_row().cells
        g = cells[0].paragraphs[0]
        if row[0].startswith('Combined'):
            run = g.add_run(row[0]); run.bold = True
        else:
            run = g.add_run(row[0]); run.italic = True
        for i in range(1, 5):
            cells[i].text = row[i]
    doc.add_paragraph('')
    doc.add_paragraph('Abbreviations: AUC = area under the receiver-operating-characteristic curve; CI = confidence interval.')
    doc.add_paragraph('Note: Confidence intervals for single genes were estimated by the DeLong method; single-gene AUCs involve no model fitting. The combined-signature row is an in-cohort refit of a logistic-regression signature within each validation cohort, reported for completeness and not used in gene selection; for perfect separation, an exact one-sided 95% lower confidence bound was computed from the rank statistic (alpha^(1/[n_case x n_control])) because the DeLong interval degenerates at AUC = 1.000. The primary multi-gene validation is the coefficient-transfer analysis (validation AUC 0.984 for ARDS and 0.708 for sepsis; see Disease Prediction Model). UPP1 showed lower discrimination in the sepsis cohort (AUC = 0.753) than in the ARDS cohort (AUC = 0.935). External validation played no role in gene selection: candidates were nominated by machine learning within discovery cohorts and screened cross-disease within discovery data (Supplementary Table S1); validation confirmed six of seven ARDS-nominated candidates and refuted MUC1 (AUC = 0.575 in GSE243066).')
    doc.save(f"{BASE}/02_tables/Table2_validation_AUC_CI.docx")
    print("Saved Table2_validation_AUC_CI.docx")

table2_to_docx()
