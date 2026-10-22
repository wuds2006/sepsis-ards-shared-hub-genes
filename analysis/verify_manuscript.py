#!/usr/bin/env python3
"""Manuscript verification for Path-A restructured 02_Manuscript_Revised.md.
Checks: section word counts, body total <=3500, abstract <=250, title <=100 chars,
keywords <=5, citation completeness 1..N with range expansion (e.g. "(6-8)"), first-occurrence
ascending order, figure/table references, em dash, "P= " spacing, banned AI words, unique T5
sentence, paired-claim co-occurrence.
Usage: python verify_manuscript.py <manuscript.md> <n_refs>
"""
import re
import sys

path = sys.argv[1] if len(sys.argv) > 1 else "../AJRCCM_submission/09_review_deAI_2026-09-25/02_Manuscript_Revised.md"
n_refs = int(sys.argv[2]) if len(sys.argv) > 2 else 22

text = open(path, encoding="utf-8").read()

# ---------- section splitting (plain-text ALL-CAPS headings) ----------
SEC_RE = re.compile(r"^(ABSTRACT|INTRODUCTION|METHODS|RESULTS|DISCUSSION|CONCLUSION|REFERENCES|DATA AVAILABILITY"
                    r"|ETHICS APPROVAL AND CONSENT TO PARTICIPATE)\s*$", re.M)
marks = [(m.start(), m.group(1)) for m in SEC_RE.finditer(text)]
sections = {}
for i, (pos, name) in enumerate(marks):
    end = marks[i + 1][0] if i + 1 < len(marks) else len(text)
    sections[name] = text[pos:end]

def wc(s):
    return len(s.split())

def wcount_between(a, b):
    for name in [a, b]:
        if name not in sections:
            return None
    return wc(sections[a][sections[a].find("\n"):]) if a == b else None

# word counts: heading content only (exclude the heading line itself)
sec_words = {}
for name, body in sections.items():
    sec_words[name] = wc(body) - 1  # minus heading token

def sec_span(a, b):
    """word count from start of section a to start of section b (a..b inclusive)"""
    ia = [i for i, (p, n) in enumerate(marks) if n == a][0]
    ib = [i for i, (p, n) in enumerate(marks) if n == b][0]
    start = marks[ia][0]
    end = marks[ib][0] if ib > ia else len(text)
    return wc(text[start:end])

body_total = None
names = [n for _, n in marks]
if all(x in names for x in ["INTRODUCTION", "METHODS", "RESULTS", "DISCUSSION"]):
    body_total = (
        sec_words.get("INTRODUCTION", 0)
        + sec_words.get("METHODS", 0)
        + sec_words.get("RESULTS", 0)
        + sec_words.get("DISCUSSION", 0)
        + sec_words.get("CONCLUSION", 0)
    )

print("== section word counts ==")
for name in ["ABSTRACT", "INTRODUCTION", "METHODS", "RESULTS", "DISCUSSION", "CONCLUSION"]:
    print(f"{name}: {sec_words.get(name)}")
print(f"[not counted in body] ETHICS: {sec_words.get('ETHICS APPROVAL AND CONSENT TO PARTICIPATE')}")
print(f"[not counted in body] DATA AVAILABILITY: {sec_words.get('DATA AVAILABILITY')}")
print(f"BODY TOTAL (Intro..Conclusion): {body_total}  (target <= 3500)")
print(f"ABSTRACT: {sec_words.get('ABSTRACT')}  (target <= 250)")

# ---------- title ----------
lines = text.splitlines()
title = lines[0].strip()
print(f"\nTITLE: {title}")
print(f"TITLE chars: {len(title)}  (target <= 100)")

# ---------- keywords ----------
kw = re.search(r"Keywords:\s*(.+)", text)
if kw:
    kws = [k.strip() for k in re.split(r"[;,]", kw.group(1)) if k.strip()]
    print(f"KEYWORDS: {len(kws)}  -> {kws}  (target <= 5)")
else:
    print("KEYWORDS: NOT FOUND")

# ---------- citations with range expansion ----------
main = None
if "REFERENCES" in sections:
    main = text[: text.find("\nREFERENCES")]
else:
    main = text

cites = []
# remove supplementary/table parentheticals? No: numeric citations only.
for m in re.finditer(r"\(([^()]*\d[^()]*)\)", main):
    inner = m.group(1)
    # citation-like content: numbers, commas, hyphens/en-dashes, spaces, P/CI contexts excluded
    for part in inner.split(";"):
        part = part.strip()
        nums = re.findall(r"\d+", part)
        if not nums:
            continue
        # a part qualifies as citation if it is composed only of digits, separators, spaces,
        # possibly prefixed by nothing; exclude things like "P=0.249" "n=24" "AUC 0.753"
        token = re.sub(r"[\d,\s\u2013\u2014-]", "", part)
        if token:  # has non-citation characters
            continue
        # expand ranges: "6-8" / "6–8" / "9—12" style two-number ranges first
        rng = re.match(r"^(\d+)\s*[\u2013\u2014-]\s*(\d+)$", part)
        if rng:
            lo, hi = int(rng.group(1)), int(rng.group(2))
            if hi - lo > 30:
                continue
            cites.extend(range(lo, hi + 1))
            continue
        segs = part.split(",")
        expanded = []
        ok = True
        for seg in segs:
            seg = seg.strip()
            sub = re.match(r"^(\d+)\s*[\u2013\u2014-]\s*(\d+)$", seg)
            if sub:
                lo, hi = int(sub.group(1)), int(sub.group(2))
                if hi - lo > 30:
                    ok = False
                    break
                expanded.extend(range(lo, hi + 1))
            else:
                expanded.extend(int(x) for x in re.findall(r"\d+", seg))
        if ok:
            cites.extend(expanded)

cites = [c for c in cites if 1 <= c <= n_refs]
missing = [c for c in range(1, n_refs + 1) if c not in cites]
# first occurrence order
seen, order = set(), []
for c in cites:
    if c not in seen:
        seen.add(c)
        order.append(c)
ascending = all(order[i] < order[i + 1] for i in range(len(order) - 1))
print(f"\n== citations (n={n_refs}) ==")
print(f"cited: {sorted(set(cites))}")
print(f"missing: {missing}")
print(f"first-occurrence order: {order}")
print(f"strictly ascending at first mention: {ascending}")

# ---------- figures / tables ----------
figs = sorted(set(int(x) for x in re.findall(r"Figure\s+(\d+)", main)))
tabs = sorted(set(int(x) for x in re.findall(r"Table\s+(\d+)", main)))
print(f"\nmain Figures cited: {figs}")
print(f"main Tables cited: {tabs}")

# ---------- style checks ----------
emdash = text.count("\u2014")
print(f"\n== style ==")
print(f"em dash count: {emdash}  (target 0)")
p_space = re.findall(r"P=\s+\d|P\s+=\s*\d", text)
print(f"'P=' followed by space: {len(p_space)}  (target 0) {p_space[:5]}")
banned = {w: len(re.findall(rf"\b{w}\b", text, re.I)) for w in ["robust", "coherent", "key barrier", "delve", "pivotal", "seamless"]}
print(f"banned AI words: {banned}")
t5 = "does constrain causal effects of ARDS on sepsis risk"
print(f"T5 sentence occurrences: {text.count(t5)}  (target 1)")

# ---------- paired-claim checks ----------
body = sections.get("RESULTS", "") + sections.get("DISCUSSION", "") + sections.get("ABSTRACT", "") + sections.get("METHODS", "") + sections.get("CONCLUSION", "")
pairs = {
    "E26+E29 (KM P=0.002 & LRT P=0.16)": (r"P=0\.002", r"P=0\.16"),
    "E32+E34 (MDOR 2.26 & 1.04)": (r"2\.26", r"1\.04"),
    "UPP1 HR 0.64 & univariate P=0.640": (r"hazard ratio per SD 0\.64", r"P=0\.136 to 0\.640"),
    "no-causality & phenotypic misclassification": (
        r"cannot distinguish true absence of causality from phenotypic misclassification",
        r"markers of disease state rather than established causal targets",
    ),
    "AUC 0.753-0.998 & disease-vs-disease 0.50": (r"0\.753 to 0\.998", r"combined six-gene signature AUC=0\.50"),
}
print("\n== paired claims ==")
for name, (a, b) in pairs.items():
    print(f"{name}: {bool(re.search(a, text)) and bool(re.search(b, text))}")
