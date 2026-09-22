#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re
from docx import Document

SRC = r"E:/workbuddy/2026-07-19-16-45-02/corrected/manuscript_GeneticEpidemiology_v9_final.docx"
doc = Document(SRC)
paras = [p.text for p in doc.paragraphs]

# Locate key headings
for i, p in enumerate(paras):
    s = p.strip()
    if re.match(r"^(abstract|keywords|introduction|methods|results|discussion|references|figure legends)$", s, flags=re.IGNORECASE):
        print(f"[{i:3}] === {s.upper()} ===")

print("\n----- ABSTRACT AREA -----")
# find abstract
ai = next(i for i,p in enumerate(paras) if p.strip().lower()=="abstract")
words=[]
for p in paras[ai+1:ai+15]:
    if re.match(r"^(keywords|introduction|methods|results|discussion)", p.strip(), flags=re.I):
        break
    if p.strip():
        words.append(p)
print("WORDCOUNT:", len(" ".join(words).split()))
for w in words:
    print("  |", w[:120])

print("\n----- GTEx enrichment % occurrences -----")
for i,p in enumerate(paras):
    if re.search(r"\d+\.?\d*\s*%", p) and re.search(r"gtex|enrich|51|48|44", p, re.I):
        print(f"[{i}] {p[:160]}")

print("\n----- References vs Figure Legends boundary -----")
ri = next(i for i,p in enumerate(paras) if p.strip().lower()=="references")
fi = next((i for i,p in enumerate(paras) if p.strip().lower()=="figure legends"), None)
print("References at", ri, "Figure Legends at", fi)
refblock = paras[ri+1:fi] if fi else paras[ri+1:]
refblock = [r for r in refblock if r.strip()]
print("Actual reference entries:", len(refblock))
no_year=[r for r in refblock if not re.search(r"\(\d{4}\)", r)]
no_doi=[r for r in refblock if "doi.org" not in r.lower() and "doi:" not in r.lower() and "pmid" not in r.lower()]
print("Missing (Year):", len(no_year))
for r in no_year: print("   Y?", r[:90])
print("Missing DOI/PMID:", len(no_doi))
for r in no_doi: print("   D?", r[:90])
