# -*- coding: utf-8 -*-
"""FIG6 RECON: 复现 03_redraw_Fig6.py 的数据抽取，与 SI Table S18 逐格对账。"""
import numpy as np
from scipy import stats
from docx import Document
from docx.table import Table
from docx.oxml.ns import qn
import paths_config as P  # 统一路径入口（2026-09-20）

AF = P.need(P.AF1, 'Additional file 1（从期刊补充材料下载后用 AF1_DOCX 指定）')
d = Document(AF)
tabs = [Table(ch, d) for ch in d.element.body.iterchildren() if ch.tag == qn('w:tbl')]
print('n tables =', len(tabs))
grid = lambda i: [[c.text.strip() for c in r.cells] for r in tabs[i].rows]
f = lambda x: (float(x) if x not in (None, '', 'NA') else np.nan)

def norm(g):
    g = (g or '').strip().lower()
    if 'non' in g and 'candidate' in g: return 'Non-candidate'
    if 'candidate' in g: return 'Candidate'
    if 't2dm' in g: return 'T2DM'
    if 'housekeeping' in g or 'hk' == g: return 'Housekeeping'
    return g

print('\n--- tab index sanity ---')
for i in range(len(tabs)):
    first = grid(i)[0][0][:48] if grid(i) else ''
    print(f'  [{i}] rows={len(tabs[i].rows)} first={first!r}')

S2 = grid(1); print('\nS2 header:', S2[0])
S1 = grid(0); print('S1 header:', S1[0])
S5 = grid(5); print('\nS5 header:', S5[0])
S17 = grid(18); print('\nS17 header:', S17[0])
S14 = grid(15); print('\nS14 header:', S14[0])
