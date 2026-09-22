# -*- coding: utf-8 -*-
"""P0 Step 6：从官方数据重绘 Figure 8（跨性状三基因集组织轴）"""
import os, csv, shutil
import numpy as np
from scipy import stats
from docx import Document
from docx.table import Table
from docx.oxml.ns import qn
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import paths_config as P  # 统一路径入口（2026-09-20）

AF = P.need(P.AF1, 'Additional file 1（从期刊补充材料下载后用 AF1_DOCX 指定）')
NEW = P.need(P.DATA_Z, '官方 MetaXcan Z 数据层')
OUT = P.OUT_MAIN
BK = os.path.join(OUT, '_backup_before_officialZ_redraw_20260917')
os.makedirs(BK, exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.linewidth': 0.7,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'ps.fonttype': 42})
C_TEST, C_HK, C_SCZ = '#C0392B', '#1E8449', '#2471A3'

d = Document(AF)
tabs = [b for ch in d.element.body.iterchildren()
        for b in ([Table(ch, d)] if ch.tag == qn('w:tbl') else [])]
grid = lambda i: [[c.text.strip() for c in r.cells] for r in tabs[i].rows]
f = lambda x: (float(x) if x not in (None, '', 'NA') else np.nan)

def rho_ci(x, y):
    m = ~(np.isnan(x) | np.isnan(y))
    a, b = x[m], y[m]
    r = stats.spearmanr(a, b)[0]
    n = len(a)
    z = np.arctanh(r); se = 1 / np.sqrt(n - 3)
    return r, n, np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)

def rdc(x, y):
    m = ~(np.isnan(x) | np.isnan(y))
    return 100 * ((x[m] > 0) == (y[m] > 0)).mean()

# --- 测试集：tissue-only（GTEx WB vs NT），全部有 Z 的基因 ---
S2 = grid(1); i2 = {v: k for k, v in enumerate(S2[0])}
G = {r[0]: r[1] for r in grid(0)[1:] if r[0] and 'andidate' in r[1] or (r[0] and r[1])}
zwb = np.array([f(r[i2['Z_Whole_Blood']]) for r in S2[1:] if len(r) >= 10 and r[0]])
znt = np.array([f(r[i2['Z_Nerve_Tibial']]) for r in S2[1:] if len(r) >= 10 and r[0]])
genes2 = [r[0] for r in S2[1:] if len(r) >= 10 and r[0]]
is_cand = np.array([('andidate' in str(G.get(g, '')) and 'Non' not in str(G.get(g, ''))) or
                    ('Non-Candidate' in str(G.get(g, ''))) or ('T2DM' in str(G.get(g, ''))) for g in genes2])
r_t, n_t, lo_t, hi_t = rho_ci(zwb[is_cand], znt[is_cand])
print('testbed tissue-only: rho=%.3f n=%d CI %.3f-%.3f  consistency=%.1f%%'
      % (r_t, n_t, lo_t, hi_t, rdc(zwb[is_cand], znt[is_cand])))

# --- housekeeping ---
S5 = grid(5); h5 = S5[0]
zi = [i for i, h in enumerate(h5) if 'Z, Nerve_Tibial' in h]
zj = [i for i, h in enumerate(h5) if 'Z, Whole_Blood' in h]
print('S5 cols NT:', zi, 'WB:', zj)
hk_nt, hk_wb = [], []
for r in S5[1:]:
    if len(r) <= max(zi + zj) or not r[0]:
        continue
    for a, b in zip(zi, zj):
        hk_nt.append(f(r[a])); hk_wb.append(f(r[b]))
hk_nt = np.array(hk_nt); hk_wb = np.array(hk_wb)
r_hk, n_hk, lo_hk, hi_hk = rho_ci(hk_wb, hk_nt)
print('housekeeping tissue-only: rho=%.3f n=%d CI %.3f-%.3f' % (r_hk, n_hk, lo_hk, hi_hk))

# --- SCZ genome-wide ---
rows = list(csv.DictReader(open(os.path.join(NEW, 'scz_z_4arm_official.csv'), encoding='utf-8-sig')))
g = lambda k: np.array([float(r[k]) if r[k] not in ('', 'NA') else np.nan for r in rows])
ze, zwb_s, znt_s = g('eqZ'), g('wbZ'), g('ntZ')
r_s, n_s, lo_s, hi_s = rho_ci(zwb_s, znt_s)
print('SCZ tissue-only: rho=%.3f n=%d CI %.3f-%.3f  consistency=%.1f%%'
      % (r_s, n_s, lo_s, hi_s, rdc(zwb_s, znt_s)))

panels = [('HOTAIR testbed\n(FinnGen R13, official Z)', r_t, n_t, lo_t, hi_t, C_TEST),
          ('Housekeeping control\n(FinnGen R13, official Z)', r_hk, n_hk, lo_hk, hi_hk, C_HK),
          ('PGC3 schizophrenia\n(genome-wide, official Z)', r_s, n_s, lo_s, hi_s, C_SCZ)]
fig, ax = plt.subplots(figsize=(5.2, 3.2))
xs = np.arange(3)
for i, (lab, r, n, lo, hi, c) in enumerate(panels):
    ax.errorbar(i, r, yerr=[[r - lo], [hi - r]], fmt='o', ms=6.5, color=c, ecolor=c, capsize=4, lw=1.2)
    ax.text(i, hi + 0.03, '$\\rho$ = %+.2f\nn = %d' % (r, n), ha='center', fontsize=6.2)
ax.set_xticks(xs); ax.set_xticklabels([p[0] for p in panels], fontsize=6.4)
ax.set_ylabel("Spearman $\\rho$ (GTEx Whole_Blood vs Nerve_Tibial)")
_hi = max(hi for _, _, _, _, hi, _ in panels)
ax.set_ylim(0.30, max(0.80, _hi + 0.11))   # 留白：避免 rho/n 标注压到标题（2026-09-20）
ax.axhline(0, color='#AAAAAA', lw=0.5)
ax.set_title('Within-GTEx tissue-context comparison across three gene sets', fontsize=7.4, pad=5)
fig.tight_layout()
for ext in ['png', 'pdf']:
    src = os.path.join(OUT, f'Fig8.{ext}')
    if os.path.exists(src) and not os.path.exists(os.path.join(BK, f'Fig8.{ext}')):
        shutil.copy2(src, os.path.join(BK, f'Fig8.{ext}'))
    fig.savefig(src, dpi=600)
plt.close(fig)
print('Fig8 redrawn')
print('backups:', sorted(os.listdir(BK)))
