# -*- coding: utf-8 -*-
"""Fig. 6 重绘（第二版）：直接由 AF1 官方表按来源分离取值"""
import os, shutil
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
OUT = P.OUT_MAIN

if not os.environ.get('FIG6_FROM_03'):
    raise SystemExit(
        '本脚本已被**整体取代**，因此不产生任何文件（唯一输出 Fig.6 已由 08 重新出具）。\n'
        '  正确做法：python 08_redraw_Fig6_labels_20260920.py   ← 官方执行顺序第 6 步，含 SI Table S18 逐格回归断言\n'
        '  被取代原因：本版 Fig.6 是 2026-09-17 之前的图形（ylim = 24、图例在绘图区内），与现行图件不一致。\n'
        '  如确需运行本旧版（仅作对照）：设 FIG6_FROM_03=1 后再运行。\n'
        '  （退出码 1 是设计如此，不是脚本坏了；详见本目录 README.md「已弃用」一节。）')

BK = os.path.join(OUT, '_backup_before_officialZ_redraw_20260917')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.linewidth': 0.7,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'ps.fonttype': 42})
C_GTEX, C_EQTL = '#C0392B', '#2471A3'
d = Document(AF)
tabs = [b for ch in d.element.body.iterchildren()
        for b in ([Table(ch, d)] if ch.tag == qn('w:tbl') else [])]
grid = lambda i: [[c.text.strip() for c in r.cells] for r in tabs[i].rows]
f = lambda x: (float(x) if x not in (None, '', 'NA') else np.nan)
def norm(g):
    g = (g or '').strip().lower()
    if 'non' in g and 'candidate' in g:
        return 'Non-candidate'
    if 'candidate' in g:
        return 'Candidate'
    if 't2dm' in g:
        return 'T2DM'
    if 'housekeeping' in g or 'hk' == g:
        return 'Housekeeping'
    return g

# GTEx: Table S2（测试集, P_ACAT_O）+ Table S1 分组；housekeeping 来自 Table S5
S2 = grid(1); i2 = {n: k for k, n in enumerate(S2[0])}
G = {r[0]: norm(r[1]) for r in grid(0)[1:] if r[0]}
GT = {}
for r in S2[1:]:
    if len(r) < 10 or not r[0]:
        continue
    p = f(r[i2['P_ACAT_O']])
    if not np.isnan(p):
        GT.setdefault((G.get(r[0], '?'), r[1]), []).append(p)
# housekeeping GTEx（Table S5: 每个基因 3 表型 ACAT-O 列）
S5 = grid(5); h5 = S5[0]
acat_cols = [i for i, h in enumerate(h5) if 'ACAT' in h.upper()]
for r in S5[1:]:
    if not r[0]:
        continue
    for ph, ci in zip(['DR', 'DN', 'DPN'], acat_cols[:3]) if len(acat_cols) >= 3 else []:
        p = f(r[ci])
        if not np.isnan(p):
            GT.setdefault(('Housekeeping', ph), []).append(p)
print('GTEx housekeeping ACAT-O 列数:', len(acat_cols))

# eQTLGen: Table S17（测试集）+ Table S14（housekeeping）
EQ = {}
for r in grid(18)[1:]:
    if len(r) < 8 or not r[0]:
        continue
    p = f(r[4])
    if not np.isnan(p):
        EQ.setdefault((norm(r[2]), r[1]), []).append(p)
for r in grid(15)[1:]:
    if len(r) < 6 or not r[0]:
        continue
    p = f(r[3])
    if not np.isnan(p):
        EQ.setdefault(('Housekeeping', r[1]), []).append(p)

GRPS = ['Candidate', 'Non-candidate', 'T2DM', 'Housekeeping']
PH = ['DR', 'DN', 'DPN']
def cp(k, n):
    lo = 0.0 if k == 0 else stats.beta.ppf(0.025, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(0.975, k + 1, n - k)
    return 100 * lo, 100 * hi

fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.0), sharey=True)
x = np.arange(len(GRPS)); w = 0.36
for pi, ph in enumerate(PH):
    ax = axes[pi]
    for src, DD, col, off in [('GTEx v8 multi-tissue ACAT-O', GT, C_GTEX, -w / 2),
                              ('eQTLGen whole blood', EQ, C_EQTL, w / 2)]:
        vals, los, his, ns = [], [], [], []
        for g in GRPS:
            ps = np.array(DD.get((g, ph), []))
            if len(ps) == 0:
                vals.append(0); los.append(0); his.append(0); ns.append((0, 0)); continue
            k = int((ps < 0.05).sum()); n = len(ps); r = 100 * k / n; lo, hi = cp(k, n)
            vals.append(r); los.append(r - lo); his.append(hi - r); ns.append((k, n))
        ax.bar(x + off, vals, w, color=col, edgecolor='black', lw=0.5, label=src if pi == 0 else None)
        ax.errorbar(x + off, vals, yerr=[los, his], fmt='none', ecolor='#222222', elinewidth=0.8, capsize=3)
        for j, v in enumerate(vals):
            ax.text(x[j] + off, v + his[j] + 0.5, '%d/%d' % ns[j], ha='center', fontsize=5.2)
        print(' %-4s %-28s' % (ph, src), dict(zip(GRPS, ['%d/%d' % t for t in ns])))
    ax.set_xticks(x); ax.set_xticklabels(['Cand', 'Non-cand', 'T2DM', 'HK'], fontsize=6.4)
    ax.set_ylim(0, 26)
    ax.set_title('(%s) %s' % ('abc'[pi], ph), fontsize=7.6, pad=4)
    if pi == 0:
        ax.set_ylabel('Genes with nominal p < 0.05 (%)')
        ax.legend(frameon=False, fontsize=6.0, loc='upper left')
fig.tight_layout()
for ext in ['png', 'pdf']:
    src = os.path.join(OUT, f'Fig6.{ext}')
    if os.path.exists(src) and not os.path.exists(os.path.join(BK, f'Fig6.{ext}')):
        shutil.copy2(src, os.path.join(BK, f'Fig6.{ext}'))
    fig.savefig(src, dpi=600)
plt.close(fig)
print('Fig6 redrawn (source-split, official Z)')
